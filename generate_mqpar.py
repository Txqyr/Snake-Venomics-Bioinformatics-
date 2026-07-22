#!/usr/bin/env python3
"""
generate_mqpar.py
=================
Generate MaxQuant mqpar.xml file(s) from either a CSV sample list or directly
from LAB_WORK_DAMIEN_VENOMICS.xlsx (no CSV needed).

Experiment names come from the Experiment_Name column in the lab workbook.
Genus/species are read from the workbook directly (or looked up from
Esquerre_Venom_Collection.xlsx as a fallback).

── Simplest usage (all samples from lab workbook, split by genus) ──────────
    python3 generate_mqpar.py \\
        --template        template_mqpar.xml \\
        --lab-database    LAB_WORK_DAMIEN_VENOMICS.xlsx \\
        --raw-list        raw_files.txt \\
        --fasta           /path/to/Serpentes.fasta \\
        --results-output  /home/shared/Venomics/MaxQuantAnalyses/GenusBatch \\
        --output-dir      xmls/ \\
        --split-by        genus \\
        --maxquant-exe    "dotnet /home/shared/software/MaxQuant/bin/MaxQuantCmd.dll"

    raw_files.txt is generated on the server with:
        find /data/raw -name "*.raw" > raw_files.txt
    then copied locally before running this script.

    With --split-by genus, each genus writes results to:
        {results-output}/{Genus}/   (e.g. .../GenusBatch/Pseudonaja/)

── Filter to a specific run or run type ────────────────────────────────────
    Add:  --filter-run  Jan2023          (or Nov2025, Feb2026, …)
          --filter-type DDA              (or DIA)

── CSV mode ────────────────────────────────────────────────────────────────
    python3 generate_mqpar.py \\
        --template        template_mqpar.xml \\
        --samples         sample_list.csv \\
        --database        Esquerre_Venom_Collection.xlsx \\
        --lab-database    LAB_WORK_DAMIEN_VENOMICS.xlsx \\
        --raw-list        raw_files.txt \\
        --fasta           /path/to/Serpentes.fasta \\
        --results-output  /home/shared/Venomics/MaxQuantAnalyses/MyRun \\
        --output-dir      xmls/

CSV columns (when --samples is used)
-------------------------------------
  raw_file      Required. .raw filename (must appear in --raw-list).
  devc_number   Required. DEVC identifier (e.g. DEVC_146 or DEVC 146).
  experiment    Optional. Overrides all auto-lookup logic.
  fraction      Optional. Integer fraction value (default 1).
"""

import argparse
import copy
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: pip install pandas openpyxl")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalize_devc(s: str) -> str:
    """Normalise any DEVC string to a canonical lookup key."""
    s = str(s).strip()
    m = re.match(r'DEVC\s*(\d+)', s, re.IGNORECASE)
    if m:
        return f"DEVC_{int(m.group(1)):03d}"
    return re.sub(r'\s+', '_', s).upper()


def sanitize(s: str) -> str:
    """Replace characters unsafe for filenames/XML text with underscores."""
    return re.sub(r'[^\w]', '_', s).strip('_')


# ---------------------------------------------------------------------------
# Raw file finder
# ---------------------------------------------------------------------------

def build_raw_file_index(root: str) -> dict[str, str]:
    """Walk *root* recursively and return a dict: lower-case filename -> full path.

    If the same filename appears in multiple subdirectories the first one found
    (depth-first) is kept and a warning is printed.
    """
    index: dict[str, str] = {}
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            key = fname.lower()
            if key not in index:
                index[key] = os.path.join(dirpath, fname)
            else:
                print(f"  [WARN] duplicate filename ignored: "
                      f"{os.path.join(dirpath, fname)}")
    return index


def build_raw_file_index_from_list(list_file: str) -> dict[str, str]:
    """Build a filename→full-path index from a pre-generated text file.

    The file should contain one absolute path per line, e.g. the output of:
        find /home/shared/Venomics/MassSpec_RawData -name "*.raw" > raw_files.txt

    This is equivalent to --raw-root but works when the raw files live on a
    remote server that the script cannot walk directly.  Generate the list
    once on the server, copy it alongside the script, and pass it with
    --raw-list.
    """
    index: dict[str, str] = {}
    with open(list_file, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            path = line.strip()
            if not path or path.startswith('#'):
                continue
            key = os.path.basename(path).lower()
            if key not in index:
                index[key] = path
            else:
                print(f"  [WARN] duplicate filename ignored: {path}")
    return index


# ---------------------------------------------------------------------------
# Database loader
# ---------------------------------------------------------------------------

def load_database(xlsx_path: str) -> dict:
    """Return dict: normalised_devc_key -> (Genus, species)."""
    df = pd.read_excel(xlsx_path, sheet_name='MASTER')
    db = {}
    for _, row in df.iterrows():
        raw_id  = str(row.get('D E Venom Collection Number', '')).strip()
        genus   = str(row.get('Genus',   '')).strip()
        species = str(row.get('Species', '')).strip()
        if not raw_id or raw_id == 'nan':
            continue
        if not genus or genus == 'nan' or not species or species == 'nan':
            continue
        db[normalize_devc(raw_id)] = (genus, species)
    return db


def load_experiment_names(xlsx_path: str) -> tuple[dict, dict]:
    """Load Experiment_Name column from LAB_WORK_DAMIEN_VENOMICS.xlsx MASTER sheet.

    Returns two dicts:
      by_raw_stem  — lower-case raw file stem (no extension) -> Experiment_Name
      by_devc      — normalised DEVC key -> Experiment_Name (first match wins)
    """
    df = pd.read_excel(xlsx_path, sheet_name='MASTER')
    by_raw_stem: dict[str, str] = {}
    by_devc:     dict[str, str] = {}

    for _, row in df.iterrows():
        raw_file = str(row.get('ORIGINAL RAW FILE NAME', '') or '').strip()
        devc_raw = str(row.get('D E Venom Collection Number', '') or '').strip()
        exp_name = str(row.get('Experiment_Name', '') or '').strip()

        if not exp_name or exp_name == 'nan':
            continue

        # Index by raw file stem (case-insensitive, no extension)
        if raw_file and raw_file != 'nan':
            stem = Path(raw_file).stem.lower()
            by_raw_stem[stem] = exp_name

        # Index by DEVC key (first occurrence wins — handles duplicate DEVCs)
        if devc_raw and devc_raw != 'nan':
            key = normalize_devc(devc_raw)
            if key not in by_devc:
                by_devc[key] = exp_name

    return by_raw_stem, by_devc


def load_samples_from_lab_database(
    xlsx_path: str,
    filter_run:  str = '',
    filter_type: str = '',
) -> list[dict]:
    """Read the MASTER sheet of LAB_WORK_DAMIEN_VENOMICS.xlsx as a sample list.

    Each row becomes a dict with keys: raw_file, devc_number, experiment,
    genus, species, fraction (always 1).

    Optional filters:
      filter_run   — keep only rows where Run == filter_run  (e.g. 'Jan2023')
      filter_type  — keep only rows where TypeOfRun == filter_type (e.g. 'DDA')
    """
    df = pd.read_excel(xlsx_path, sheet_name='MASTER')
    samples = []

    for _, row in df.iterrows():
        raw_file = str(row.get('ORIGINAL RAW FILE NAME', '') or '').strip()
        if not raw_file or raw_file == 'nan':
            continue

        # Optional filters
        if filter_run:
            row_run = str(row.get('Run', '') or '').strip()
            if row_run.lower() != filter_run.lower():
                continue
        if filter_type:
            row_type = str(row.get('TypeOfRun', '') or '').strip()
            if row_type.lower() != filter_type.lower():
                continue

        devc_raw = str(row.get('D E Venom Collection Number', '') or '').strip()
        exp_name = str(row.get('Experiment_Name', '') or '').strip()
        genus    = str(row.get('Genus',   '') or '').strip()
        species  = str(row.get('Species', '') or '').strip()

        run_time    = str(row.get('Run time',    '') or '').strip()
        mcwo_filter = str(row.get('MCWO Filter', '') or '').strip()
        run_val     = str(row.get('Run',          '') or '').strip()

        samples.append(dict(
            raw_file    = raw_file,
            devc_number = devc_raw     if devc_raw     != 'nan' else '',
            experiment  = exp_name     if exp_name     != 'nan' else '',
            genus       = genus        if genus        != 'nan' else 'Unknown',
            species     = species      if species      != 'nan' else 'unknown',
            fraction    = '',          # default 1 in resolve_rows
            # Pass-through fields used by --split-by
            run_time    = run_time     if run_time     != 'nan' else '',
            mcwo_filter = mcwo_filter  if mcwo_filter  != 'nan' else '',
            run         = run_val      if run_val      != 'nan' else '',
        ))

    return samples


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------

def _set_list_children(parent: ET.Element, tag: str, values: list):
    """Replace all children of *parent* with <tag>value</tag> elements.
    Preserves the element's tail (whitespace after closing tag) so that
    sibling elements stay on separate lines.
    """
    tail = parent.tail          # e.g. '\n   ' — whitespace before next sibling
    parent.clear()
    parent.tail = tail          # restore so </filePaths>\n   <experiments> is kept
    parent.text = '\n  '
    for i, v in enumerate(values):
        child = ET.SubElement(parent, tag)
        child.text = str(v)
        child.tail = '\n  ' if i < len(values) - 1 else '\n'


def _fix_xml_string(xml_str: str) -> str:
    """Post-process the serialised XML string to match MaxQuant's expected format.

    ElementTree introduces several deviations from the template:
      1. Single-quoted XML declaration  -> double-quoted
      2. Strips xmlns attributes from root element
      3. Escapes '>' as '&gt;' in text content (unnecessary but valid XML)
      4. Collapses empty elements to self-closing <tag /> form
    """
    # 1. Double-quote the XML declaration
    xml_str = xml_str.replace("<?xml version='1.0' encoding='utf-8'?>",
                               '<?xml version="1.0" encoding="utf-8"?>')

    # 2. Restore xmlns attributes on root element
    xml_str = xml_str.replace(
        '<MaxQuantParams>',
        '<MaxQuantParams xmlns:xsd="http://www.w3.org/2001/XMLSchema" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">',
        1,
    )

    # 3. Unescape &gt; back to > (safe in element text content)
    xml_str = xml_str.replace('&gt;', '>')

    # 4. Expand self-closing empty tags:  <tag />  ->  <tag></tag>
    #    Also handles <tag/> (no space)
    xml_str = re.sub(r'<(\w+)\s*/>', r'<\1></\1>', xml_str)

    return xml_str


# ---------------------------------------------------------------------------
# XML builder
# ---------------------------------------------------------------------------

def _update_sdrf_data(
    root: ET.Element,
    raw_files: list,
    experiments: list,
    fractions: list,
) -> None:
    """Ensure <sdrfData> exists and has one row per raw file.

    MaxQuant 2.6+ crashes with a NullReferenceException if <sdrfData> is
    absent (templates created with older versions don't have it).

    Strategy:
      • If the template already has <sdrfData>, keep its header/column
        structure but replace its <rows> so the count matches the new files.
      • If absent, inject a minimal block sufficient to prevent the crash.
    """
    _SDRF_HEADERS = [
        'source name',
        'technology type',
        'comment[label]',
        'comment[fraction identifier]',
        'comment[instrument]',
        'comment[data file]',
    ]

    sdrf = root.find('sdrfData')
    if sdrf is None:
        # Inject element directly before <parameterGroups> (or at end)
        sdrf = ET.Element('sdrfData')
        ref = root.find('parameterGroups')
        if ref is not None:
            ref.addprevious(sdrf) if hasattr(ref, 'addprevious') else root.append(sdrf)
        else:
            root.append(sdrf)

    # Ensure <headers> block exists
    headers_elem = sdrf.find('headers')
    if headers_elem is None:
        headers_elem = ET.SubElement(sdrf, 'headers')
        for h in _SDRF_HEADERS:
            s = ET.SubElement(headers_elem, 'string')
            s.text = h

    # Determine header list for row generation
    header_texts = [s.text or '' for s in headers_elem.findall('string')]
    if not header_texts:
        header_texts = _SDRF_HEADERS

    # Rebuild <rows> — one <List_string> per file
    rows_elem = sdrf.find('rows')
    if rows_elem is not None:
        sdrf.remove(rows_elem)
    rows_elem = ET.SubElement(sdrf, 'rows')

    for raw, exp, frac in zip(raw_files, experiments, fractions):
        raw_basename = os.path.basename(raw)
        row_el = ET.SubElement(rows_elem, 'List_string')
        for h in header_texts:
            s = ET.SubElement(row_el, 'string')
            hl = h.lower()
            if hl == 'source name':
                s.text = exp
            elif hl == 'technology type':
                s.text = 'proteomic profiling by mass spectrometry'
            elif 'label' in hl:
                s.text = 'label free sample'
            elif 'fraction' in hl:
                s.text = str(frac)
            elif 'data file' in hl:
                s.text = raw_basename
            else:
                s.text = ''


def build_mqpar(
    template_tree: ET.ElementTree,
    raw_files: list,
    experiments: list,
    fractions: list,
    fasta_path: str,
    threads: int,
    session_name: str,
    results_output: str,
) -> ET.ElementTree:
    """Return a modified deep copy of template_tree.

    results_output is written to <combinedFolder>, <fixedSearchFolder>, and
    <customTxtFolder> so MaxQuant writes its txt results there directly.
    """
    tree = copy.deepcopy(template_tree)
    root = tree.getroot()
    n = len(raw_files)

    # FASTA
    for elem in root.findall('.//fastaFilePath'):
        elem.text = fasta_path

    # Session name
    name_elem = root.find('name')
    if name_elem is not None:
        name_elem.text = session_name

    # Threads
    t_elem = root.find('numThreads')
    if t_elem is not None:
        t_elem.text = str(threads)

    # Output folders (root level)
    for tag in ('fixedSearchFolder', 'combinedFolder'):
        elem = root.find(tag)
        if elem is not None:
            elem.text = results_output

    # customTxtFolder — directs the txt results to the results_output directory
    txt_elem = root.find('customTxtFolder')
    if txt_elem is None:
        txt_elem = ET.SubElement(root, 'customTxtFolder')
    txt_elem.text = results_output

    # Output folders inside parameterGroup
    for pg in root.findall('.//parameterGroup'):
        cf = pg.find('combinedFolder')
        if cf is not None:
            cf.text = results_output
        tf = pg.find('tempFolder')
        if tf is not None:
            tf.text = os.path.join(results_output, 'tmp')

    # Single-sample safety: MaxQuant hangs (no log output, never terminates)
    # if matchBetweenRuns or LFQ are on with only one raw file because the
    # alignment / LFQ-normalisation steps require ≥2 samples.
    # iBAQ is unaffected and stays on.
    if n == 1:
        mbr_elem = root.find('matchBetweenRuns')
        if mbr_elem is not None:
            mbr_elem.text = 'False'
        for pg in root.findall('.//parameterGroup'):
            lfq = pg.find('lfqMode')
            if lfq is not None:
                lfq.text = '0'
            fast_lfq = pg.find('fastLfq')
            if fast_lfq is not None:
                fast_lfq.text = 'False'
        print(f"      [single-sample] {session_name}: disabled "
              f"matchBetweenRuns / lfqMode / fastLfq.")

    # File paths
    fp_elem = root.find('filePaths')
    if fp_elem is not None:
        _set_list_children(fp_elem, 'string', raw_files)

    # Experiment labels
    exp_elem = root.find('experiments')
    if exp_elem is not None:
        _set_list_children(exp_elem, 'string', experiments)

    # Fractions (per-file values passed in)
    frac_elem = root.find('fractions')
    if frac_elem is not None:
        _set_list_children(frac_elem, 'short', [str(f) for f in fractions])

    # Per-file repeated elements
    ptm_elem = root.find('ptms')
    if ptm_elem is not None:
        _set_list_children(ptm_elem, 'boolean', ['False'] * n)

    pgi_elem = root.find('paramGroupIndices')
    if pgi_elem is not None:
        _set_list_children(pgi_elem, 'int', ['0'] * n)

    rc_elem = root.find('referenceChannel')
    if rc_elem is not None:
        _set_list_children(rc_elem, 'string', [''] * n)

    # MaxQuant 2.6+ requires a <sdrfData> block with one row per file.
    # Rebuild it so templates saved with older MQ versions don't crash.
    _update_sdrf_data(root, raw_files, experiments, fractions)

    return tree


def write_xml(tree: ET.ElementTree, path: str):
    """Serialise *tree* to *path* with MaxQuant-compatible formatting."""
    xml_str = ET.tostring(tree.getroot(), encoding='unicode', xml_declaration=False)
    xml_str = '<?xml version=\'1.0\' encoding=\'utf-8\'?>\n' + xml_str
    xml_str = _fix_xml_string(xml_str)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(xml_str)


# ---------------------------------------------------------------------------
# Row resolver
# ---------------------------------------------------------------------------

def resolve_rows(csv_rows, db, exp_by_stem, exp_by_devc,
                 raw_dir: str = '', raw_index: dict = {}):
    """Parse and resolve every row into a dict with all fields filled in.

    Accepts rows from either a CSV file or load_samples_from_lab_database().
    When genus/species are already present in the row (lab-database mode) the
    Esquerre_Venom_Collection lookup is skipped.

    raw_index: dict built by build_raw_file_index() for recursive path lookup.
               Takes precedence over raw_dir when present.

    Returns (resolved_list, skipped_list).
    Each resolved entry is a dict with keys:
      raw_path, exp_label, fraction, genus, species, devc_key, devc_raw
    """
    resolved = []
    skipped  = []

    def _clean_path(p: str) -> str:
        """Strip shell-quoting artifacts from a path.

        Handles paths that were accidentally stored with surrounding quotes
        (single or double) or bash backslash-escapes (e.g. 'Jan\\ 2023').
        """
        p = p.strip()
        # Strip surrounding quotes (single or double)
        if len(p) >= 2 and p[0] == p[-1] and p[0] in ('"', "'"):
            p = p[1:-1]
        # Unescape bash backslash sequences: '\ ' → ' ', '\(' → '(', etc.
        import re as _re
        p = _re.sub(r'\\(.)', r'\1', p)
        return p

    for i, row in enumerate(csv_rows, 1):
        raw_file  = row.get('raw_file',   '').strip()
        devc_raw  = row.get('devc_number','').strip()
        exp_label = row.get('experiment', '').strip()
        frac_str  = str(row.get('fraction', '') or '').strip()
        # Pre-resolved genus/species (lab-database mode; blank in CSV mode)
        genus_pre   = row.get('genus',   '')
        species_pre = row.get('species', '')

        if not raw_file:
            print(f"  Row {i}: [SKIP] empty raw_file.")
            skipped.append(f"Row {i}: empty raw_file")
            continue

        # Strip any shell-quoting / backslash-escaping from the stored filename
        raw_file = _clean_path(raw_file)

        # Resolve path (priority: raw_index > raw_dir > as-is)
        if raw_index:
            found = raw_index.get(raw_file.lower())
            if found:
                raw_path = found
            else:
                print(f"  Row {i}: [WARN] '{raw_file}' not found in raw file index (--raw-root/--raw-list).")
                raw_path = raw_file   # keep original; MaxQuant will error if missing
        elif raw_dir and not os.path.isabs(raw_file):
            raw_path = os.path.join(raw_dir, raw_file)
        else:
            raw_path = _clean_path(raw_file)

        fraction_val = int(frac_str) if frac_str else 1

        # Genus/species: use pre-resolved values from lab sheet, else look up
        if genus_pre and genus_pre not in ('', 'Unknown'):
            genus, species = genus_pre, species_pre
            devc_key = normalize_devc(devc_raw) if devc_raw else ''
        elif devc_raw:
            devc_key = normalize_devc(devc_raw)
            lookup   = db.get(devc_key) if db else None
            if not lookup:
                print(f"  Row {i}: [WARN] '{devc_raw}' not found in database.")
                genus, species = 'Unknown', 'unknown'
            else:
                genus, species = lookup
        else:
            devc_key = ''
            genus, species = 'Unknown', 'unknown'

        # Resolve experiment label (priority order):
        #   1. Explicit value in row (CSV column or lab-sheet Experiment_Name)
        #   2. Lab database match by raw file stem
        #   3. Lab database match by DEVC number
        #   4. Auto-generate from DEVC + genus + species + stem
        stem = Path(raw_file).stem
        source = 'sheet' if exp_label else ''
        if not exp_label and exp_by_stem:
            exp_label = exp_by_stem.get(stem.lower(), '')
            if exp_label:
                source = 'lab-db/filename'
        if not exp_label and exp_by_devc and devc_raw:
            exp_label = exp_by_devc.get(devc_key, '')
            if exp_label:
                source = 'lab-db/devc'
        if not exp_label:
            devc_prefix = sanitize(devc_raw.upper().replace(' ', '_')) if devc_raw else ''
            suffix_part = (stem[len(devc_prefix):].lstrip('_')
                           if devc_prefix and stem.upper().startswith(devc_prefix.upper())
                           else stem)
            parts = []
            if devc_raw:
                parts.append(sanitize(devc_raw.upper().replace(' ', '_')))
            parts += [sanitize(genus), sanitize(species)]
            if suffix_part:
                parts.append(sanitize(suffix_part))
            exp_label = '_'.join(parts)
            source = 'auto'

        frac_label = "independent" if fraction_val == 32767 else f"fraction {fraction_val}"
        print(f"  Row {i}: {raw_path}")
        print(f"           exp [{source:16s}]: {exp_label}")
        print(f"           fraction         : {frac_label}")

        resolved.append(dict(
            raw_path=raw_path,
            exp_label=exp_label,
            fraction=fraction_val,
            genus=genus,
            species=species,
            devc_key=devc_key,
            devc_raw=devc_raw,
            # Pass-through metadata for compound split keys
            run_time=str(row.get('Run time') or row.get('run_time') or '').strip(),
            mcwo_filter=str(row.get('MCWO Filter') or row.get('mcwo_filter') or '').strip(),
            run=str(row.get('Run') or row.get('run') or '').strip(),
        ))

    return resolved, skipped


# ---------------------------------------------------------------------------
# Batch runner script generator
# ---------------------------------------------------------------------------

def write_run_script(xml_paths: list[Path], maxquant_exe: str, output_dir: Path):
    """Write a bash shell script that runs MaxQuant on each XML sequentially.

    All paths in the script are relative to the script's own directory so the
    entire output folder can be copied to a remote server and run as-is.
    """
    script_path = output_dir / 'run_all.sh'
    lines = [
        '#!/usr/bin/env bash',
        '# Auto-generated by generate_mqpar.py',
        '# Runs MaxQuant on each mqpar.xml in sequence.',
        '# Copy this folder (all *_mqpar.xml files + run_all.sh) to the server,',
        '# then run:  bash run_all.sh',
        '',
        '# Resolve the directory this script lives in (works when copied anywhere)',
        'SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"',
        # Use a bash array so "dotnet /path/to/MaxQuantCmd.dll" splits correctly
        # even when the exe string contains spaces (dotnet + DLL path).
        f'MAXQUANT=({maxquant_exe})',
        'LOG_DIR="$SCRIPT_DIR/logs"',
        'mkdir -p "$LOG_DIR"',
        '',
        'overall_start=$(date +%s)',
        'failures=0',
        '',
    ]
    for xml in xml_paths:
        name = xml.stem   # e.g. Acanthophis_mqpar
        # Use just the filename — the XML sits next to run_all.sh
        xml_filename = xml.name
        lines += [
            f'echo "================================================================"',
            f'echo "[$(date +\'%Y-%m-%d %H:%M:%S\')] Starting {name}"',
            f'echo "================================================================"',
            f'"${{MAXQUANT[@]}}" "$SCRIPT_DIR/{xml_filename}" '
            f'2>&1 | tee "$LOG_DIR/{name}.log"',
            f'if [ ${{PIPESTATUS[0]}} -ne 0 ]; then',
            f'  echo "[WARN] {name} exited with an error — continuing to next run."',
            f'  failures=$((failures + 1))',
            f'fi',
            '',
        ]
    lines += [
        'overall_end=$(date +%s)',
        'elapsed=$(( overall_end - overall_start ))',
        'echo "================================================================"',
        'echo "All runs complete in $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"',
        'echo "Failed runs: $failures / ' + str(len(xml_paths)) + '"',
        'echo "Logs written to: $LOG_DIR"',
        'echo "================================================================"',
    ]
    script_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    script_path.chmod(0o755)
    return script_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate MaxQuant mqpar.xml file(s) from a CSV sample list.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--template',      required=True,
                        help='Template mqpar.xml')
    parser.add_argument('--samples',       default='',
                        help='CSV file (one row per .raw file). '
                             'If omitted, --lab-database is used as the sample source.')
    parser.add_argument('--database',      default='',
                        help='Esquerre_Venom_Collection.xlsx (for genus/species lookup). '
                             'Not needed when --lab-database is the sample source.')
    parser.add_argument('--lab-database',  default='',
                        help='LAB_WORK_DAMIEN_VENOMICS.xlsx. '
                             'When --samples is omitted this is used as the sample source.')
    parser.add_argument('--filter-run',    default='',
                        help='Only include rows where Run == this value '
                             '(e.g. Jan2023, Nov2025, Feb2026). Lab-database mode only.')
    parser.add_argument('--filter-type',   default='',
                        help='Only include rows where TypeOfRun == this value '
                             '(e.g. DDA, DIA). Lab-database mode only.')
    parser.add_argument('--raw-list',       default='',
                        help='Text file with one absolute .raw/.mzML path per line. '
                             'Generate on the server with: find /data -name "*.raw" > raw_files.txt '
                             'then copy locally and pass here. '
                             'Use --raw-root instead if running the script on the same machine as the data.')
    parser.add_argument('--raw-root',       default='',
                        help='Root directory to search recursively for .raw files by name. '
                             'Only useful when the raw files are on the same machine as this script. '
                             'Prefer --raw-list for remote/server workflows.')
    parser.add_argument('--fasta',         default='',
                        help='Path to FASTA file')
    parser.add_argument('--output-dir',    default='.',
                        help='Where to write the output XML(s)')
    parser.add_argument('--session-name',  default='',
                        help='Session name / XML stem (single-file mode only)')
    parser.add_argument('--results-output', default='',
                        help='Base directory where MaxQuant writes results. '
                             'With --split-by, each genus gets its own sub-folder here '
                             '(e.g. /home/shared/Venomics/MaxQuantAnalyses/GenusBatch). '
                             'Sets <combinedFolder>, <fixedSearchFolder>, and <customTxtFolder> in the XML.')
    parser.add_argument('--threads',       type=int, default=None,
                        help='CPU threads')
    parser.add_argument('--split-by',      nargs='+',
                        choices=['genus', 'species', 'run_time', 'mcwo_filter', 'run'],
                        default=None,
                        help='Generate one XML per unique combination of these fields and write run_all.sh. '
                             'Multiple values allowed, e.g.: --split-by genus run_time mcwo_filter. '
                             'Choices: genus, species, run_time, mcwo_filter, run.')
    parser.add_argument('--maxquant-exe',  default='',
                        help='MaxQuant executable for run_all.sh '
                             '(e.g. "dotnet /path/to/MaxQuantCmd.dll"). '
                             'Required when --split-by is used.')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # -- Validate arguments --
    if not args.samples and not args.lab_database:
        sys.exit("Either --samples or --lab-database must be provided.")
    if args.samples and not args.lab_database and not args.database:
        sys.exit("--database (Esquerre_Venom_Collection.xlsx) is required in CSV mode.")

    # -- Load databases --
    print("[1/3] Loading databases ...")
    db: dict = {}
    if args.database:
        db = load_database(args.database)
        print(f"      Esquerre collection: {len(db)} entries loaded.")

    exp_by_stem: dict[str, str] = {}
    exp_by_devc: dict[str, str] = {}
    if args.lab_database:
        exp_by_stem, exp_by_devc = load_experiment_names(args.lab_database)
        print(f"      Lab workbook: {len(exp_by_stem)} raw-file / "
              f"{len(exp_by_devc)} DEVC experiment-name entries.")
    print()

    # -- Parse template --
    print(f"[2/3] Parsing template: {args.template}")
    ET.register_namespace('xsd', 'http://www.w3.org/2001/XMLSchema')
    ET.register_namespace('xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    template_tree = ET.parse(args.template)
    root_elem = template_tree.getroot()
    template_fasta    = root_elem.findtext('.//fastaFilePath') or ''
    template_threads  = int(root_elem.findtext('numThreads') or 16)
    template_combined = (root_elem.findtext('combinedFolder') or
                         root_elem.findtext('.//combinedFolder') or '')
    print(f"      FASTA  : {template_fasta}")
    print(f"      Threads: {template_threads}\n")

    resolved_fasta   = args.fasta or template_fasta
    resolved_threads = args.threads or template_threads
    base_combined    = args.results_output or template_combined

    # -- Load sample rows --
    if args.samples:
        # CSV mode
        with open(args.samples, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            input_rows = list(reader)
        print(f"[3/3] Processing {len(input_rows)} row(s) from CSV: {args.samples}\n")
    else:
        # Lab-database mode — read MASTER sheet directly
        filters = []
        if args.filter_run:
            filters.append(f"Run={args.filter_run}")
        if args.filter_type:
            filters.append(f"TypeOfRun={args.filter_type}")
        filter_desc = f" (filtered: {', '.join(filters)})" if filters else " (all samples)"
        print(f"[3/3] Reading samples from lab workbook{filter_desc} ...")
        input_rows = load_samples_from_lab_database(
            args.lab_database,
            filter_run=args.filter_run,
            filter_type=args.filter_type,
        )
        print(f"      {len(input_rows)} samples loaded.\n")

    # Build raw-file index: --raw-list (pre-generated) or --raw-root (live walk)
    raw_index: dict[str, str] = {}
    if args.raw_list:
        print(f"      Loading raw file list from: {args.raw_list} ...")
        raw_index = build_raw_file_index_from_list(args.raw_list)
        print(f"      {len(raw_index)} files indexed.\n")
    elif args.raw_root:
        print(f"      Indexing raw files under: {args.raw_root} ...")
        raw_index = build_raw_file_index(args.raw_root)
        print(f"      {len(raw_index)} files indexed.\n")

    resolved, skipped = resolve_rows(
        input_rows, db, exp_by_stem, exp_by_devc,
        raw_dir='', raw_index=raw_index)

    if not resolved:
        sys.exit("No valid files found — nothing to write.")

    print()

    # ── Build groups ──────────────────────────────────────────────────────────
    # Map --split-by token → resolved row key
    SPLIT_KEY_MAP = {
        'genus':       'genus',
        'species':     'species',
        'run_time':    'run_time',
        'mcwo_filter': 'mcwo_filter',
        'run':         'run',
    }
    if args.split_by:
        groups: dict[str, list] = {}
        for r in resolved:
            parts = []
            for k in args.split_by:
                val = sanitize(r.get(SPLIT_KEY_MAP[k], '') or 'unknown')
                # When splitting by species, always prefix with genus so filenames
                # are unambiguous (e.g. Pseudonaja_aspidorhyncha, not just aspidorhyncha).
                # Skip the prefix only if genus is already an explicit split key
                # (it would already appear as its own part).
                if k == 'species' and 'genus' not in args.split_by:
                    genus_val = sanitize(r.get('genus', '') or 'unknown')
                    if genus_val:
                        val = f"{genus_val}_{val}"
                parts.append(val)
            key = '_'.join(p for p in parts if p)
            groups.setdefault(key, []).append(r)
        print(f"  Split by {' + '.join(args.split_by)}: {len(groups)} group(s) — "
              f"{', '.join(sorted(groups))}\n")
    else:
        # Single XML: one group with the full set
        if args.session_name:
            group_key = sanitize(args.session_name)
        else:
            devc_labels = [sanitize(r['devc_raw'].upper().replace(' ', '_'))
                           for r in resolved if r['devc_raw']]
            unique = list(dict.fromkeys(devc_labels))
            group_key = (unique[0] if len(unique) == 1
                         else f"{unique[0]}_to_{unique[-1]}"
                         if unique else 'maxquant_batch')
        groups = {group_key: resolved}

    # ── Generate one XML per group ────────────────────────────────────────────
    generated_xmls: list[Path] = []

    for group_key, rows in sorted(groups.items()):
        session_name = group_key

        # Each group gets its own sub-folder under results_output so runs don't
        # overwrite each other's results.
        if args.split_by and base_combined:
            results_output = os.path.join(base_combined, session_name)
        else:
            results_output = base_combined

        xml_path = output_dir / f"{session_name}_mqpar.xml"

        modified = build_mqpar(
            template_tree=template_tree,
            raw_files=[r['raw_path']  for r in rows],
            experiments=[r['exp_label'] for r in rows],
            fractions=[r['fraction']  for r in rows],
            fasta_path=resolved_fasta,
            threads=resolved_threads,
            session_name=session_name,
            results_output=results_output,
        )
        write_xml(modified, str(xml_path))
        generated_xmls.append(xml_path)
        print(f"  Written: {xml_path.name}  ({len(rows)} file(s))")

    # ── run_all.sh ────────────────────────────────────────────────────────────
    run_script = None
    if args.split_by and len(generated_xmls) > 1:
        exe = args.maxquant_exe or 'MaxQuantCmd.exe'
        run_script = write_run_script(generated_xmls, exe, output_dir)

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print(f"Mode         : {'split by ' + ' + '.join(args.split_by) if args.split_by else 'single XML'}")
    print(f"XMLs written : {len(generated_xmls)}")
    print(f"Output dir   : {output_dir.resolve()}")
    if run_script:
        print(f"Run script   : {run_script.resolve()}")
        print(f"  → bash {run_script.resolve()}")
    if skipped:
        print(f"Skipped rows : {len(skipped)}")
        for s in skipped:
            print(f"  - {s}")
    print("=" * 60)


if __name__ == '__main__':
    main()
