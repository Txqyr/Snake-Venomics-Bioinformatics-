#!/usr/bin/env python3
"""
venomics_pipeline.py
====================
One-shot venomics LC-MS/MS processing pipeline.

Stage 0 first removes flagged proteins (Reverse decoys, Potential
contaminants, Only identified by site) from every MaxQuant
proteinGroups*.txt under --input-dir, writing cleaned copies into
<out>/FilteredProteinGroups/ (or --filtered-dir) with the same folder
structure and filenames as the input. The pipeline then walks THIS
filtered directory for the remaining stages, so the masterfile, summaries,
and pie charts are all built from the cleaned data. Use --skip-filter to
walk --input-dir unfiltered instead.

After Stage 0, the pipeline walks the filtered directory for MaxQuant
proteinGroups.txt files (in any subdirectory), annotates each against a
master toxin classification table, merges them into a single masterfile,
and produces:

  * <out>/FilteredProteinGroups/<species>/proteinGroups.txt  Stage 0 cleaned MaxQuant runs
  * <out>/AllProteinGroups_anotated.csv            merged annotated catalogue
  * <out>/protein_groups_masterfile.csv            all iBAQ columns + annotations
  * <out>/MassSpecData.csv                         transposed format (like R script expects)
  * <out>/SummaryResults.csv                          per-sample family % + Venom_pct column
  * <out>/SummaryResults_AveragedBySpeciesAndSampleType.csv  species/form means
  * <out>/PieCharts/*.pdf                                all-proteins pie charts
  * <out>/PieCharts_filteredOver1percent/*.pdf           all-proteins ≥1% (--filtered-pies)
  * <out>/PieCharts_VenomOnly/*.pdf                      venom-only pie charts (re-normalised)
  * <out>/PieCharts_VenomOnly_filteredOver1percent/*.pdf venom-only ≥1% (--filtered-pies)
  * <out>/PieCharts_VenomOnly_AveragePerSpecies/*.pdf    (skipped with --skip-species-pies)
  * <out>/unclassified_for_review.csv                    proteins not classified by any method

Pie chart colours come verbatim from AnalyseCleanData_2.R.

Every stage has an on/off toggle so you can re-run only the parts you need.

ANNOTATION STRATEGY (three tiers, applied in order):
  Tier 0  Exact string match against --master-annotation CSV
  Tier 1  Hard-coded venom keyword rules (applied to full FASTA header text)
  Tier 2  Frequency-weighted scoring from --unique-words CSV
  Tier 0 overrides all; Tier 1 before Tier 2 for unmatched rows.

USAGE
-----
Full pipeline (with auto-classification of new proteins):
    python3 venomics_pipeline.py \\
        --input-dir         /path/to/MaxQuant_runs_root \\
        --master-annotation AllProteinGroups_anotated3.csv \\
        --unique-words      unique_words__Family.csv \\
        --collection        Esquerre_Venom_Collection.xlsx \\
        --output-dir        results/

This first writes filtered MaxQuant runs to results/FilteredProteinGroups/
and then builds the masterfile, summaries, and pie charts from those
filtered files.

Skip Stage 0 and use the raw, unfiltered MaxQuant output instead:
    python3 venomics_pipeline.py \\
        --input-dir         /path/to/MaxQuant_runs_root \\
        --master-annotation AllProteinGroups_anotated3.csv \\
        --collection        Esquerre_Venom_Collection.xlsx \\
        --output-dir        results/ \\
        --skip-filter

Regenerate plots from an existing masterfile only:
    python3 venomics_pipeline.py \\
        --masterfile results/protein_groups_masterfile.csv \\
        --collection Esquerre_Venom_Collection.xlsx \\
        --output-dir results/
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Optional

try:
    import pandas as pd
    import numpy as np
except ImportError:
    sys.exit("pandas and numpy required: pip install pandas numpy openpyxl")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
except ImportError:
    sys.exit("matplotlib required: pip install matplotlib")


# ═════════════════════════════════════════════════════════════════════════════
# SECTION A — Classification constants (Tiers 1 & 2)
# ═════════════════════════════════════════════════════════════════════════════

# Colour palette verbatim from AnalyseCleanData_2.R lines 127-153.
COLOR_MAPPING: dict[str, str] = {
    '3FTx':                '#F47B5B',
    'PLA2':                '#B5EFB5',
    'SVMP':                '#FBE426',
    'SVSP':                '#2ED9FF',
    'CTL':                 '#AA0DFE',
    'CRiSP':               '#AAF400',
    'Kun':                 '#3283FE',
    'LAAO':                '#F7E1A0',
    'Prothrombin activator': '#E580DB',
    '5N':                  '#F8A19F',
    'AChE':                '#F6222E',
    'AmPep':               '#DD6F91',
    'Cys':                 '#F45366',
    'Dis':                 '#1C7F93',
    'Hyal':                '#1C8356',
    'NGF':                 '#B10DA1',
    'NP':                  '#1CBE4F',
    'Oha-Vesp':            '#7ED7D1',
    'PDE':                 '#C075A6',
    'PLB':                 '#FC1CBF',
    'PLC':                 '#FA0087',
    'VEGF':                '#BDCDFF',
    'Venom factor':        '#822E1C',
    'Venom peroxiredoxin': '#C4451C',
    'Waprin':              '#782AB6',
    'Kazal':               '#9FC2CC',
    'Other':               '#A2A2A5',
    'Unidentified':        '#D6D6D6',
}

FAMILY_ORDER: list[str] = [
    '3FTx', 'PLA2', 'SVMP', 'SVSP', 'CTL', 'CRiSP', 'Kun', 'LAAO',
    'Prothrombin activator',
    '5N', 'AChE', 'AmPep', 'Cys', 'Dis', 'Hyal', 'NGF', 'NP',
    'Oha-Vesp', 'PDE', 'PLB', 'PLC', 'VEGF',
    'Venom factor', 'Venom peroxiredoxin', 'Waprin',
    'Other',
    'Kazal',
]

COLLECTION_META_COLS: list[str] = [
    'D E Venom Collection Number',
    'Family',
    'Genus',
    'Species',
    'Species2',
    'Form_short',
    'Sex',
    'Date collected',
    'SVL',
    'Tail length',
]

ANNOT_COLS: list[str] = ['Protein', 'Family', 'Group', 'Subgroup',
                         'Function', 'Notes', 'UniProt_ID']

# ── Tier 1: hard-coded venom keyword rules ────────────────────────────────────
# Each tuple: (Family, Function, [keyword_list])
# Matched case-insensitively against the full FASTA header text.
# Rules are evaluated in order; first match wins.

VENOM_RULES: list[tuple[str, str, list[str]]] = [
    ('3FTx',  'Toxin', ['alpha-bungarotoxin', 'alpha bungarotoxin',
                         'alpha-cobratoxin', 'alpha cobratoxin',
                         'alpha-elapitoxin', 'alpha elapitoxin',
                         'muscarinic toxin', 'cardiotoxin', 'cytotoxin',
                         'fasciculin', 'neurotoxin homolog', 'dendrotoxin',
                         'weak neurotoxin', 'short neurotoxin',
                         'long neurotoxin', 'three-finger toxin',
                         'threefinger toxin', '3ftx',
                         'cobrotoxin', 'erabutoxin',
                         'alpha-neurotoxin', 'beta-neurotoxin']),
    ('PLA2',  'Toxin', ['phospholipase a2', 'phospholipase-a2',
                         'phospholipase a 2', 'pla2', 'pla-2',
                         'phospholipase a(2)',
                         'notexin', 'taipoxin', 'paradoxin',
                         'venom phospholipase',
                         'group ii phospholipase', 'group ia phospholipase',
                         'pseudarin', 'pseudexin',
                         'agkistrodotoxin', 'crotoxin',
                         'beta-bungarotoxin']),
    ('SVMP',  'Toxin', ['snake venom metalloproteinase',
                         'svmp', 'adamalysin', 'reprolysin',
                         'zinc metalloproteinase', 'zinc metalloprotease',
                         'metalloproteinase-disintegrin',
                         'bothropasin', 'jararhagin', 'atrolysin',
                         'lebetase', 'svmp-hop', 'svmp-aca', 'svmp-pse',
                         'peptidase m12b',
                         'metalloproteinase type iii',
                         'metalloproteinase type i',
                         'metalloproteinase type ii']),
    ('SVSP',  'Toxin', ['snake venom serine protease',
                         'thrombin-like enzyme', 'thrombin-like protease',
                         'ancrod', 'batroxobin', 'cerastobin',
                         'venombin', 'reptilase', 'brevinase',
                         'venom serine proteinase', 'serine protease svsp',
                         'stephensease', 'scutellarase',
                         'pseudarinase', 'textarinase']),
    ('CTL',   'Toxin', ['c-type lectin-like', 'c-type lectin like',
                         'snaclec', 'botrocetin', 'vipecetin',
                         'convulxin', 'factor ix/x-binding',
                         'echicetin', 'agkisacucetin',
                         'snake venom lectin']),
    ('CRiSP', 'Toxin', ['cysteine-rich secretory protein',
                         'cysteine rich secretory protein',
                         'crisp', 'natrin', 'pseudechetoxin',
                         'helothermine', 'ablomin', 'triflin',
                         'latisemin', 'ophanin']),
    ('Kun',   'Toxin', ['kunitz-type serine protease inhibitor',
                         'kunitz type serine protease inhibitor',
                         'kunitz-type protease inhibitor',
                         'bpti-kunitz', 'bpti/kunitz',
                         'textilinin', 'calcicludine', 'venom kunitz',
                         'kp-aca', 'kp-hem', 'kp-pse', 'kp-den',
                         'kp-tri', 'kp-oxy', 'kp-aus', 'kp-not',
                         'kp-sus', 'kp-hop']),
    ('LAAO',  'Toxin', ['l-amino acid oxidase', 'l-amino-acid oxidase',
                         'l amino acid oxidase', 'laao',
                         'venom l-amino']),
    ('Hyal',  'Toxin', ['hyaluronidase', 'hyaluronate glycanohydrolase',
                         'venom hyaluronidase']),
    ('5N',    'Toxin', ['snake venom 5-nucleotidase',
                         'snake venom 5 nucleotidase',
                         'venom 5-nucleotidase',
                         'ecto-5-nucleotidase', 'ecto 5-nucleotidase',
                         '5-nucleotidase ecto']),
    ('PDE',   'Toxin', ['venom phosphodiesterase',
                         'snake venom phosphodiesterase',
                         'phosphodiesterase i']),
    ('PLB',   'Toxin', ['phospholipase b', 'phospholipase-b',
                         'lysophospholipase', 'venom plb']),
    ('PLC',   'Toxin', ['phospholipase c', 'phospholipase-c']),
    ('NGF',   'Toxin', ['nerve growth factor', 'venom ngf',
                         'snake venom nerve growth']),
    ('NP',    'Toxin', ['natriuretic peptide', 'natriuretic protein',
                         'cnp-like', 'dnp-like']),
    ('VEGF',  'Toxin', ['vascular endothelial growth factor',
                         'vegf-like', 'venom vegf',
                         'snake venom vascular']),
    ('AChE',  'Toxin', ['acetylcholinesterase',
                         'snake venom cholinesterase']),
    ('AmPep', 'Toxin', ['venom aminopeptidase',
                         'snake venom aminopeptidase']),
    ('Waprin','Toxin', ['waprin', 'omwaprin', 'lachewaprin',
                         'supwaprin', 'nawaprin',
                         'ku-wap-fusin', 'kuwapfusin', 'wap-fusin']),
    ('Cys',   'Toxin', ['cystatin', 'cystatin b', 'venom cystatin',
                         'snake venom cystatin']),
    ('Dis',   'Toxin', ['disintegrin', 'bitistatin', 'trigramin',
                         'echistatin', 'kistrin', 'rhodostomin',
                         'obtustatin', 'vipera disintegrin']),
    ('Oha-Vesp','Toxin',['vespryn', 'ohanin', 'ohanin-like']),
    ('Venom factor',        'Toxin', ['cobra venom factor', 'cvf-like',
                                       'venom complement-depleting']),
    ('Venom peroxiredoxin', 'Toxin', ['venom peroxiredoxin',
                                       'snake venom peroxiredoxin']),
    ('Prothrombin activator', 'Toxin', ['prothrombin activator',
                              'pseutarin', 'oscutarin',
                              'trocarin', 'notecarin',
                              'coagulation factor xa',
                              'coagulation factor x-activating',
                              'coagulation factor v-activating']),
    ('Kazal',  'Toxin', ['kazal-type serine protease inhibitor',
                          'kazal type serine protease inhibitor',
                          'kazal-type inhibitor', 'kazal type inhibitor',
                          'kazal domain', 'kazal protease inhibitor']),
    # Broader fallbacks (less specific — placed after narrow rules)
    ('SVSP',  'Toxin', ['snake venom serine', 'venom serine protease']),
    ('SVMP',  'Toxin', ['metalloproteinase']),
    ('Kun',   'Toxin', ['kunitz', 'kunitz-type', '/bpti']),
    ('5N',    'Toxin', ['5-nucleotidase', '5 nucleotidase', 'nucleotidase']),
    ('CTL',   'Toxin', ['c-type lectin']),
    ('CRiSP', 'Toxin', ['crisp isoform']),
]

# ── Non-venom / housekeeping rules ────────────────────────────────────────────
NON_VENOM_RULES: list[tuple[str, str, list[str]]] = [
    ('Hemoglobin', 'Other', ['hemoglobin', 'globin', 'myoglobin']),
    ('Structural', 'Other', ['collagen', 'keratin', 'fibronectin', 'laminin',
                              'elastin', 'fibulin', 'microfibril', 'lamin',
                              'titin', 'dermatopontin', 'fibrillin',
                              'extracellular matrix protein', 'decorin',
                              'proteoglycan', 'tenascin']),
    ('Muscular',   'Other', ['troponin', 'tropomyosin', 'myosin', 'actin',
                              'dystrophin', 'myomegalin', 'myomesin',
                              'myozenin', 'sarcomere']),
    ('Metabolic',  'Other', ['pyruvate kinase', 'malate dehydrogenase',
                              'glyceraldehyde', 'phosphogluco',
                              'enolase', 'aldolase',
                              'lactate dehydrogenase',
                              'isocitrate dehydrogenase',
                              'aspartate aminotransferase',
                              'alanine aminotransferase',
                              'glutamate dehydrogenase',
                              'succinate dehydrogenase',
                              'adenylate kinase']),
    ('Regulatory', 'Other', ['calreticulin', 'calmodulin',
                              'immunoglobulin', 'fibrinogen',
                              'complement', 'annexin',
                              'ubiquitin', 'proteasome',
                              'heat shock protein', 'chaperone',
                              'ferritin', 'transferrin', 'albumin', 'serpin']),
    ('Other',      'Other', ['ribosomal protein', 'ribosomal',
                              'histone', 'tubulin', 'vimentin',
                              'cytochrome', 'mitochondrial',
                              'atpase', 'atp synthase', 'rna polymerase',
                              'helicase', 'actin', 'glucose-regulated protein',
                              'peptidyl-prolyl isomerase',
                              'peptidylprolyl isomerase',
                              'phosphoglucomutase', 'serotransferrin',
                              'pyruvate kinase', 'glutathione',
                              'sulfhydryl oxidase', 'protein disulfide-isomerase',
                              'splicing factor', 'zinc finger protein',
                              'f-box protein', 'ring finger protein',
                              'snrna', 'coiled-coil', 'lamin b',
                              'intrinsic factor', 'sodium/potassium',
                              'aquaglyceroporin', 'ubiquitin', 'polycomb']),
]

_VENOM_FAM_SET = {r[0] for r in VENOM_RULES}


# ═════════════════════════════════════════════════════════════════════════════
# SECTION B — Classification helpers
# ═════════════════════════════════════════════════════════════════════════════

def build_keyword_scorer(uw_path: Path,
                          top_n: int = 50,
                          min_n: int = 3) -> dict[str, list[tuple[str, float]]]:
    """
    Build per-family keyword scoring lists from unique_words__Family.csv.
    Weight = n_headers_in_group × word_length (longer phrases = more specific).
    Returns {family: [(keyword_lower, weight), ...]}
    """
    uw = pd.read_csv(uw_path)
    uw['Word'] = uw['Word'].astype(str).str.lower().str.strip()
    uw = uw[uw['n_headers_in_group'] >= min_n]
    uw['word_len'] = uw['Word'].str.split().str.len()
    uw['weight'] = uw['n_headers_in_group'] * uw['word_len']

    scorer: dict[str, list[tuple[str, float]]] = {}
    for fam, grp in uw.groupby('Group'):
        top = grp.nlargest(top_n, 'weight')[['Word', 'weight']].values.tolist()
        scorer[str(fam)] = [(w, float(wt)) for w, wt in top]
    return scorer


def _apply_rules(text_lower: str,
                 rules: list[tuple[str, str, list[str]]]
                 ) -> tuple[str, str] | None:
    """Return (Family, Function) for first matching rule, or None."""
    for family, function, keywords in rules:
        for kw in keywords:
            if kw in text_lower:
                return family, function
    return None


def _score_text(text_lower: str,
                scorer: dict[str, list[tuple[str, float]]]) -> dict[str, float]:
    scores: dict[str, float] = {}
    for fam, kw_list in scorer.items():
        s = sum(wt for kw, wt in kw_list if kw in text_lower)
        if s > 0:
            scores[fam] = s
    return scores


def _full_text(fasta_headers: str) -> str:
    """Lowercase all semicolon-joined fasta headers for matching."""
    return str(fasta_headers).lower()


def classify_row(text_lower: str,
                 scorer: Optional[dict[str, list[tuple[str, float]]]],
                 min_score: float = 3.0
                 ) -> tuple[str, str, str]:
    """
    Classify a single protein by text.
    Returns (Family, Function, tier_label).
    """
    # Tier 1: hard venom rules
    hit = _apply_rules(text_lower, VENOM_RULES)
    if hit:
        return hit[0], hit[1], 'tier1_venom'

    # Tier 1b: non-venom rules
    hit = _apply_rules(text_lower, NON_VENOM_RULES)
    if hit:
        return hit[0], hit[1], 'tier1_nonvenom'

    # Tier 2: unique-words scoring
    if scorer:
        scores = _score_text(text_lower, scorer)
        if scores:
            v_scores = {k: v for k, v in scores.items()
                        if k in _VENOM_FAM_SET}
            best_any   = max(scores, key=scores.get)
            best_v     = max(v_scores, key=v_scores.get) if v_scores else None
            best_v_sc  = v_scores.get(best_v, 0) if best_v else 0
            best_any_sc = scores[best_any]

            if best_any_sc >= min_score:
                if (best_v and best_v_sc >= min_score and
                        best_v_sc >= 0.4 * best_any_sc):
                    return best_v, 'Toxin', 'tier2_venom'
                else:
                    fn = 'Toxin' if best_any in _VENOM_FAM_SET else 'Other'
                    return best_any, fn, 'tier2_nonvenom'

    return 'Other', 'Other', 'unmatched'


# ═════════════════════════════════════════════════════════════════════════════
# SECTION C — Pipeline stages
# ═════════════════════════════════════════════════════════════════════════════

COLLECTION_META_COLS: list[str] = [
    'D E Venom Collection Number', 'Family', 'Genus', 'Species', 'Species2',
    'Form_short', 'Sex', 'Date collected', 'SVL', 'Tail length',
]
ANNOT_COLS: list[str] = ['Protein', 'Family', 'Group', 'Subgroup',
                         'Function', 'Notes', 'UniProt_ID']


# ── iBAQ column helper ────────────────────────────────────────────────────────

# MaxQuant writes an "iBAQ peptides" column (peptide count, not intensity) that
# starts with "iBAQ " and must be excluded from sample-iBAQ selection.
_IBAQ_EXCLUDE = {'iBAQ peptides'}


def ibaq_cols_of(df: pd.DataFrame) -> list[str]:
    """Return sample iBAQ intensity column names, excluding MaxQuant meta-columns."""
    return [c for c in df.columns
            if c.startswith('iBAQ ') and c not in _IBAQ_EXCLUDE]


# ── Stage 0: filter flagged proteins ──────────────────────────────────────────
#
# Removes proteins flagged by MaxQuant as decoy hits, lab contaminants, or
# identifications supported only by a modification site. Runs before the
# walk/annotate/merge stages so every downstream output (masterfile,
# MassSpecData, summaries, pie charts) is built from the cleaned data.

FILTER_FLAG_COLUMNS: list[str] = [
    'Reverse',
    'Potential contaminant',
    'Only identified by site',
]


def find_raw_proteingroups_files(root: Path) -> list[Path]:
    """
    Find proteinGroups*.txt files under root, excluding any that already
    look like filter-stage output (so re-running the pipeline on an
    output directory doesn't re-filter already-filtered files).
    """
    if not root.exists():
        raise FileNotFoundError(f"Input directory not found: {root}")
    return sorted(p for p in root.rglob('proteinGroups*.txt')
                   if p.is_file() and '_filtered' not in p.stem)


def filter_proteingroups_file(filepath: Path, output_root: Path) -> dict:
    """
    Remove rows where Reverse, Potential contaminant, or Only identified by
    site is '+', and save the cleaned table under output_root, mirroring the
    species subfolder structure with the same filename
    (output_root/<species>/proteinGroups.txt) so output_root can be walked
    by find_proteingroups_files() exactly like a normal MaxQuant run root.
    """
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    total_before = len(df)

    breakdown: dict[str, int | str] = {}
    for col in FILTER_FLAG_COLUMNS:
        if col not in df.columns:
            breakdown[col] = 'column not found'
            continue
        mask = df[col].astype(str).str.strip() == '+'
        breakdown[col] = int(mask.sum())
        df = df[~mask]

    total_after = len(df)

    # Name the output folder after the species/run folder (parent of the
    # proteinGroups file). Fall back to the file's own stem if the parent
    # has no useful name (e.g. file sits directly in the input root).
    species_name = filepath.resolve().parent.name
    species_name = species_name if species_name and species_name != '.' else filepath.stem
    species_out_dir = output_root / species_name
    species_out_dir.mkdir(parents=True, exist_ok=True)
    out_path = species_out_dir / filepath.name

    df.to_csv(out_path, sep='\t', index=False)

    return {
        'file': filepath.name,
        'source': str(filepath.parent),
        'before': total_before,
        'after': total_after,
        'removed': total_before - total_after,
        'breakdown': breakdown,
        'saved_to': str(out_path),
    }


def run_filter_stage(input_dir: Path, output_root: Path) -> Path:
    """
    Stage 0 — remove flagged proteins from every proteinGroups*.txt under
    input_dir, writing cleaned copies into output_root with the same folder
    structure and filenames. Returns output_root, which downstream stages
    walk in place of input_dir.
    """
    files = find_raw_proteingroups_files(input_dir)
    if not files:
        sys.exit(f"[filter] No proteinGroups*.txt files found under {input_dir}")

    print(f"[filter] found {len(files)} proteinGroups file(s) under {input_dir}")
    print(f"[filter] removing flagged proteins: {', '.join(FILTER_FLAG_COLUMNS)}")

    total_removed = 0
    total_kept = 0
    for f in files:
        result = filter_proteingroups_file(f, output_root)
        total_removed += result['removed']
        total_kept += result['after']
        details = ', '.join(
            f"{col} {n}" for col, n in result['breakdown'].items()
            if isinstance(n, int) and n > 0
        )
        suffix = f" ({details})" if details else ''
        print(f"[filter]   {result['source']}/{result['file']}: "
              f"{result['before']} -> {result['after']} proteins{suffix}")

    print(f"[filter] done — {total_kept} proteins kept, {total_removed} removed "
          f"across {len(files)} file(s)")
    print(f"[filter] filtered runs saved to: {output_root}")
    return output_root


# ── Stage 1: walk ─────────────────────────────────────────────────────────────

def find_proteingroups_files(root: Path) -> list[Path]:
    if not root.exists():
        raise FileNotFoundError(f"Input directory not found: {root}")
    return sorted(p for p in root.rglob('proteinGroups*.txt') if p.is_file())


def label_for_file(path: Path) -> str:
    parts = path.parts
    parent = path.parent.name
    grandparent = path.parent.parent.name if len(parts) >= 3 else ''
    if parent.lower() in ('txt', 'combined', 'uploads', '.', ''):
        if grandparent and grandparent.lower() not in ('combined', 'uploads', '.', ''):
            return grandparent
    return parent or path.stem


# ── Stage 2: annotate (Tier 0 — exact match + accession fallback) ────────────

def _lead_accession(fasta_header: str) -> str:
    """Return the first UniProt accession (|XXXXXX|) from a Fasta header."""
    m = re.search(r'\|([A-Z0-9]+)\|', str(fasta_header))
    return m.group(1) if m else ''


def _clean_fasta_header(h: str) -> str:
    """Strip leading/trailing whitespace and stray semicolons from a header."""
    return str(h).strip().lstrip(';').strip()


def load_master_annotation(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load the master annotation CSV and return two lookup tables:
      primary   — indexed by exact (cleaned) Fasta header string
      secondary — indexed by lead UniProt accession (fallback when header
                  strings differ between MaxQuant runs due to group composition
                  changes or stray leading semicolons in the annotation file)
    """
    df = pd.read_csv(path, low_memory=False)

    # Normalise Fasta headers: strip whitespace and leading semicolons so that
    # entries like ";tr|XXXX|..." match the clean "tr|XXXX|..." from MaxQuant.
    df['Fasta headers'] = df['Fasta headers'].astype(str).apply(_clean_fasta_header)
    df = df.drop_duplicates(subset='Fasta headers', keep='first')

    primary = df.set_index('Fasta headers')

    # Build secondary lookup keyed by lead UniProt accession
    df2 = df.copy()
    df2['_acc'] = df2['Fasta headers'].apply(_lead_accession)
    df2 = df2[df2['_acc'] != ''].drop_duplicates(subset='_acc', keep='first')
    secondary = df2.set_index('_acc')

    return primary, secondary


def annotate_proteingroups(pg_file: Path,
                           master_annot: pd.DataFrame,
                           master_accession: pd.DataFrame) -> pd.DataFrame:
    pg = pd.read_csv(pg_file, sep='\t', low_memory=False)
    rename_map = {c: 'Fasta headers' for c in pg.columns
                  if c.strip().lower().replace('.', ' ') == 'fasta headers'}
    pg = pg.rename(columns=rename_map)
    if 'Fasta headers' not in pg.columns:
        raise ValueError(f"No 'Fasta headers' column in {pg_file}")

    ibaq_cols = ibaq_cols_of(pg)
    pg = pg[['Fasta headers'] + ibaq_cols].copy()
    pg = pg[~pg['Fasta headers'].astype(str).str.startswith(('REV_', 'CON_'))]

    # ── Tier 0a: exact Fasta header match ─────────────────────────────────────
    for col in ANNOT_COLS:
        pg[col] = (pg['Fasta headers'].map(master_annot[col])
                   if col in master_annot.columns else '')

    # ── Tier 0b: fallback — match on lead UniProt accession ───────────────────
    # Catches cases where the protein-group composition differs between MaxQuant
    # runs (e.g. extra proteins appended, or a leading ";" in the annotation
    # file) so the exact string doesn't match even though it's the same protein.
    blank = pg['Family'].isna() | (pg['Family'].astype(str).str.strip() == '')
    if blank.any() and not master_accession.empty:
        pg['_lead_acc'] = pg['Fasta headers'].apply(_lead_accession)
        for col in ANNOT_COLS:
            if col not in master_accession.columns:
                continue
            fallback_vals = pg.loc[blank, '_lead_acc'].map(master_accession[col])
            filled = blank & fallback_vals.notna()
            if filled.any():
                pg.loc[filled, col] = fallback_vals[filled]
        n_rescued = int((blank & (pg['Family'].notna() & (pg['Family'].astype(str).str.strip() != ''))).sum())
        if n_rescued:
            print(f"           Tier 0b: rescued {n_rescued} proteins via accession fallback")
        pg = pg.drop(columns=['_lead_acc'], errors='ignore')

    return pg


# ── Stage 3: merge ────────────────────────────────────────────────────────────

def merge_all(annotated_dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Merge all annotated DataFrames into one masterfile.

    iBAQ columns (unique per file) are joined column-wise.
    Annotation columns (Family, Protein, …) are merged row-wise so that
    every protein keeps its annotation regardless of which file it came from.
    Using pd.concat(axis=1) + duplicated() would silently drop annotations
    for proteins that don't appear in the first file in the list.
    """
    t0 = time.time()
    ibaq_frames:  list[pd.DataFrame] = []
    annot_rows:   list[pd.DataFrame] = []
    seen_samples: set[str] = set()
    _EMPTY = {'', 'nan', 'None', 'NaN', 'NaT'}

    for label, df in annotated_dfs.items():
        ibaq_cols = ibaq_cols_of(df)
        rename_map = {}
        for c in ibaq_cols:
            if c in seen_samples:
                rename_map[c] = f'{c} [{label}]'
            else:
                seen_samples.add(c)
        if rename_map:
            df = df.rename(columns=rename_map)
            seen_samples.update(rename_map.values())

        df = df.drop_duplicates(subset='Fasta headers', keep='first')

        # iBAQ: each file contributes unique sample columns → concat column-wise
        current_ibaq = list(ibaq_cols_of(df))
        ibaq_frames.append(df[['Fasta headers'] + current_ibaq]
                           .set_index('Fasta headers'))

        # Annotation: collect all rows → will be coalesced row-wise below
        annot_present = [c for c in ANNOT_COLS if c in df.columns]
        if annot_present:
            annot_rows.append(df[['Fasta headers'] + annot_present])

    print(f"[merge]   merging {len(ibaq_frames)} file(s), "
          f"~{sum(len(d) for d in ibaq_frames)} rows ...", flush=True)

    # ── iBAQ: column-wise join ────────────────────────────────────────────────
    ibaq_merged = pd.concat(ibaq_frames, axis=1, sort=False, copy=False)

    # ── Annotation: row-wise coalesce ─────────────────────────────────────────
    # For each unique Fasta header, take the first non-null/non-empty value
    # across all files for each annotation column.
    if annot_rows:
        all_annot = pd.concat(annot_rows, axis=0, ignore_index=True)
        for col in ANNOT_COLS:
            if col in all_annot.columns:
                all_annot[col] = all_annot[col].astype(str).replace(_EMPTY, None)
                all_annot[col] = all_annot[col].where(all_annot[col].notna(), None)
        # groupby preserves insertion order; first() picks the first non-null value
        annot_merged = (all_annot
                        .groupby('Fasta headers', sort=False)
                        .first())
    else:
        annot_merged = pd.DataFrame(index=ibaq_merged.index)

    # ── Combine ───────────────────────────────────────────────────────────────
    merged = pd.concat([ibaq_merged, annot_merged], axis=1, sort=False)
    merged = merged.reset_index()

    ibaq_cols = list(ibaq_cols_of(merged))
    annot_present = [c for c in ANNOT_COLS if c in merged.columns]
    merged = merged[['Fasta headers'] + ibaq_cols + annot_present]

    print(f"[merge]   done in {time.time()-t0:.1f}s  "
          f"({len(merged)} rows × {merged.shape[1]} cols)", flush=True)
    return merged


# ── Stage 3.25: re-annotate masterfile from master annotation ─────────────────
# Needed when --masterfile is used: the old masterfile may have stale
# 'auto:unmatched' Notes (and wrong Family=Other) for proteins that ARE in the
# master annotation.  Running Tier 0 again fixes both Family and Notes.

def reannotate_from_master(master: pd.DataFrame,
                           master_annot: pd.DataFrame,
                           master_accession: pd.DataFrame) -> tuple[pd.DataFrame, int, int]:
    """
    For every row in *master* whose Notes starts with 'auto:', check whether
    the protein is actually present in *master_annot* (Tier 0a) or
    *master_accession* (Tier 0b).  If found, copy the annotation columns from
    the master and clear the stale 'auto:' Note.

    Returns (updated_df, n_fixed_exact, n_fixed_accession).
    """
    df = master.copy()
    if 'Notes' not in df.columns:
        return df, 0, 0

    stale = df['Notes'].fillna('').astype(str).str.startswith('auto:')
    if not stale.any():
        return df, 0, 0

    # Ensure annotation columns are object dtype so string assignment doesn't warn
    for col in ANNOT_COLS:
        if col in df.columns and df[col].dtype != object:
            df[col] = df[col].astype(object)

    # Tier 0a: exact Fasta header match
    exact_match = stale & df['Fasta headers'].isin(master_annot.index)
    n_exact = int(exact_match.sum())
    if n_exact:
        for col in ANNOT_COLS:
            if col == 'Notes' or col not in master_annot.columns:
                continue
            vals = df.loc[exact_match, 'Fasta headers'].map(master_annot[col])
            df.loc[exact_match, col] = vals.values
        # Clear stale Notes; keep blank (master's Notes for these rows is usually null)
        df.loc[exact_match, 'Notes'] = master_annot['Notes'].reindex(
            df.loc[exact_match, 'Fasta headers'].values
        ).fillna('').values if 'Notes' in master_annot.columns else ''

    # Tier 0b: lead UniProt accession fallback (for still-stale rows)
    still_stale = stale & ~exact_match
    if still_stale.any() and not master_accession.empty:
        df['_lead_acc'] = df['Fasta headers'].apply(_lead_accession)
        acc_match = still_stale & df['_lead_acc'].isin(master_accession.index)
        n_acc = int(acc_match.sum())
        if n_acc:
            for col in ANNOT_COLS:
                if col == 'Notes' or col not in master_accession.columns:
                    continue
                vals = df.loc[acc_match, '_lead_acc'].map(master_accession[col])
                df.loc[acc_match, col] = vals.values
            note_vals = master_accession['Notes'].reindex(
                df.loc[acc_match, '_lead_acc'].values
            ).fillna('').values if 'Notes' in master_accession.columns else [''] * n_acc
            df.loc[acc_match, 'Notes'] = note_vals
        df = df.drop(columns=['_lead_acc'], errors='ignore')
    else:
        n_acc = 0

    return df, n_exact, n_acc


# ── Stage 3.5: auto-classify proteins missing a Family ────────────────────────

def auto_classify_masterfile(
    master: pd.DataFrame,
    scorer: Optional[dict[str, list[tuple[str, float]]]],
    min_score: float = 3.0,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """
    For every row where Family is blank/NaN apply Tier 1 + Tier 2 classification.
    Rows already annotated by exact-match (Tier 0) are left untouched.

    Returns the updated masterfile and a dict with classification counts.
    """
    df = master.copy()

    # Identify rows that need classification
    blank_mask = df['Family'].isna() | (df['Family'].astype(str).str.strip() == '')

    counts: dict[str, int] = {
        'already_annotated': int((~blank_mask).sum()),
        'tier1_venom':    0,
        'tier1_nonvenom': 0,
        'tier2_venom':    0,
        'tier2_nonvenom': 0,
        'unmatched':      0,
    }

    if blank_mask.sum() == 0:
        print("[classify] all proteins already have a Family — nothing to do.")
        return df, counts

    print(f"[classify] {blank_mask.sum()} proteins lack a Family → running "
          f"Tier 1 (keyword rules) + Tier 2 (unique-words scorer) ...", flush=True)

    rows_to_classify = df[blank_mask].index.tolist()

    families, functions, tiers, notes = [], [], [], []
    for idx in rows_to_classify:
        fh    = str(df.at[idx, 'Fasta headers'])
        text  = _full_text(fh)
        fam, fn, tier = classify_row(text, scorer, min_score)
        families.append(fam)
        functions.append(fn)
        tiers.append(tier)
        notes.append(f'auto:{tier}')
        counts[tier] = counts.get(tier, 0) + 1

    df.loc[blank_mask, 'Family']   = families
    df.loc[blank_mask, 'Function'] = functions
    df.loc[blank_mask, 'Notes']    = notes
    # Leave Group / Subgroup / UniProt_ID blank for auto-classified rows
    # (Group = same as Family for venom toxins, but we keep it blank so the
    #  user can manually curate if needed)

    return df, counts


# ── Family-name normalisation ─────────────────────────────────────────────────

# Maps legacy / misspelled / abbreviated family names to current canonical names.
# Applied once after the masterfile is loaded so Tier 0 rows are also covered.
FAMILY_RENAME: dict[str, str] = {
    'Procoaculant':        'Prothrombin activator',
    'Procoagulant':        'Prothrombin activator',
    'Procoagulant toxin':  'Prothrombin activator',
    'Kaz':                 'Kazal',
}


def normalise_family_names(df: pd.DataFrame) -> pd.DataFrame:
    """Rename legacy Family values to their current canonical names."""
    if 'Family' not in df.columns:
        return df
    df = df.copy()
    df['Family'] = df['Family'].replace(FAMILY_RENAME)
    return df


# ── Stage 4: MassSpecData.csv ─────────────────────────────────────────────────

def load_collection(path: Path) -> pd.DataFrame:
    """
    Load sample metadata from a collection spreadsheet (MASTER sheet).

    Supports two formats:

    Old format  (e.g. Esquerre_Venom_Collection.xlsx)
        Must have an 'iBAQ_name' column containing the full iBAQ column name,
        e.g. 'iBAQ DEVC_002_Acanthophis_antarcticus_MilkedVenom_Nov2021_...'

    New format  (e.g. LAB_WORK_DAMIEN_VENOMICS.xlsx)
        Has an 'Experiment_Name' column containing the same string WITHOUT the
        leading 'iBAQ ' prefix, e.g.
        'DEVC_002_Acanthophis_antarcticus_MilkedVenom_Nov2021_...'
        The loader prepends 'iBAQ ' automatically.
        Column mapping:
          Experiment_Name → iBAQ_name index
          Family          → Tax_Family  (taxonomic family, e.g. Elapidae)
          Genus, Species  → used directly; Species2 derived if absent
          Form_short, Date collected, D E Venom Collection Number → used directly
          Sex, SVL, Tail length → filled with '' if absent in sheet
    """
    xl = pd.read_excel(path, sheet_name='MASTER', dtype=str)

    if 'iBAQ_name' in xl.columns:
        # ── Old format ────────────────────────────────────────────────────────
        xl = xl[xl['iBAQ_name'].notna()].copy()
        xl = xl.set_index('iBAQ_name')

    elif 'Experiment_Name' in xl.columns:
        # ── New lab-work format ───────────────────────────────────────────────
        xl = xl[xl['Experiment_Name'].notna()].copy()
        xl['iBAQ_name'] = 'iBAQ ' + xl['Experiment_Name'].str.strip()
        # Drop duplicate iBAQ names (same experiment listed twice in the sheet)
        xl = xl.drop_duplicates(subset='iBAQ_name')
        xl = xl.set_index('iBAQ_name')

    else:
        raise ValueError(
            f"MASTER sheet in {path} has neither an 'iBAQ_name' column nor an "
            "'Experiment_Name' column.  Cannot build the sample-metadata lookup.\n"
            "  Old-format files need: iBAQ_name\n"
            "  New-format files need: Experiment_Name"
        )

    xl = xl.rename(columns={'Family': 'Tax_Family'})

    # Derive Species2 (binomial) when it is absent
    if 'Species2' not in xl.columns:
        xl['Species2'] = (
            xl.get('Genus',   pd.Series(dtype=str)).fillna('') + ' ' +
            xl.get('Species', pd.Series(dtype=str)).fillna('')
        ).str.strip()

    # Guarantee that all columns expected downstream are present
    for col in ('Sex', 'SVL', 'Tail length'):
        if col not in xl.columns:
            xl[col] = ''

    return xl


def build_massspecdata(master: pd.DataFrame,
                       collection: pd.DataFrame) -> pd.DataFrame:
    ibaq_cols = ibaq_cols_of(master)
    fasta_headers = master['Fasta headers'].tolist()
    n_proteins = len(fasta_headers)
    n_meta = len(COLLECTION_META_COLS)

    header_row = ['Fasta headers'] + COLLECTION_META_COLS + fasta_headers
    annot_rows = []
    for col in ANNOT_COLS:
        vals = (master[col].fillna('').astype(str).replace('nan', '').tolist()
                if col in master.columns else [''] * n_proteins)
        annot_rows.append([col] + [''] * n_meta + vals)

    meta_map = {
        'D E Venom Collection Number': 'D E Venom Collection Number',
        'Family': 'Tax_Family', 'Genus': 'Genus', 'Species': 'Species',
        'Species2': 'Species2', 'Form_short': 'Form_short', 'Sex': 'Sex',
        'Date collected': 'Date collected', 'SVL': 'SVL',
        'Tail length': 'Tail length',
    }
    sample_rows = []
    for ibaq_col in ibaq_cols:
        ibaq_vals = master[ibaq_col].tolist()
        key = ibaq_col.split(' [')[0]
        if key in collection.index:
            row = collection.loc[key]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            meta_vals = [row.get(meta_map[c], '') for c in COLLECTION_META_COLS]
        else:
            devc_m = re.search(r'DEVC[_\s](\d+)', ibaq_col, re.IGNORECASE)
            devc = f"DEVC {int(devc_m.group(1)):03d}" if devc_m else ''
            meta_vals = [devc] + [''] * (n_meta - 1)
        sample_rows.append([ibaq_col] + list(meta_vals) + ibaq_vals)

    return pd.DataFrame([header_row] + annot_rows + sample_rows)


# ── Stage 5: iBAQ proportions ─────────────────────────────────────────────────

def compute_proportions(master: pd.DataFrame,
                        venom_only: bool,
                        no_pla2_inhibitor: bool
                        ) -> tuple[pd.DataFrame, list[str]]:
    df = master.copy()
    if venom_only:
        # Keep rows annotated as venom by either convention used in the codebase.
        df = df[df['Function'].fillna('').astype(str).str.strip().isin(
            {'Venom', 'Toxin'})]
    if no_pla2_inhibitor:
        df = df[df['Family'] != 'PLA2 inhibitor']

    df['Family'] = df['Family'].fillna('Other').replace('', 'Other')
    ibaq_cols = ibaq_cols_of(df)

    present = df['Family'].dropna().unique().tolist()
    ordered = [f for f in FAMILY_ORDER if f in present]
    extras  = sorted([f for f in present if f not in ordered and f != 'Other'])
    families = ordered + extras
    if 'Other' in present and 'Other' not in families:
        families.append('Other')

    fam_totals = df.groupby('Family')[ibaq_cols].apply(
        lambda g: g.apply(pd.to_numeric, errors='coerce').fillna(0).sum()
    )
    col_totals = fam_totals.sum(axis=0).replace(0, np.nan)
    pct = fam_totals.div(col_totals, axis=1) * 100

    out = pct.T.reset_index().rename(columns={'index': 'ibaq_name'})
    if 'ibaq_name' not in out.columns:
        out.columns = ['ibaq_name'] + list(out.columns[1:])
    out = out.dropna(subset=families, how='all').fillna(0)
    return out[['ibaq_name'] + families], families


def compute_venom_percent(master: pd.DataFrame,
                          no_pla2_inhibitor: bool = False) -> dict[str, float]:
    """
    For each iBAQ sample column, return the percentage of total iBAQ signal
    that comes from venom (Toxin) proteins.  Non-toxin proteins (Function ==
    'Other', e.g. Structural, Metabolic, Hemoglobin …) are the denominator
    complement.
    """
    ibaq_cols = ibaq_cols_of(master)
    df = master.copy()
    if no_pla2_inhibitor:
        df = df[df['Family'] != 'PLA2 inhibitor']
    for c in ibaq_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)

    # Accept both 'Venom' (master annotation CSV convention) and 'Toxin'
    # (auto-classifier convention).  Anything else (Other, blank, AMBIGUOUS)
    # is treated as non-venom.
    venom_mask = df['Function'].fillna('').astype(str).str.strip().isin(
        {'Venom', 'Toxin'})
    totals      = df[ibaq_cols].sum()
    venom_sums  = df.loc[venom_mask, ibaq_cols].sum()

    result: dict[str, float] = {}
    for c in ibaq_cols:
        t = float(totals[c])
        result[c] = round(float(venom_sums[c]) / t * 100, 2) if t > 0 else 0.0
    return result


# ── Stage 6: summary tables ───────────────────────────────────────────────────

def _year_only(val) -> str:
    """Extract a 4-digit year from a date string or return the value as-is."""
    s = str(val).strip()
    m = re.search(r'\b(1[89]\d\d|20\d\d)\b', s)
    return m.group(1) if m else (s if s not in ('', 'nan', 'NaT', 'None') else '')


# Canonical ordering of run-metadata columns inserted into the summary.
RUN_META_COLS = [
    'Sample_type', 'Run_date', 'Pilot_run', 'Run_time_min', 'Filter_kDa', 'Acquisition',
]


def _parse_run_metadata(ibaq_name: str) -> dict:
    """
    Extract LC-MS/MS run parameters encoded in the iBAQ column name.

    Expected name pattern (underscore-delimited after the 'iBAQ ' prefix):
        DEVC_NNN_Genus_species_<SampleType>_<MonYYYY>[_PilotRun]_<N>min_<K>kDa_<ACQ>

    Returns a dict with keys matching RUN_META_COLS.
    """
    s = ibaq_name.replace('iBAQ ', '').strip()

    # Sample type
    sample_type = ''
    for kw in ('MilkedVenom', 'PreservedFixedGland', 'FrozenEthanolGland'):
        if kw in s:
            sample_type = kw
            break

    # Run date  e.g.  Nov2021,  Jan2023,  Feb2026
    m_date = re.search(
        r'(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)(\d{4})', s)
    run_date = (m_date.group(1) + ' ' + m_date.group(2)) if m_date else ''

    # Pilot run flag
    pilot_run = 'Yes' if 'PilotRun' in s else 'No'

    # Gradient run time  e.g.  30min  120min
    m_time = re.search(r'(\d+)min', s)
    run_time_min = int(m_time.group(1)) if m_time else ''

    # Molecular-weight filter  e.g.  10kDa  30kDa
    m_filter = re.search(r'(\d+)kDa', s)
    filter_kda = int(m_filter.group(1)) if m_filter else ''

    # Data-acquisition mode  e.g.  DDA  DIA  PRM  (may appear at end of string)
    m_acq = re.search(r'(?<![A-Za-z])(DDA|DIA|PRM|SWATH)(?![A-Za-z])', s, re.IGNORECASE)
    acquisition = m_acq.group(1).upper() if m_acq else ''

    return {
        'Sample_type':   sample_type,
        'Run_date':      run_date,
        'Pilot_run':     pilot_run,
        'Run_time_min':  run_time_min,
        'Filter_kDa':    filter_kda,
        'Acquisition':   acquisition,
    }


def build_summary_table(prop_df, families, collection,
                        venom_pct: dict[str, float] | None = None,
                        sample_universe_df: pd.DataFrame | None = None):
    """
    Build the per-sample summary table.

    Parameters
    ----------
    prop_df : DataFrame
        Venom-only family proportions (rows = samples that have ≥1 venom protein).
        May be empty for samples/runs with no classified venom proteins.
    families : list[str]
        Toxin-family column names present in prop_df.
    collection : DataFrame
        MASTER collection metadata, indexed by iBAQ_name key.
    venom_pct : dict[str, float] | None
        Fraction of total iBAQ that is venom, keyed by iBAQ column name.
    sample_universe_df : DataFrame | None
        Proportions computed from ALL proteins (prop_all).  When supplied this
        is used as the definitive sample list so that samples with zero venom
        proteins (absent from prop_df) still appear in the summary — they will
        have 0 for every family column and Venom_pct drawn from venom_pct dict.
        If None, prop_df is used as the sample source (original behaviour).
    """
    # Determine which DataFrame drives the sample list.
    # sample_universe_df (prop_all) always contains every sample; prop_df
    # (prop_venom) only contains samples with ≥1 venom protein.
    iter_df = sample_universe_df if (sample_universe_df is not None
                                     and not sample_universe_df.empty) else prop_df

    if iter_df.empty:
        # Nothing to report at all — return an empty but well-formed DataFrame.
        meta_cols = ['Sample', 'D E Venom Collection Number', 'Family',
                     'Genus', 'Species', 'Species2', 'Form_short', 'Sex',
                     'Date collected', 'SVL', 'Tail length'] + RUN_META_COLS
        if venom_pct is not None:
            meta_cols.append('Venom_pct')
        return pd.DataFrame(columns=meta_cols + families)

    meta_rows = []
    for _, row in iter_df.iterrows():
        ibaq_name = row['ibaq_name']
        key = ibaq_name.split(' [')[0]
        meta: dict = {'Sample': ibaq_name}
        if key in collection.index:
            c = collection.loc[key]
            if isinstance(c, pd.DataFrame):
                c = c.iloc[0]
            meta['D E Venom Collection Number'] = c.get('D E Venom Collection Number', '')
            meta['Family']       = c.get('Tax_Family', '')
            meta['Genus']        = c.get('Genus', '')
            meta['Species']      = c.get('Species', '')
            meta['Species2']     = c.get('Species2', '')
            meta['Form_short']   = c.get('Form_short', '')
            meta['Sex']          = c.get('Sex', '')
            meta['Date collected'] = _year_only(c.get('Date collected', ''))
            meta['SVL']          = c.get('SVL', '')
            meta['Tail length']  = c.get('Tail length', '')
        else:
            m = re.search(r'DEVC[_\s](\d+)', ibaq_name, re.IGNORECASE)
            meta['D E Venom Collection Number'] = (
                f"DEVC {int(m.group(1)):03d}" if m else '')
            for k in ('Family','Genus','Species','Species2','Form_short',
                      'Sex','Date collected','SVL','Tail length'):
                meta[k] = ''
        # Run-level metadata parsed from the iBAQ column name itself
        meta.update(_parse_run_metadata(ibaq_name))
        meta_rows.append(meta)

    meta_df = pd.DataFrame(meta_rows)

    # Build family-proportion lookup from prop_df (venom-only).
    # Samples absent from prop_df (no venom proteins) will be left-joined as
    # NaN and then zero-filled — correctly representing 0 % for every toxin family.
    if not prop_df.empty and families:
        fam_df = prop_df.rename(columns={'ibaq_name': 'Sample'})[['Sample'] + families]
    else:
        fam_df = pd.DataFrame(columns=['Sample'] + families)

    out = pd.merge(meta_df, fam_df, on='Sample', how='left')
    # Samples with no venom proteins get 0 % for every toxin-family column.
    out[families] = out[families].fillna(0)

    # Insert Venom_pct right after the last run-metadata column, before family
    # columns.  Fall back to end-of-DataFrame when families is empty.
    if venom_pct is not None:
        cols = list(out.columns)
        if families and families[0] in cols:
            insert_pos = cols.index(families[0])
        elif RUN_META_COLS and RUN_META_COLS[-1] in cols:
            insert_pos = cols.index(RUN_META_COLS[-1]) + 1
        else:
            insert_pos = len(cols)
        out.insert(
            insert_pos,
            'Venom_pct',
            out['Sample'].map(venom_pct).fillna(0).round(2),
        )
    return out


def build_averaged_summary(summary, families):
    s = summary.copy()
    s['Group_key'] = (s['Species2'].fillna('').astype(str) + ' | '
                      + s['Form_short'].fillna('').astype(str))
    agg = {f: 'mean' for f in families}
    agg['Sample'] = 'count'
    out = (s.groupby('Group_key', sort=True)
             .agg({**agg, 'Family': 'first', 'Genus': 'first',
                   'Species': 'first', 'Species2': 'first',
                   'Form_short': 'first'})
             .rename(columns={'Sample': 'n_samples'})
             .reset_index())
    front = ['Group_key','Family','Genus','Species','Species2','Form_short','n_samples']
    return out[front + families]


# ── Stage 7: pie charts ───────────────────────────────────────────────────────

def _get_color(family: str) -> str:
    return COLOR_MAPPING.get(family, '#CCCCCC')


def plot_pie(values, labels, title, out_path):
    """
    Render a single publication-quality pie chart.

    All percentage information is displayed in the legend (colour swatch +
    family name + value).  No in-wedge text is drawn: autopct labels inside
    narrow wedges produce crossing centre-lines in matplotlib and are
    unnecessary given the legend already carries every value.

    Parameters
    ----------
    values   : list[float]  – percentage values for each family
    labels   : list[str]    – family names matching COLOR_MAPPING keys
    title    : str          – chart title; may contain matplotlib LaTeX math,
                              e.g. r'$\\it{Genus\\ species}$\\nDetail line'
    out_path : Path | str   – output file (.pdf recommended for publication)
    """
    pairs  = [(v, l) for v, l in zip(values, labels) if v > 0]
    if not pairs:
        return
    vf     = [p[0] for p in pairs]
    lf     = [p[1] for p in pairs]
    colors = [_get_color(l) for l in lf]

    # ── Canvas ────────────────────────────────────────────────────────────────
    # 5 × 5 in is a comfortable single-column panel size for most journals.
    fig, ax = plt.subplots(figsize=(5, 5))

    # ── Pie wedges ────────────────────────────────────────────────────────────
    # labels=None and no autopct → clean wedges with no floating text.
    # linewidth=0 removes the wedge border lines entirely; white edges on thin
    # wedges converge at the pie centre and produce a visible spoke/crossing
    # artefact regardless of wedge size.
    ax.pie(
        vf,
        labels=None,
        colors=colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(linewidth=0),
    )

    # ── Legend ────────────────────────────────────────────────────────────────
    # Placed below the chart in two columns so it never overlaps the pie.
    # Each entry shows a colour swatch, the family name, and its percentage.
    legend_labels = [f'{l}  ({v:.1f}%)' for l, v in zip(lf, vf)]
    patches = [mpatches.Patch(facecolor=c, edgecolor='#cccccc', linewidth=0.5)
               for c in colors]
    ax.legend(
        patches, legend_labels,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.20),
        ncol=2,
        fontsize=8,
        frameon=True,
        framealpha=0.9,
        edgecolor='#cccccc',
    )

    # ── Title ─────────────────────────────────────────────────────────────────
    # linespacing gives breathing room when the title spans two lines
    # (species name on line 1, sample details on line 2).
    ax.set_title(title, fontsize=9, pad=10, linespacing=1.6)

    # ── Save ──────────────────────────────────────────────────────────────────
    # dpi=300 is the standard minimum for print-quality raster exports;
    # for PDF (vector) it has no effect on the curves themselves but some
    # journal submission systems inspect the dpi metadata.
    fig.savefig(out_path, format='pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)


def pie_title(ibaq_name, collection):
    key = ibaq_name.split(' [')[0]
    if key in collection.index:
        row = collection.loc[key]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        g = str(row.get('Genus','') or '').strip()
        s = str(row.get('Species','') or '').strip()
        f = str(row.get('Form_short','') or '').strip()
        d = str(row.get('D E Venom Collection Number','') or '').strip()
        return f'{g} {s} {f} {d}'.strip()
    return ibaq_name.replace('iBAQ ', '')


def _display_title(ibaq_name: str, collection: pd.DataFrame) -> str:
    """
    Build a matplotlib-formatted display title for a pie chart.

    Line 1 — Genus + Species italicised via matplotlib's LaTeX math renderer:
              r'$\\it{Genus\\ species}$'
    Line 2 — Form_short and DEVC collection number separated by ' | '.

    Falls back to the plain pie_title() string when metadata is missing,
    so the chart always has *some* title.
    """
    key = ibaq_name.split(' [')[0]
    if key in collection.index:
        row = collection.loc[key]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        g = str(row.get('Genus',   '') or '').strip()
        s = str(row.get('Species', '') or '').strip()
        f = str(row.get('Form_short', '') or '').strip()
        d = str(row.get('D E Venom Collection Number', '') or '').strip()

        species_part = f'{g} {s}'.strip()
        detail_part  = ' | '.join(x for x in [f, d] if x)

        if species_part:
            # Spaces inside LaTeX math mode must be escaped as '\ '
            latex_species = r'$\it{' + species_part.replace(' ', r'\ ') + r'}$'
            return f'{latex_species}\n{detail_part}' if detail_part else latex_species

        return detail_part or ibaq_name.replace('iBAQ ', '')

    return pie_title(ibaq_name, collection)  # plain-text fallback


def _safe_fn(s):
    return re.sub(r'[^\w\-_\. ]', '_', s).strip() or 'sample'


def generate_all_pies(prop_df, families, collection, out_dir,
                       subfolder: str = 'PieCharts',
                       filtered: bool = False,
                       species_avg_df=None):
    """
    Render per-sample pie charts into out_dir/subfolder/.
    If filtered=True also renders ≥1% versions into out_dir/subfolder_filteredOver1percent/.
    If species_avg_df is given, renders averaged pies into out_dir/subfolder_AveragePerSpecies/.
    """
    full_dir = out_dir / subfolder
    full_dir.mkdir(parents=True, exist_ok=True)
    filt_dir = out_dir / f'{subfolder}_filteredOver1percent' if filtered else None
    if filt_dir:
        filt_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    for _, row in prop_df.iterrows():
        ibaq_name     = row['ibaq_name']
        plain_title   = pie_title(ibaq_name, collection)       # plain text → filename
        display_title = _display_title(ibaq_name, collection)  # LaTeX → chart title
        values        = [float(row.get(f, 0)) for f in families]
        plot_pie(values, families, display_title,
                 full_dir / (_safe_fn(plain_title) + '_pie_chart.pdf'))
        n += 1
        if filt_dir:
            keep = [(v, f) for v, f in zip(values, families) if v >= 1.0]
            if keep:
                plot_pie([x[0] for x in keep], [x[1] for x in keep],
                         display_title,
                         filt_dir / (_safe_fn(plain_title) + '_filtered_pie_chart.pdf'))

    if species_avg_df is not None:
        avg_dir = out_dir / f'{subfolder}_AveragePerSpecies'
        avg_dir.mkdir(parents=True, exist_ok=True)
        for _, row in species_avg_df.iterrows():
            sp2  = str(row.get('Species2',   '') or '').strip()
            form = str(row.get('Form_short', '') or '').strip()
            plain_title = f'{sp2} {form}'.strip()
            # Italicise the species name for display; form info on a second line
            if sp2:
                latex_sp      = r'$\it{' + sp2.replace(' ', r'\ ') + r'}$'
                display_title = f'{latex_sp}\n{form}' if form else latex_sp
            else:
                display_title = plain_title
            values = [float(row.get(f, 0)) for f in families]
            plot_pie(values, families, display_title,
                     avg_dir / (_safe_fn(plain_title) + '_species_avg_pie.pdf'))
    return n


# ── Stage 8: unclassified report ─────────────────────────────────────────────

def flag_unclassified(master: pd.DataFrame, out_path: Path) -> int:
    """
    Write proteins that could not be classified by any method to out_path.

    Includes:
      • Rows where Family is still blank/NaN (should be rare after auto-classify)
      • Rows where Notes == 'auto:unmatched' — went through all three tiers
        (Tier 0 exact match, Tier 1 keyword rules, Tier 2 unique-words scoring)
        and still couldn't be matched to a specific family; assigned 'Other'
        as a fallback. These are the proteins that most need manual review.

    Excludes proteins that were successfully identified by Tier 1 or Tier 2
    (Notes == 'auto:tier1_venom', 'auto:tier1_nonvenom', 'auto:tier2_venom',
    'auto:tier2_nonvenom') even if their Family ended up as 'Other'.
    """
    notes = master['Notes'].astype(str).str.strip()
    blank_mask    = master['Family'].isna() | (master['Family'].astype(str).str.strip() == '')
    unmatched_mask = notes == 'auto:unmatched'
    mask = blank_mask | unmatched_mask

    annot_cols = [c for c in ['Fasta headers', 'Protein', 'Family', 'Notes', 'UniProt_ID']
                  if c in master.columns]
    unc = master.loc[mask, annot_cols].copy()
    unc = unc.drop_duplicates(subset='Fasta headers')

    # Filter out empty MaxQuant placeholder rows: NaN, blank, or semicolons-only
    fh = unc['Fasta headers'].astype(str).str.strip()
    unc = unc[~(fh.isin({'', 'nan', 'NaN', 'None'}) | fh.str.fullmatch(r'[;]+'))]

    if unc.empty:
        return 0

    def short_desc(h):
        first = str(h).split(';')[0]
        m = re.search(r'\|\w+\s+(.+?)\s+OS=', first)
        return m.group(1) if m else first[:120]

    unc.insert(1, 'Description', unc['Fasta headers'].apply(short_desc))
    unc.to_csv(out_path, index=False)
    return len(unc)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION D — Main orchestration
# ═════════════════════════════════════════════════════════════════════════════

def main() -> None:
    p = argparse.ArgumentParser(
        description='Unified venomics LC-MS/MS pipeline: '
                    'walk → annotate (exact match) → '
                    'auto-classify (keyword + unique-words) → '
                    'merge → summarise → plot.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # ── Inputs ────────────────────────────────────────────────────────────────
    p.add_argument('--input-dir',
                   help='Root directory to walk for proteinGroups*.txt.')
    p.add_argument('--master-annotation',
                   help='Master classification CSV '
                        '(e.g. AllProteinGroups_anotated3.csv).')
    p.add_argument('--unique-words',
                   help='unique_words__Family.csv — used for Tier 2 '
                        'keyword scoring of proteins not in the master '
                        'annotation. Optional but recommended.')
    p.add_argument('--collection', required=True,
                   help='Sample metadata spreadsheet. Accepts two formats: '
                        '(1) old format with an iBAQ_name column '
                        '(e.g. Esquerre_Venom_Collection.xlsx), or '
                        '(2) new lab-work format with an Experiment_Name column '
                        '(e.g. LAB_WORK_DAMIEN_VENOMICS.xlsx). '
                        'Both must have a MASTER sheet.')
    p.add_argument('--masterfile',
                   help='Existing protein_groups_masterfile.csv. '
                        'Providing this skips walk/annotate/merge/classify.')

    # ── Output ────────────────────────────────────────────────────────────────
    p.add_argument('--output-dir', default='.')
    p.add_argument('--filtered-dir',
                   help='Where to save Stage 0 filtered proteinGroups runs '
                        '(Reverse / Potential contaminant / Only identified '
                        'by site removed), mirroring the --input-dir folder '
                        'structure with proteinGroups.txt filenames unchanged. '
                        'Default: <output-dir>/FilteredProteinGroups. '
                        'All downstream stages (including pie charts) read '
                        'from this filtered data, not --input-dir.')

    # ── Stage toggles ─────────────────────────────────────────────────────────
    p.add_argument('--skip-filter',          action='store_true',
                   help='Skip Stage 0 protein filtering and walk --input-dir '
                        'directly (unfiltered, including Reverse/contaminant/'
                        'site-only rows).')
    p.add_argument('--skip-walk',            action='store_true')
    p.add_argument('--skip-annotate',        action='store_true')
    p.add_argument('--skip-classify',        action='store_true',
                   help='Skip Tier 1 + Tier 2 auto-classification step.')
    p.add_argument('--skip-merge',           action='store_true')
    p.add_argument('--skip-massspecdata',    action='store_true')
    p.add_argument('--skip-summary',         action='store_true')
    p.add_argument('--skip-averaged',        action='store_true')
    p.add_argument('--skip-plots',           action='store_true')
    p.add_argument('--skip-species-pies',    action='store_true')
    p.add_argument('--skip-unclassified',    action='store_true')

    # ── Analysis options ──────────────────────────────────────────────────────
    p.add_argument('--filtered-pies',        action='store_true',
                   help='Also generate ≥1%% filtered versions of every pie chart.')
    p.add_argument('--no-pla2-inhibitor',    action='store_true',
                   help='Exclude PLA2 inhibitor proteins from all proportions.')
    p.add_argument('--classify-min-score',   type=float, default=3.0,
                   help='Min Tier 2 score to accept a family (default 3.0).')
    p.add_argument('--classify-top-n',       type=int, default=50,
                   help='Max keywords per family for Tier 2 (default 50).')

    args = p.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --masterfile auto-skips walk/annotate/merge
    if args.masterfile:
        args.skip_walk = args.skip_annotate = args.skip_merge = True

    # ── Load unique-words scorer (optional) ───────────────────────────────────
    scorer: Optional[dict] = None
    if args.unique_words and not args.skip_classify:
        print(f"[classify] loading unique-words scorer: {args.unique_words}")
        scorer = build_keyword_scorer(
            Path(args.unique_words),
            top_n=args.classify_top_n,
        )
        print(f"           {len(scorer)} families in scorer")

    # ── Stage 0: filter flagged proteins (Reverse / contaminant / site-only) ──
    # Runs before walk so every downstream output — masterfile, MassSpecData,
    # summaries, and pie charts — is built from the cleaned data.
    walk_root: Optional[Path] = None
    if not args.skip_walk:
        if not args.input_dir:
            sys.exit("--input-dir required unless --skip-walk or --masterfile given")
        input_root = Path(args.input_dir)

        if not args.skip_filter:
            filtered_dir = (Path(args.filtered_dir) if args.filtered_dir
                           else out_dir / 'FilteredProteinGroups')
            walk_root = run_filter_stage(input_root, filtered_dir)
        else:
            print("[filter] skipped — walking --input-dir unfiltered")
            walk_root = input_root

    # ── Stage 1: walk ─────────────────────────────────────────────────────────
    files: list[Path] = []
    annotated: dict[str, pd.DataFrame] = {}
    if not args.skip_walk:
        root = walk_root
        files = find_proteingroups_files(root)
        print(f"[walk] found {len(files)} proteinGroups file(s) under {root}")
        if not files:
            sys.exit("No proteinGroups*.txt files found.")

    # ── Stage 2: annotate (Tier 0 exact match) ────────────────────────────────
    if not args.skip_annotate and files:
        if not args.master_annotation:
            sys.exit("--master-annotation required unless --skip-annotate")
        print(f"[annotate] loading: {args.master_annotation}")
        master_annot, master_accession = load_master_annotation(Path(args.master_annotation))
        print(f"           {len(master_annot)} annotated fasta headers in master "
              f"({len(master_accession)} unique accessions for fallback)")

        for f in files:
            label = label_for_file(f)
            print(f"[annotate]   {label}  ({f.name})")
            annotated[label] = annotate_proteingroups(f, master_annot, master_accession)

    # ── Stage 3: merge ────────────────────────────────────────────────────────
    if not args.skip_merge and annotated:
        print("[merge] combining into protein_groups_masterfile.csv")
        master = merge_all(annotated)
    elif args.masterfile:
        print(f"[load] reading masterfile: {args.masterfile}")
        master = pd.read_csv(args.masterfile, low_memory=False)
    else:
        mf = out_dir / 'protein_groups_masterfile.csv'
        print(f"[load] reading masterfile from output dir: {mf}")
        master = pd.read_csv(mf, low_memory=False)

    ibaq_cols = ibaq_cols_of(master)
    print(f"         {len(master)} protein groups, {len(ibaq_cols)} samples")

    # ── Stage 3.25: re-annotate stale 'auto:' rows from master annotation ────
    # Applies when --masterfile is used: old auto:unmatched notes may persist
    # for proteins that ARE in the master annotation (just not exact-matched
    # in the previous run).  Re-running Tier 0 corrects Family + clears Notes.
    if args.master_annotation and 'Notes' in master.columns:
        stale_count = master['Notes'].fillna('').astype(str).str.startswith('auto:').sum()
        if stale_count > 0:
            print(f"[reannotate] {stale_count} rows have stale auto: notes — "
                  f"re-checking against master annotation ...")
            # Load master annotation if not already loaded (i.e. --masterfile path)
            if 'master_annot' not in dir() or master_annot is None:
                master_annot, master_accession = load_master_annotation(
                    Path(args.master_annotation))
            master, n_exact, n_acc = reannotate_from_master(
                master, master_annot, master_accession)
            print(f"[reannotate] fixed {n_exact} via exact match, "
                  f"{n_acc} via accession fallback")
            still_stale = master['Notes'].fillna('').astype(str).str.startswith('auto:').sum()
            print(f"[reannotate] {still_stale} rows genuinely absent from master annotation")

    # ── Stage 3.5: auto-classify proteins missing a Family ───────────────────
    if not args.skip_classify:
        master, cls_counts = auto_classify_masterfile(
            master, scorer, min_score=args.classify_min_score
        )
        total_new = sum(v for k, v in cls_counts.items()
                        if k != 'already_annotated')
        print(f"[classify] results:")
        print(f"           already annotated (Tier 0 exact match): "
              f"{cls_counts['already_annotated']}")
        print(f"           Tier 1 venom rules:    {cls_counts.get('tier1_venom', 0)}")
        print(f"           Tier 1 non-venom:       {cls_counts.get('tier1_nonvenom', 0)}")
        if scorer:
            print(f"           Tier 2 keyword venom:  {cls_counts.get('tier2_venom', 0)}")
            print(f"           Tier 2 keyword other:  {cls_counts.get('tier2_nonvenom', 0)}")
        print(f"           Unmatched (→ Other):    {cls_counts.get('unmatched', 0)}")
        print(f"           Total newly classified: {total_new}")
    else:
        print("[classify] skipped")

    # Normalise legacy / abbreviated family names to current canonical names
    # (e.g. Procoaculant → Prothrombin activator, Kaz → Kazal).
    # Runs after classify so Tier 0 exact-match rows are also corrected.
    master = normalise_family_names(master)

    # Save AllProteinGroups_anotated.csv — annotation catalogue for proteins
    # found in these runs, with all three classification tiers applied.
    # (Saved here rather than at Stage 2 so auto-classifications are included.)
    if not args.skip_annotate or args.masterfile:
        annot_present = [c for c in ANNOT_COLS if c in master.columns]
        cat = master[['Fasta headers'] + annot_present].drop_duplicates(
            subset='Fasta headers')
        cat_path = out_dir / 'AllProteinGroups_anotated.csv'
        cat.to_csv(cat_path, index=False)
        print(f"[annotate] saved catalogue: {cat_path} "
              f"({len(cat)} unique proteins, all tiers applied)")

    # Save masterfile (with all classifications applied).
    # Always write to out_dir so Tier 1/2 classifications are not lost,
    # even when --masterfile was supplied (the source file is never overwritten
    # unless out_dir happens to be the same directory).
    master_path = out_dir / 'protein_groups_masterfile.csv'
    if args.masterfile and Path(args.masterfile).resolve() == master_path.resolve():
        # Same path: only write back if we actually changed something
        changed = not args.skip_classify
        if changed:
            master.to_csv(master_path, index=False)
            print(f"[classify] updated masterfile in-place → {master_path}")
        else:
            print(f"[classify] no changes — masterfile unchanged: {master_path}")
    else:
        master.to_csv(master_path, index=False)
        print(f"[merge/classify] saved masterfile → {master_path}")

    # ── Load collection metadata ──────────────────────────────────────────────
    print(f"[collection] loading {args.collection}")
    collection = load_collection(Path(args.collection))
    print(f"             {len(collection)} samples in MASTER sheet")

    # ── Stage 4: MassSpecData.csv ─────────────────────────────────────────────
    if not args.skip_massspecdata:
        print("[massspecdata] building MassSpecData.csv ...")
        msd = build_massspecdata(master, collection)
        msd_path = out_dir / 'MassSpecData.csv'
        msd.to_csv(msd_path, index=False, header=False)
        print(f"[massspecdata] saved → {msd_path} ({msd.shape[0]} × {msd.shape[1]})")

    # ── Stage 5: proportions (both all-proteins and venom-only) ──────────────
    need_props = (not args.skip_summary or
                  not args.skip_averaged or
                  not args.skip_plots)
    prop_all, families_all     = None, []
    prop_venom, families_venom = None, []
    venom_pct: dict[str, float] = {}

    if need_props:
        no_pla2 = args.no_pla2_inhibitor
        print("[proportions] computing family iBAQ % — all proteins ...")
        prop_all, families_all = compute_proportions(
            master, venom_only=False, no_pla2_inhibitor=no_pla2)
        print(f"              {len(prop_all)} samples, {len(families_all)} families")

        print("[proportions] computing family iBAQ % — venom-only ...")
        prop_venom, families_venom = compute_proportions(
            master, venom_only=True, no_pla2_inhibitor=no_pla2)
        print(f"              {len(prop_venom)} samples, {len(families_venom)} families")

        print("[proportions] computing venom fraction per sample ...")
        venom_pct = compute_venom_percent(master, no_pla2_inhibitor=no_pla2)

    # ── Stage 6: per-sample summary ───────────────────────────────────────────
    # Two variants are written:
    #   SummaryResults_VenomOnly.csv   — family columns = % of total VENOM iBAQ
    #                                    (i.e. relative toxin-family composition)
    #   SummaryResults_AllProteins.csv — family columns = % of TOTAL iBAQ
    #                                    (venom + non-venom; shows full sample
    #                                    protein landscape incl. Other)
    # Both include the run-level metadata (sample type, run time, filter, etc.)
    # and Venom_pct so the two scales can always be reconciled.
    # The legacy filename SummaryResults.csv is kept as an alias for VenomOnly
    # so existing downstream scripts are not broken.
    summary = None
    summary_all_proteins = None
    if not args.skip_summary:
        print("[summary] building per-sample summary table — venom-only proportions ...")
        summary = build_summary_table(prop_venom, families_venom, collection,
                                      venom_pct=venom_pct,
                                      sample_universe_df=prop_all)
        path_venom = out_dir / 'SummaryResults_VenomOnly.csv'
        summary.to_csv(path_venom, index=False)
        # Legacy alias — keeps any existing downstream scripts working
        (out_dir / 'SummaryResults.csv').write_text(path_venom.read_text())
        print(f"[summary] saved → {path_venom}")

        print("[summary] building per-sample summary table — all-protein proportions ...")
        summary_all_proteins = build_summary_table(prop_all, families_all, collection,
                                                   venom_pct=venom_pct,
                                                   sample_universe_df=prop_all)
        path_all = out_dir / 'SummaryResults_AllProteins.csv'
        summary_all_proteins.to_csv(path_all, index=False)
        print(f"[summary] saved → {path_all}")

    # ── Stage 7: averaged summary ─────────────────────────────────────────────
    species_avg = None
    if not args.skip_averaged:
        if summary is None:
            summary = build_summary_table(prop_venom, families_venom, collection,
                                          venom_pct=venom_pct,
                                          sample_universe_df=prop_all)
        if summary_all_proteins is None:
            summary_all_proteins = build_summary_table(prop_all, families_all, collection,
                                                       venom_pct=venom_pct,
                                                       sample_universe_df=prop_all)
        print("[averaged] averaging by species + form — venom-only ...")
        species_avg = build_averaged_summary(summary, families_venom)
        path = out_dir / 'SummaryResults_VenomOnly_AveragedBySpeciesAndSampleType.csv'
        species_avg.to_csv(path, index=False)
        # Legacy alias
        (out_dir / 'SummaryResults_AveragedBySpeciesAndSampleType.csv').write_text(
            path.read_text())
        print(f"[averaged] saved → {path}  ({len(species_avg)} groups)")

        print("[averaged] averaging by species + form — all proteins ...")
        species_avg_all = build_averaged_summary(summary_all_proteins, families_all)
        path_all = out_dir / 'SummaryResults_AllProteins_AveragedBySpeciesAndSampleType.csv'
        species_avg_all.to_csv(path_all, index=False)
        print(f"[averaged] saved → {path_all}  ({len(species_avg_all)} groups)")

    # ── Stage 8: pie charts (two sets by default) ─────────────────────────────
    if not args.skip_plots:
        print("[pies] rendering pie charts — all proteins → PieCharts/ ...")
        n_all = generate_all_pies(
            prop_all, families_all, collection, out_dir,
            subfolder='PieCharts',
            filtered=args.filtered_pies,
        )
        print(f"[pies] {n_all} all-proteins pie(s) written")

        print("[pies] rendering pie charts — venom only → PieCharts_VenomOnly/ ...")
        n_venom = generate_all_pies(
            prop_venom, families_venom, collection, out_dir,
            subfolder='PieCharts_VenomOnly',
            filtered=args.filtered_pies,
            species_avg_df=species_avg if not args.skip_species_pies else None,
        )
        print(f"[pies] {n_venom} venom-only pie(s) written")

    # ── Stage 9: unclassified report ─────────────────────────────────────────
    if not args.skip_unclassified:
        path = out_dir / 'unclassified_for_review.csv'
        n = flag_unclassified(master, path)
        if n:
            print(f"[unclassified] {n} proteins with no Family → {path}")
        else:
            print("[unclassified] all proteins have a Family assigned.")

    print("\nDONE.")


if __name__ == '__main__':
    main()
