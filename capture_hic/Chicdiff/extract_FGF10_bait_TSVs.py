#!/usr/bin/env python3
"""Export reproducible FGF10 promoter/proximal Chicdiff results.

Requires Python 3 standard library only. Does NOT rerun statistics.

Default inputs:
  file.path/IMR_vs_Lymph_full_results.tsv
  file.path/fgf10_anchor_audit/
      06_all_candidate_fragment_chicago_scores_counts.tsv

The exact-pair source is optional: it is original per-replicate CHiCAGO
counts/scores, not Chicdiff region-level significance.
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path('file.path/chicdiff_results/differential')
DEFAULT_EXACT = Path('file.path/annotations/fgf10_anchor_audit/06_all_candidate_fragment_chicago_scores_counts.tsv')
SNP_POS = 43619485  # rs78368856, hg38, one-based genomic position
SNP_FRAGMENT = 282709
PROMOTER = {282949, 282950}
BAITS = {282949, 282950, 282953}
BAIT_INFO = [
    ('282949', 'FGF10 promoter-overlapping bait A', 'chr5', '44384622', '44388806',
     'Overlaps assayed FGF10 promoter interval chr5:44,388,434-44,390,361 at its left edge'),
    ('282950', 'FGF10 promoter-overlapping bait B', 'chr5', '44388807', '44392024',
     'Overlaps assayed FGF10 promoter interval chr5:44,388,434-44,390,361'),
    ('282953', 'FGF10 promoter-proximal bait', 'chr5', '44398564', '44401751',
     'Located downstream of the specified promoter interval; not promoter-overlapping'),
]


def read_tsv(path):
    with path.open('r', newline='', encoding='utf-8-sig') as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        if not rd.fieldnames:
            raise ValueError(f'Empty or headerless TSV: {path}')
        return list(rd), list(rd.fieldnames)


def write_tsv(path, headers, rows):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f'  {path.name}: {len(rows)} rows')


def as_int(value):
    return int(value.strip())


def overlaps(row):
    # SNP is interior to the source fragment, so BED/1-based edge conventions
    # do not affect this specific SNP-overlap classification.
    return (row['OEchr'] == 'chr5'
            and as_int(row['OEstart']) <= SNP_POS
            and as_int(row['OEend']) >= SNP_POS)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', type=Path, default=ROOT/'IMR_vs_Lymph_full_results.tsv')
    ap.add_argument('--exact', type=Path, default=DEFAULT_EXACT,
                    help='Original per-replicate CHiCAGO score/count TSV; optional')
    ap.add_argument('--outdir', type=Path, default=ROOT/'FGF10_tables')
    args = ap.parse_args()
    if not args.results.is_file():
        ap.error(f'Chicdiff results not found: {args.results}')
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows, cols = read_tsv(args.results)
    required = {'baitID', 'OEchr', 'OEstart', 'OEend', 'regionID',
                'log2FoldChange', 'pvalue', 'padj', 'weighted_padj'}
    if not required.issubset(cols):
        ap.error('Missing Chicdiff columns: ' + ', '.join(sorted(required-set(cols))))

    selected = [r for r in rows if as_int(r['baitID']) in BAITS]
    selected.sort(key=lambda r: (as_int(r['baitID']), as_int(r['OEstart']), as_int(r['regionID'])))
    promoter = [r for r in selected if as_int(r['baitID']) in PROMOTER]
    proximal = [r for r in selected if as_int(r['baitID']) == 282953]
    promoter_snp = [r for r in promoter if overlaps(r)]
    proximal_snp = [r for r in proximal if overlaps(r)]
    all_snp = [r for r in selected if overlaps(r)]

    print('Exporting FGF10 Chicdiff TSV files:')
    outputs = (
        ('FGF10_all_3_baits_all.tsv', selected),
        ('FGF10_promoter_baits_all.tsv', promoter),
        ('FGF10_promoter_baits_SNP_overlap.tsv', promoter_snp),
        ('FGF10_proximal_bait282953_all.tsv', proximal),
        ('FGF10_proximal_bait282953_SNP_overlap.tsv', proximal_snp),
        ('FGF10_all_3_baits_SNP_overlap.tsv', all_snp),
    )
    for name, content in outputs:
        write_tsv(args.outdir/name, cols, content)

    annot_cols = ['baitID', 'annotation', 'chr', 'baitStart', 'baitEnd', 'description']
    annot_rows = [dict(zip(annot_cols, r)) for r in BAIT_INFO]
    annot_rows.append(dict(zip(annot_cols, (
        '282709', 'SNP-containing OTHER-END fragment (not a captured bait)',
        'chr5', '43616691', '43619914',
        'Contains rs78368856 at chr5:43,619,485 (GRCh38); compare to promoter and proximal baits'
    ))))
    write_tsv(args.outdir/'FGF10_bait_and_SNP_fragment_annotations.tsv', annot_cols, annot_rows)

    if args.exact.is_file():
        original, ocols = read_tsv(args.exact)
        necessary = {'sample','baitID','otherEndID','N','score'}
        if not necessary.issubset(ocols):
            ap.error('Exact-pair file lacks columns: ' + ', '.join(sorted(necessary-set(ocols))))
        match = [r for r in original if as_int(r['baitID']) in BAITS
                 and as_int(r['otherEndID']) == SNP_FRAGMENT]
        rank = {'IMR1': 0, 'IMR2': 1, 'Lymph1': 2, 'Lymph2': 3}
        match.sort(key=lambda r: (as_int(r['baitID']), rank.get(r['sample'], 9)))
        write_tsv(args.outdir/'FGF10_exact_SNP_fragment_pairs_counts_scores.tsv', ocols, match)
        if len(match) != 12:
            print(f'WARNING: expected 12 exact-pair rows (3 baits x 4 samples), found {len(match)}', file=sys.stderr)
    else:
        print(f'NOTE: Exact-pair source absent ({args.exact}); skipped sample-specific scores/counts.', file=sys.stderr)
    print('\nSNP-overlapping Chicdiff regions: ', len(all_snp), ' (overlapping tests are NOT independent loops)', sep='')
    print('Statistical contrasts: verify group direction before describing sign of log2FoldChange.')


if __name__ == '__main__':
    main()
