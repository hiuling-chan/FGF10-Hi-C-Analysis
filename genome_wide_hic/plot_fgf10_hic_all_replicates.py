#!/usr/bin/env python3
"""Generate chr5 genome-wide Hi-C and Capture Hi-C figures (GRCh38/hg38).

- Close view: chr5:42,796,905-44,441,081; balanced Hi-C at 20-kb resolution.
- Wide view: chr5:42,300,000-46,000,000; balanced Hi-C at 50-kb resolution.
- Insulation: 20-kb input, 300-kb window.
- Capture Hi-C arcs: all finite-score cis <=2-Mb pairs from the supplied BEDPE
  files whose BOTH anchors overlap the plotted interval; not bait-restricted.
  Blue: score <5 (including <3); orange: 5 to <7; red: >=7.
- GENCODE v50 gene spans and public CTCF tracks.

Required inputs are NOT distributed with this code. This script does not run
alignment, Cooler balancing, CHiCAGO, or insulation calling; it visualizes their
precomputed outputs. The original project-specific input names are documented
below so a user can supply them without exposing laboratory storage paths.

Expected project layout:
  PROJECT_ROOT/hic_results/contact_maps/cool/IMR-1.20000_balanced.cool
  PROJECT_ROOT/hic_results/contact_maps/cool/IMR-1.50000_balanced.cool
  (analogous files for IMR-2, Lymph-1, and Lymph-2)
  PROJECT_ROOT/hic_results/tads/**/{sample}.20000_balanced_insulation.tsv
  PROJECT_ROOT/annotations/gencode.v50.basic.annotation.gtf.gz

Expected Capture Hi-C layout:
  CAPTURE_ROOT/{IMR1,IMR2,Lymph1,Lymph2}/chicago/data/
    {sample}.cis_le_2Mb.bedpe

Examples:
  python plot_fgf10_hic_all_replicates.py --project-root PROJECT_ROOT --capture-root CAPTURE_ROOT --png
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import shutil
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

# ---------------- Figure and public-reference defaults ----------------
GENCODE_FILENAME = 'gencode.v50.basic.annotation.gtf.gz'
GENCODE_URL = (
    'https://ftp.ebi.ac.uk/pub/databases/gencode/'
    'Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz'
)
CHROM = 'chr5'
SNP_POS_1_BASED = 43_619_485
SAMPLES = {
    'IMR1': ('IMR-1', 'IMR-90 replicate 1'),
    'IMR2': ('IMR-2', 'IMR-90 replicate 2'),
    'Lymph1': ('Lymph-1', 'Lymphoblast replicate 1'),
    'Lymph2': ('Lymph-2', 'Lymphoblast replicate 2'),
}
# gene labels to keep on the plot
GENES = ('ANXA2R', 'ZNF131', 'NIM1K', 'HMGCS1', 'CCL28',
         'TMEM267', 'PAIP1', 'NNT', 'FGF10')
# ENCODE CTCF tracks
ENCODE_TRACKS = (
    ('ENCODE CTCF\nlung fibroblast', 'ENCFF777ODE', ('bigBed', 'bb', 'bigWig', 'bw')),
    ('ENCODE CTCF\nIMR-90', 'ENCFF307XFM', ('bigBed', 'bb', 'bigWig', 'bw')),
)

@dataclass(frozen=True)
class View:
    name: str
    start: int
    end: int
    resolution: int
    depth: int
    hic_height: float
    arc_height: float
    width_cm: float

VIEWS = (
    View('close', 42_796_905, 44_441_081, 20_000, 1_000_000, 7.0, 1.55, 40),
    View('wide', 42_300_000, 46_000_000, 50_000, 2_000_000, 6.7, 1.40, 42),
)

# ---------------- Utility helpers ----------------
def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def is_overlap(a: int, b: int, start: int, end: int) -> bool:
    return b > start and a < end


def chrom_normalized(s: str) -> str:
    return s if s.startswith('chr') else 'chr' + s


def safe_download(url: str, dest: Path, timeout: int = 60) -> bool:
    if dest.is_file() and dest.stat().st_size > 0:
        return True
    mkdir(dest.parent)
    tmp = dest.with_name(dest.name + '.download')
    try:
        with urllib.request.urlopen(url, timeout=timeout) as src, tmp.open('wb') as out:
            shutil.copyfileobj(src, out, length=1024*1024)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise IOError('downloaded file is empty')
        tmp.replace(dest)
        return True
    except Exception as e:
        tmp.unlink(missing_ok=True)
        print(f'  Download attempt failed ({url}): {e}', file=sys.stderr)
        return False


def stage_encode_ctcf(accession: str, label: str, cache_dir: Path, *, allow_download: bool = True):
    out_bed = cache_dir / f'{accession}.chr5.peaks.bed'
    if out_bed.exists() and out_bed.stat().st_size > 0:
        return out_bed
    if not allow_download:
        print(f'  ENCODE {accession}: not cached; skipping (--offline).', file=sys.stderr)
        return None
    urls = [
        (f'https://www.encodeproject.org/files/{accession}/@@download/{accession}.bed.gz', 'gz'),
        (f'https://www.encodeproject.org/files/{accession}/@@download/{accession}.bed', 'plain'),
        (f'https://www.encodeproject.org/files/{accession}/@@download/{accession}.narrowPeak.gz', 'gz'),
    ]
    for url, kind in urls:
        compressed = cache_dir / (accession + ('.bed.gz' if kind == 'gz' else '.bed'))
        if not safe_download(url, compressed):
            continue
        tmp = out_bed.with_name(out_bed.name + '.tmp')
        n_total = 0
        n_chr5 = 0
        try:
            opener = gzip.open if kind == 'gz' else open
            with opener(compressed, 'rt') as src, tmp.open('w') as dest:
                for line in src:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3 or parts[0].startswith(('#', 'track', 'browser')):
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    n_total += 1
                    if chrom_normalized(parts[0]) != CHROM or not is_overlap(start, end, VIEWS[-1].start, VIEWS[-1].end):
                        continue
                    name = parts[3] if len(parts) > 3 else label.replace(' ', '_')
                    score = parts[4] if len(parts) > 4 else '0'
                    dest.write(f'{CHROM}\t{start}\t{end}\t{name}\t{score}\t.\n')
                    n_chr5 += 1
            if n_total == 0:
                raise ValueError('not a BED file with intervals')
            tmp.replace(out_bed)
            print(f'  ENCODE {accession}: extracted {n_chr5} chr5 peaks from {n_total} total peaks')
            return out_bed
        except Exception as e:
            print(f'  Could not parse {compressed} as BED: {e}', file=sys.stderr)
            tmp.unlink(missing_ok=True)
            compressed.unlink(missing_ok=True)
    print(f'WARNING: ENCODE {accession} unavailable; omitting {label} track.', file=sys.stderr)
    return None


def stage_external_tracks(cache_dir: Path, *, allow_download: bool = True) -> dict[str, tuple[Path, str, str]]:
    staged: dict[str, tuple[Path, str, str]] = {}
    mkdir(cache_dir)
    for label, url, file_type, color in ATAC_TRACKS:
        local = cache_dir / Path(url).name
        if local.is_file() and local.stat().st_size > 0:
            staged[label] = (local, file_type, color)
        elif allow_download and safe_download(url, local):
            staged[label] = (local, file_type, color)
        else:
            print(f'  Optional ATAC-seq track unavailable: {label!r}; skipping.', file=sys.stderr)
    for label, accession, _exts in ENCODE_TRACKS:
        peaks = stage_encode_ctcf(accession, label, cache_dir, allow_download=allow_download)
        if peaks:
            staged[label] = (peaks, 'bed', '#444444')
    return staged


def download_annotation(gtf: Path, *, allow_download: bool = True) -> None:
    if gtf.is_file() and gtf.stat().st_size > 0:
        return
    if not allow_download:
        raise FileNotFoundError(f'Missing GENCODE GTF: {gtf}. Supply --gtf or disable --offline.')
    mkdir(gtf.parent)
    tmp = gtf.with_name(gtf.name + '.download')
    print(f'Downloading GENCODE v50 annotation to {gtf} ...')
    try:
        urllib.request.urlretrieve(GENCODE_URL, tmp)
        with gzip.open(tmp, 'rb') as h:
            if not h.readline():
                raise ValueError('Annotation download was empty')
        tmp.replace(gtf)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def generate_gene_bed(gtf: Path, out: Path) -> int:
    """Generate simple BED6 gene-span track (not exon-level models).

    This makes the gene track visually cleaner than full transcripts.
    """
    wanted = set(GENES)
    parsed = {}
    opener = gzip.open if gtf.suffix == '.gz' else open
    attr = re.compile(r'(?:^|;)\s*gene_name\s+"([^"]+)"')
    with opener(gtf, 'rt') as src:
        for raw in src:
            if raw.startswith('#'):
                continue
            parts = raw.rstrip('\n').split('\t')
            if len(parts) != 9 or chrom_normalized(parts[0]) != CHROM or parts[2] != 'gene':
                continue
            m = attr.search(parts[8])
            if not m or m.group(1) not in wanted:
                continue
            name = m.group(1)
            start, end = int(parts[3]) - 1, int(parts[4])
            if not is_overlap(start, end, VIEWS[-1].start, VIEWS[-1].end):
                continue
            if name in parsed:
                s, e, strand = parsed[name]
                parsed[name] = (min(start, s), max(end, e), strand)
            else:
                parsed[name] = (start, end, parts[6])
    if not parsed:
        raise RuntimeError(f'No named GENCODE gene records found in {gtf}')
    with out.open('w') as dest:
        for name, (start, end, strand) in sorted(parsed.items(), key=lambda i: i[1][0]):
            dest.write(f'{CHROM}\t{start}\t{end}\t{name}\t0\t{strand}\n')
    return len(parsed)


def generate_snp_bed(out: Path) -> Path:
    marker = out / 'candidate_modifier_snp.bed'
    marker.write_text(f'{CHROM}\t{SNP_POS_1_BASED-1}\t{SNP_POS_1_BASED}\tSNP\n')
    return marker


def find_insulation(project: Path, cool_sample: str) -> Path:
    name = f'{cool_sample}.20000_balanced_insulation.tsv'
    base = project / 'hic_results' / 'tads'
    hits = sorted(base.rglob(name))
    if not hits:
        raise FileNotFoundError(f'Missing insulation file {name}; searched under {base}')
    return hits[0]


def generate_insulation_bedgraph(tsv: Path, output: Path, prefix: str, window_bp: int = 300_000) -> Path:
    bedgraph = output / f'{prefix}_insulation_{window_bp//1000}kb.bedgraph'
    value_key = f'log2_insulation_score_{window_bp}'
    count = 0
    with tsv.open() as f, bedgraph.open('w') as bg:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or value_key not in reader.fieldnames:
            raise ValueError(f'Expected insulation field {value_key} not found in {tsv}')
        for row in reader:
            if chrom_normalized(row['chrom']) != CHROM:
                continue
            start, end = int(row['start']), int(row['end'])
            if not is_overlap(start, end, VIEWS[-1].start, VIEWS[-1].end):
                continue
            try:
                value = float(row[value_key])
            except (ValueError, TypeError):
                continue
            if math.isfinite(value):
                bg.write(f'{CHROM}\t{start}\t{end}\t{value}\n')
                count += 1
    if count == 0:
        raise RuntimeError(f'No finite chr5 insulation values found in {tsv}')
    return bedgraph


def filter_capture_bedpe(source: Path, out_prefix: Path, view: View) -> tuple[dict[str, Path], dict[str, int]]:
    """Keep all cis<=2 Mb interactions overlapping the view and split into 3 score bins.

    Creates three BEDPE files suitable for separate arc tracks overlaid in the same panel:
      - scorelt5: score < 5
      - score5to7: 5 <= score < 7
      - score7plus: score >= 7

    Returns ({bin_name: path}, {bin_name: count}).
    """
    if not source.is_file():
        raise FileNotFoundError(f'Capture BEDPE not found: {source}')

    bins = {
        'scorelt5': out_prefix.with_name(out_prefix.name + '_scorelt5.bedpe'),
        'score5to7': out_prefix.with_name(out_prefix.name + '_score5to7.bedpe'),
        'score7plus': out_prefix.with_name(out_prefix.name + '_score7plus.bedpe'),
    }
    counts = {k: 0 for k in bins}
    handles = {k: v.open('w') for k, v in bins.items()}
    try:
        with source.open() as inp:
            for line in inp:
                line = line.strip()
                if not line or line.startswith(('#', 'track', 'browser')):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    c1, start1, end1 = chrom_normalized(parts[0]), int(parts[1]), int(parts[2])
                    c2, start2, end2 = chrom_normalized(parts[3]), int(parts[4]), int(parts[5])
                except ValueError:
                    continue
                if c1 != CHROM or c2 != CHROM:
                    continue
                if not (is_overlap(start1, end1, view.start, view.end) and is_overlap(start2, end2, view.start, view.end)):
                    continue
                if len(parts) < 7:
                    raise ValueError(
                        f'BEDPE record lacks a CHiCAGO score (seventh column): {source}: {line}'
                    )
                try:
                    score = float(parts[6])
                except ValueError as exc:
                    raise ValueError(f'Non-numeric CHiCAGO score in {source}: {line}') from exc
                if not math.isfinite(score):
                    raise ValueError(f'Non-finite CHiCAGO score in {source}: {line}')
                row = f'{CHROM}\t{start1}\t{end1}\t{CHROM}\t{start2}\t{end2}\t{score}\n'
                if score >= 7:
                    key = 'score7plus'
                elif score >= 5:
                    key = 'score5to7'
                else:
                    key = 'scorelt5'
                handles[key].write(row)
                counts[key] += 1
    finally:
        for fh in handles.values():
            fh.close()
    return bins, counts


def make_ini(ini_path: Path, *, sample_key: str, sample_label: str, view: View,
             cool: Path, arcs_by_bin: dict[str, Path], genes: Path, ins_bg: Path,
             snp_bed: Path, external_tracks: dict[str, tuple[Path, str, str]]) -> None:
    sample_short = SAMPLES[sample_key][0]
    is_imr = sample_key.startswith('IMR')

    def title(*lines: str) -> str:
        # Indented continuation lines: valid INI, parsed as one title value.
        return 'title = ' + lines[0] + ''.join('\n    ' + line for line in lines[1:])

    sections: list[str] = []
    sections.append(f'''[hic matrix]
file = {cool}
file_type = hic_matrix
{title(sample_short, f'Genome-wide Hi-C ({view.resolution//1000} kb)')}
height = {view.hic_height}
depth = {view.depth}
transform = log1p
colormap = Reds
show_masked_bins = false
''')
    arc_colors = {
        'scorelt5': '#4477CC',  # blue: all scores < 5
        'score5to7': '#E6A044',  # orange: 5 <= score < 7
        'score7plus': '#D74A4A', # red: score >= 7
    }
    ordered_bins = ('scorelt5', 'score5to7', 'score7plus')
    active_bins = [key for key in ordered_bins
                   if (key in arcs_by_bin and arcs_by_bin[key].exists()
                       and arcs_by_bin[key].stat().st_size > 0)]
    for i, key in enumerate(active_bins):
        first = (i == 0)
        arc_path = arcs_by_bin[key]
        title_line = (title(sample_label, 'Capture Hi-C (cis <=2 Mb)')
                      if first else 'title = ')
        overlay_line = '' if first else 'overlay_previous = share-y\n'
        height_line = f'height = {view.arc_height}\n' if first else ''
        sections.append(f'''[capture_interactions_{key}]
file = {arc_path}
file_type = links
{title_line}
links_type = arcs
use_middle = true
compact_arcs_level = 2
line_width = 0.60
alpha = 0.80
color = {arc_colors[key]}
{height_line}{overlay_line}''')
    sections.append(f'''[insulation]
file = {ins_bg}
file_type = bedgraph
{title('300 kb insulation', 'window')}
type = line:1.20
height = 1.45
color = steelblue
show_data_range = false
''')
    sections.append(f'''[gene annotations]
file = {genes}
file_type = bed
{title('GENCODE v50', '(hg38)')}
color = navy
border_color = navy
style = UCSC
labels = true
fontsize = 7
height = 1.05
max_labels = 20
''')
    ordered = (
        'ENCODE CTCF\nlung fibroblast', 'ENCODE CTCF\nIMR-90')
    for i, label in enumerate(ordered):
        item = external_tracks.get(label)
        if not item:
            continue
        path, file_type, color = item
        short = label.split('\n')
        if file_type == 'bigwig':
            sections.append(f'''[external_bigwig_{i}]
file = {path}
file_type = bigwig
{title(*short)}
height = 0.80
color = {color}
type = fill
min_value = 0
nans_to_zeros = true
number_of_bins = 800
show_data_range = false
''')
        else:
            sections.append(f'''[external_bed_{i}]
file = {path}
file_type = bed
{title(*short)}
height = 0.58
color = #4d4d4d
border_color = #4d4d4d
labels = false
display = collapsed
''')
    if is_imr:
        sections.append(f'''[SNP position]
file = {snp_bed}
type = vlines
color = black
alpha = 0.95
line_style = dotted
line_width = 0.9
zorder = 20
''')
    sections.append('[x-axis]\nfontsize = 7\n')
    ini_path.write_text('\n'.join(sections))
    import configparser
    test = configparser.ConfigParser(interpolation=None, strict=True)
    with ini_path.open() as fh:
        test.read_file(fh)


def render(ini: Path, pdf_path: Path, view: View, png: bool) -> None:
    try:
        from pygenometracks.tracksClass import PlotTracks
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise RuntimeError('pyGenomeTracks is not available in this Python environment. '
                           'Activate your pygenometracks Conda environment and rerun.') from e

    plotter = PlotTracks(
        tracks_file=str(ini),
        fig_width=view.width_cm,
        fontsize=7,
        dpi=220,
        track_label_width=0.13,
        plot_regions=[(CHROM, view.start, view.end)],
    )
    # left axis / main panel / right-side titles
    plotter.width_ratios = (0.055, 0.855, 0.09)
    fig = plotter.plot(str(pdf_path), CHROM, view.start, view.end)
    import configparser
    from matplotlib.lines import Line2D
    cfg = configparser.ConfigParser(interpolation=None)
    cfg.read(str(ini))
    has_arcs = any(section.startswith('capture_interactions_') for section in cfg.sections())
    if has_arcs:
        handles = [
            Line2D([0], [0], color='#4477CC', lw=1.2, label='Score <5'),
            Line2D([0], [0], color='#E6A141', lw=1.2, label='Score 5 to <7'),
            Line2D([0], [0], color='#D62728', lw=1.2, label='Score >=7'),
        ]
        main_axes = [ax for ax in fig.axes if ax.get_position().width > 0.35
                     and abs(ax.get_xlim()[0] - view.start) < 1
                     and abs(ax.get_xlim()[1] - view.end) < 1]
        if len(main_axes) >= 2:
            main_axes[1].legend(
                handles=handles, loc='upper right', ncol=3, fontsize=6,
                frameon=True, framealpha=0.85, facecolor='white',
                edgecolor='none', handlelength=1.7, columnspacing=1.0,
                handletextpad=0.35, borderpad=0.2, labelspacing=0.2,
                bbox_to_anchor=(0.995, 0.995),
            )
        else:
            fig.legend(handles=handles, loc='lower center', ncol=3,
                       fontsize=6, frameon=False)

    fig.savefig(pdf_path, bbox_inches='tight', pad_inches=0.22, dpi=220)
    if png:
        fig.savefig(pdf_path.with_suffix('.png'), bbox_inches='tight', pad_inches=0.22, dpi=220)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--project-root', type=Path, required=True,
                   help='Project directory containing hic_results and annotations')
    p.add_argument('--capture-root', type=Path, required=True,
                   help='Directory containing sample-specific CHiCAGO BEDPE outputs')
    p.add_argument('--outdir', type=Path, default=None,
                   help='Default: <project-root>/figures/chr5_all_samples_v7')
    p.add_argument('--gtf', type=Path, default=None,
                   help='GENCODE v50 basic GTF or GTF.gz; downloads if absent')
    p.add_argument('--samples', nargs='+', choices=list(SAMPLES), default=list(SAMPLES))
    p.add_argument('--png', action='store_true', help='Save PNG previews in addition to PDFs')
    p.add_argument('--prepare-only', action='store_true',
                   help='Generate configuration and filtered tracks without plotting')
    p.add_argument('--annotation-cache', type=Path, default=None,
                   help='Persistent public-annotation cache (default: PROJECT_ROOT/annotations/external_tracks)')
    p.add_argument('--offline', action='store_true',
                   help='Never download GENCODE, ENCODE, or ATAC-seq files; use cached/local files')
    p.add_argument('--skip-external-tracks', action='store_true',
                   help='Omit optional ENCODE CTCF and ATAC-seq tracks')
    return p.parse_args()



def main() -> int:
    args = parse_args()
    project = args.project_root.expanduser().resolve()
    caproot = args.capture_root.expanduser().resolve()
    if not project.is_dir():
        raise FileNotFoundError(f'Project root does not exist: {project}')
    if not caproot.is_dir():
        raise FileNotFoundError(f'Capture Hi-C root does not exist: {caproot}')
    output = mkdir((args.outdir or (project / 'figures' / 'chr5_all_samples_v7')).expanduser().resolve())
    gtf = (args.gtf or project / 'annotations' / GENCODE_FILENAME).expanduser().resolve()
    cache_dir = (args.annotation_cache or project / 'annotations' / 'external_tracks').expanduser().resolve()

    download_annotation(gtf, allow_download=not args.offline)
    staged_tracks = {}
    if not args.skip_external_tracks:
        staged_tracks = stage_external_tracks(mkdir(cache_dir), allow_download=not args.offline)

    genes = output / 'chr5_GENCODE_v50_selected_gene_spans.bed'
    n_genes = generate_gene_bed(gtf, genes)
    snp_bed = generate_snp_bed(output)
    print(f'GENCODE: {n_genes} selected gene spans')
    print(f'SNP marker position (1-based): {CHROM}:{SNP_POS_1_BASED:,}')
    if staged_tracks:
        print('External tracks staged:')
        for k, (v, ft, _c) in staged_tracks.items():
            print(f'  - {k}: {v.name} [{ft}]')

    nfig = 0
    failed = []
    for sample_key in args.samples:
        cool_id, label = SAMPLES[sample_key]
        print(f'\n=== {sample_key} / {cool_id} ===')
        sample_out = mkdir(output / sample_key)
        cis_bedpe = caproot / sample_key / 'chicago' / 'data' / f'{sample_key}.cis_le_2Mb.bedpe'
        try:
            ins_tsv = find_insulation(project, cool_id)
            ins_bg = generate_insulation_bedgraph(ins_tsv, sample_out, sample_key, window_bp=300_000)
            for view in VIEWS:
                cool = project / 'hic_results' / 'contact_maps' / 'cool' / f'{cool_id}.{view.resolution}_balanced.cool'
                if not cool.is_file():
                    raise FileNotFoundError(f'Missing balanced Hi-C matrix: {cool}')
                arc_prefix = sample_out / f'{sample_key}_chr5_{view.name}_cis_le_2Mb'
                arcs_by_bin, counts = filter_capture_bedpe(cis_bedpe, arc_prefix, view)
                narc = sum(counts.values())
                print(f'    {view.name}: {narc} arcs in region '
                      f'(<5: {counts["scorelt5"]}, '
                      f'5-<7: {counts["score5to7"]}, >=7: {counts["score7plus"]})')
                if narc == 0:
                    print(f'    WARNING: no arcs within {view.name} '
                          f'view in {cis_bedpe}', file=sys.stderr)
                ini = sample_out / f'{sample_key}_chr5_{view.name}.ini'
                pdf = sample_out / f'{sample_key}_chr5_{view.name}.pdf'
                make_ini(ini, sample_key=sample_key, sample_label=label, view=view,
                         cool=cool, arcs_by_bin=arcs_by_bin, genes=genes,
                         ins_bg=ins_bg, snp_bed=snp_bed,
                         external_tracks=staged_tracks)
                if args.prepare_only:
                    print(f'    Prepared: {ini}')
                else:
                    render(ini, pdf, view, args.png)
                    print(f'    Saved: {pdf}')
                    nfig += 1
        except Exception as exc:
            print(f'ERROR for {sample_key}: {exc}', file=sys.stderr)
            failed.append(sample_key)

    (output / 'README.txt').write_text(
        'chr5 multi-sample Hi-C + Capture Hi-C figure set (hg38)\n\n'
        f'Close region: {CHROM}:{VIEWS[0].start}-{VIEWS[0].end}; 20 kb genome-wide Hi-C bins\n'
        f'Wide region: {CHROM}:{VIEWS[1].start}-{VIEWS[1].end}; 50 kb genome-wide Hi-C bins\n'
        'Capture Hi-C arcs come from sample-specific *.cis_le_2Mb.bedpe files and are\n'
        'colored in three bins: blue <5, orange 5-<7, red >=7.\n'
        'All finite-score arcs are plotted, including scores under 3.\n'
        'This is not restricted to a high-confidence subset.\n\n'
        f'SNP marker: {CHROM}:{SNP_POS_1_BASED} (1-based), shown only on IMR1/IMR2\n'
        'as a dotted black vertical line, without a separate label track.\n\n'
        'Gene track: GENCODE v50 (hg38) simplified selected gene spans.\n'
        'Optional public external tracks are added when available; missing ones\n'
        'are skipped with warnings so the script can still complete.\n\n'
        'Inputs are not distributed with this code.\n'
    )

    print(f'\nOutput directory: {output}')
    if failed:
        print('Not completed: ' + ', '.join(failed), file=sys.stderr)
        return 1
    print('Finished: ' + ('prepared configurations only' if args.prepare_only
                         else f'{nfig} figures ({len(args.samples)} samples x 2 views)'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
