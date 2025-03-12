#!/usr/bin/env python
"""
Predict simple translocations

Usage:
    predict_traslocation.py [options] <chrom_sizes> <model> <seq_dir> <atac> <ctcf> <SV>...

Arguments:
    chrom_sizes  Path to tab-separated chromosome size file
    model        Path to model weights to instantiate Corigami model
    seq_dir      Path to directory with DNA sequences for Corigami model
    atac         Path to ATAC signal bigwig file
    ctcf         Path to CTCF signal bigwig file
    SV           String defining the translocation see README for details.

Options:
    -a CODE --assembly=CODE          Genome assembly code used to generate data [default: hg38]
    -o DIR --outdir=DIR              Directory to write results [default: ./]
    -r BINSIZE --resolution=BINSIZE  Bin size measured in base pairs [default: 8192]
    -w WINSIZE --window=WINSIZE      Size of window that Corigami predicted [default: 2097152]
    -n N --nsteps=N                  Number of steps per window used to shift the prediction region [default: 4]
    --atac_norm=FUN                  Function used to normalize ATAC signal [default: log]
    --ctcf_norm=FUN                  Function used to normalize CTCF signal
    --skip_cool                      Do not write results in cool format
    --skip_figures                   Do not plot figures of translocation
    --skip_csv                       Do not write results in csv format
"""
import argparse
import os
import re
from collections import namedtuple
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
import cooler

import torch
from torch.utils.data import Dataset, DataLoader

from corigami.inference.utils.model_utils import load_default
from corigami.model.corigami_models import ConvTransModel
from corigami.data.data_feature import GenomicFeature, SequenceFeature

import matplotlib.pyplot as plt


class Interval(namedtuple('Interval', ['chrom', 'start', 'end'])):
    __slots__ = ()

    def __repr__(self):
        strand = '+' if self[2] > self[1] else '-'
        return '{chrom}:{strand}:{start}-{end}'.format(chrom=self[0], strand=strand, start=self[1], end=self[2])

    def __len__(self):
        return abs(self.end - self.start)

    @property
    def strand(self):
        return '+' if self.end > self.start else '-'


class Translocation:
    def __init__(self, interval1, interval2):
        self.interval = (interval1, interval2)

    def __repr__(self):
        sep = ':=:' if self.interval[0].strand == self.interval[1].strand else ':x:'
        return str(self.interval[0]) + sep + str(self.interval[1])

    def __getitem__(self, idx):
        return self.interval[idx]

    def __iter__(self):
        return iter(self.interval)

    def __len__(self):
        return sum(len(i) for i in self)


class Features(namedtuple('Features', ['seq', 'atac', 'ctcf'])):
    __slots__ = ()

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return Features(self.seq[idx,:].copy(), self.atac[idx].copy(), self.ctcf[idx].copy())
        return super().__getitem__(idx)

    @property
    def DataFrame(self):
        data = np.concatenate([self.seq[:,:4], self.atac[:, np.newaxis], self.ctcf[:, np.newaxis]], axis=1)
        return pd.DataFrame(data, columns=['A', 'T', 'C', 'G', 'ATAC', 'CTCF'])

    @property
    def input(self):
        seq = torch.tensor(self.seq).unsqueeze(0)
        atac = torch.tensor(self.atac)
        ctcf = torch.tensor(np.nan_to_num(self.ctcf, 0))
        features = [ctcf, atac]
        features = torch.cat([feat.unsqueeze(0).unsqueeze(2) for feat in features], dim = 2)
        inputs = torch.cat([seq, features], dim = 2)
        return inputs

    @property
    def size(self):
        return len(self.atac)


class GenomeData:
    def __init__(self, seq_prefix, ctcf, atac, ctcf_norm=None, atac_norm='log', chromosomes=None):
        seq_prefix = Path(seq_prefix)
        if chromosomes is None:
            chromosomes = [path.name.removesuffix('.fa.gz') for path in seq_prefix.glob("*.fa.gz")]
        self.seq = {chrom: SequenceFeature(path=seq_prefix / f'{chrom}.fa.gz') for chrom in chromosomes}
        self.ctcf = GenomicFeature(path=ctcf, norm=ctcf_norm)
        self.atac = GenomicFeature(path=atac, norm=atac_norm)

    def get_dataset(self, region, window, step):
        # TODO: check that final length is fixed
        if isinstance(region, Interval):
            features = Features(*self._get_interval(region))
        elif isinstance(region, Translocation):
            features = zip(self._get_interval(region[0]), self._get_interval(region[1]))
            features = Features(*(np.concatenate(x, axis=0) for x in features))
        else:
            raise TypeError(f'Unknown type of input: {type(region)}')
        return ContigData(str(region), features, step, window)

    def _get_interval(self, interval):
        if interval.end > interval.start:
            seq = self.seq[interval.chrom].get(interval.start, interval.end)
            atac = self.atac.get(interval.chrom, interval.start, interval.end)
            ctcf = self.ctcf.get(interval.chrom, interval.start, interval.end)
        else:
            seq = self.seq[interval.chrom].get(interval.end, interval.start)
            atac = self.atac.get(interval.chrom, interval.end, interval.start)
            ctcf = self.ctcf.get(interval.chrom, interval.end, interval.start)
            seq = np.flip(seq, 0)[:, [1, 0, 3, 2, 4]]  # reverse & complement
            atac = np.flip(atac, 0)
            ctcf = np.flip(ctcf, 0)
        return seq, atac, ctcf


class ContigData(Dataset):
    def __init__(self, contig, features, step, window):
        self.contig = contig
        self.features = features
        self.step = step
        self.window = window

    def __getitem__(self, idx):
        s = self.step * idx
        t = s + self.window
        features = self.features[s:t]
        return Interval(self.contig, s, t), features.input

    def __len__(self):
        return 1 + (self.features.size - self.window) // self.step

    def __iter__(self):
        for k in range(len(self)):
            yield self[k]

    def make_bins(self, res):
        bins = pd.DataFrame({'chrom': self.contig, 'start': np.arange(0, self.size, res)})
        bins['end'] = np.minimum(bins['start'] + res, self.size)
        return bins

    def summarize_features(self, resolution):
        bins = self.make_bins(resolution)
        feat = self.features.DataFrame
        feat = feat.groupby(feat.index // resolution).mean()
        if 'G' in feat.columns and 'C' in feat.columns:
            feat['GC'] = feat['G'] + feat['C']
            feat.drop(columns=['G', 'C'], inplace=True)
        if 'A' in feat.columns and 'T' in feat.columns:
            feat['AT'] = feat['A'] + feat['T']
            feat.drop(columns=['A', 'T'], inplace=True)
        feat = pd.concat([bins, feat], axis=1)
        return feat

    @property
    def size(self):
        return self.features.size


def parse_alt(extended_alt):
    extended_alt = extended_alt.strip()

    # Build regex to match pattern
    pos = r'(chr[1-9XY][0-9]?):([0-9]+)'
    delim = r'([\[\]])?'
    regex = delim + pos + delim + pos + delim

    # Match pattern
    match = re.match(regex, extended_alt)
    if not match:
        raise ValueError(f'Cannot parse ALT: {extended_alt}')
    delim1, chrom1, end1, delim2, chrom2, start2, delim3 = match.groups()

    # Convert to int
    end1, start2 = int(end1), int(start2)

    return delim1, chrom1, end1, delim2, chrom2, start2, delim3


def parse_sv(extended_alt, winsize, chrom_sizes):

    # Parse inputs
    delim1, chrom1, end1, delim2, chrom2, start2, delim3 = parse_alt(extended_alt)
    R = winsize  # window radius

    # Main logic:
    # Translocation contigs
    start1 = end1 + R if delim1 == '[' else end1 - R
    end2 = start2 - R if delim2 == ']' else start2 + R
    # Left reference contig
    end1_1 = end1 + R if end1 > start1 else end1 - R
    # Right reference contrig
    start2_2 = start2 + R if start2 > end2 else start2 - R

    # Sanity checks
    assert 0 <= start1 < chrom_sizes[chrom1], f"Breakpoint specified exceeds {chrom1} dimension"
    assert 0 <=  end1  < chrom_sizes[chrom1], f"Breakpoint specified exceeds {chrom1} dimension"
    assert 0 <= end1_1 < chrom_sizes[chrom1], f"Breakpoint specified exceeds {chrom1} dimension"

    assert 0 <=  start2  < chrom_sizes[chrom2], f"Breakpoint specified exceeds {chrom2} dimension"
    assert 0 <=   end2   < chrom_sizes[chrom2], f"Breakpoint specified exceeds {chrom2} dimension"
    assert 0 <= start2_2 < chrom_sizes[chrom2], f"Breakpoint specified exceeds {chrom2} dimension"

    return (Translocation(Interval(chrom1, start1, end1), Interval(chrom2, start2, end2)),
            Interval(chrom1, start1, end1_1),
            Interval(chrom2, start2_2, end2))


def run_contig_inference(genome_data, contig, res, step, window, model, device):
    dataset = genome_data.get_dataset(contig, window=window, step=step)
    feat = dataset.summarize_features(res)
    pred = run_inference(dataset, model, res, device)
    return feat, pred


def run_inference(dataset, model, res, device):
    pred = {}
    with torch.no_grad():
        for reg, x in dataset:
            mat = model(x.to(device))[0].detach().cpu().numpy()
            bin0 = reg.start // res
            pred[tuple(reg)] = mat2df(bin0, mat)

    lvl_names = {f'level_{i}': x for i, x in enumerate(('chrom', 'start', 'end'))}
    pred = (pd.concat(pred)
            .reset_index(level=[0,1,2])
            .reset_index(drop=True)
            .rename(columns=lvl_names))
    pred = combine_overlaps(pred, res)
    return pred


def mat2df(bin0, mat, k=1):
    sym = (mat + mat.T) * .5  # symmetrize matrix
    rows, cols = np.triu_indices_from(sym, k=k)
    df = pd.DataFrame({
        'row': rows + bin0,
        'col': cols + bin0,
        'val': sym[rows, cols]
    })
    return df


def combine_overlaps(pred, res):
    # combine prediction from different runs
    pred['MID'] = 0.5*(pred['start'] + pred['end']) / res  # mid-point of tile
    pred['DIST'] = abs(0.5*(pred['row'] + pred['col']) - pred['MID'])  # dist of pair from mid-tile
    pred['WEIGHT'] = 1. / (pred['DIST'] + 1.)  # harmonic weights

    def weighted_mean(df):
        return (df['val'] * df['WEIGHT']).sum() / df['WEIGHT'].sum()
    pred = pred.groupby(['chrom', 'row', 'col']).apply(weighted_mean, include_groups=False).reset_index(name='count')

    return pred


def plot_sv(results, fname, tracks=['ATAC', 'CTCF'], h=12, w=9, dpi=300):
    fig, ax = plt.subplots(nrows=6, ncols=1, sharex=True, gridspec_kw={'height_ratios':[4,1] * 3}, layout='compressed')
    fig.set_size_inches(h, w)
    plot_prediction(*results[1], ax[0], ax[1], feat_cols=tracks)
    plot_prediction(*results[0], ax[2], ax[3], feat_cols=tracks)
    plot_prediction(*results[2], ax[4], ax[5], feat_cols=tracks)
    fig.savefig(fname, dpi=dpi)


def plot_prediction(track, hic, ax_hic, ax_track,
                    feat_cols=None, scale=True, point_size=.2, cmap='Reds',
                    brkline={'color': 'black', 'linestyle':'dotted'}):

    contig = track['chrom'].values[0]
    MID = (track['start'] + track['end']).values * 0.5
    bkpt = track['end'].values[track.shape[0] // 2]
    quart = bkpt / 2

    x = (MID[hic['col']] + MID[hic['row']]) * 0.5
    y = (MID[hic['col']] - MID[hic['row']]) * 0.5
    ax_hic.scatter(x, y, c=hic['count'].values, s=point_size, cmap=cmap)
    ax_hic.plot([bkpt, bkpt+quart], [0, quart], **brkline)
    ax_hic.plot([bkpt, bkpt-quart], [0, quart], **brkline)
    ax_hic.set_title(contig)
    #ax_hic.set_aspect('equal')
    ax_hic.axis('off')

    x = MID

    if feat_cols is None:
        feat_cols = [x for x in track.columns if x not in set(['chrom', 'start', 'end'])]

    Nfeats = len(feat_cols)
    for n, f in enumerate(feat_cols):
        y = track[f].values
        if scale: y /= y.max()
        ax_track.plot(x, y + n)

    ax_track.set_ylim(-.1, Nfeats+.1)
    ax_track.set_yticks(np.arange(Nfeats), feat_cols)
    ax_track.spines['left'].set_visible(False)
    ax_track.spines['top'].set_visible(False)
    ax_track.spines['right'].set_visible(False)


def write_cooler(features, interactions, chrom_sizes, fname, assembly='hg38', metadata=None):

    # make bins
    interactions = chrom_index_to_genome(features, interactions)
    coord = get_feature_coord(features)
    bins = features_to_bins(coord, metadata['resolution'], chrom_sizes)

    # make pixels
    coord = pd.merge(coord, bins.reset_index(), how='left')
    interactions['bin1_id'] = coord.loc[interactions['row'], 'index'].values
    interactions['bin2_id'] = coord.loc[interactions['col'], 'index'].values
    pixels = normalize_interactions(interactions)

    cooler.create_cooler(str(fname), bins, pixels, dtypes={'count': np.float64},
                         assembly=assembly, ordered=True, metadata=metadata)


def chrom_index_to_genome(features, interactions):

    # df := chrom, cix = chrom index, index=genome index
    df = features.reset_index()[['chrom', 'index']]
    df['cix'] = df.groupby('chrom').cumcount()

    def chrom2genome(pixels, bins, anchor):
        tmp = pd.merge(pixels, df, how='left', left_on=['chrom', anchor], right_on=['chrom', 'cix'])
        return tmp['index'].values

    interactions['row'] = chrom2genome(interactions, df, 'row')
    interactions['col'] = chrom2genome(interactions, df, 'col')
    return interactions


def get_feature_coord(features):

    def chimeric2ref(coord):
        trans_regex = r'(?P<chrom1>chr[0-9XYM]?):(?P<strand1>[-+]):(?P<start1>[0-9]+)-(?P<end1>[0-9]+)'
        trans_regex += r'(?P<linkage>:[=x]:)'
        trans_regex += r'(?P<chrom2>chr[0-9XYM]?):(?P<strand2>[-+]):(?P<start2>[0-9]+)-(?P<end2>[0-9]+)'
        bedpe = coord['chrom'].str.extract(trans_regex)
        bedpe['start1'] = pd.to_numeric(bedpe['start1'])
        bedpe['end1'] = pd.to_numeric(bedpe['end1'])
        bedpe['width1'] = np.abs(bedpe['end1'] - bedpe['start1'])
        bedpe['strand1'] = bedpe['strand1'].map({'-': -1, '+': 1})
        bedpe['start2'] = pd.to_numeric(bedpe['start2'])
        #bedpe['end2'] = pd.to_numeric(bedpe['end2'])
        bedpe['strand2'] = bedpe['strand2'].map({'-': -1, '+': 1})
        bedpe[['start', 'end']] = coord[['start', 'end']]
        bedpe['overflow'] = bedpe['start'] >= bedpe['width1']
        bedpe['start'] = bedpe.apply(
            lambda row: (row['start'] - row['width1']) * row['strand2'] + row['start2']
            if row['overflow']
            else row['start'] * row['strand1'] + row['start1'],
            axis=1
        )
        bedpe['end'] = bedpe.apply(
            lambda row: (row['end'] - row['width1']) * row['strand2'] + row['start2']
            if row['overflow']
            else row['end'] * row['strand1'] + row['start1'],
            axis=1
        )

        bedpe.rename(columns={'chrom1': 'chrom', 'strand1': 'strand'}, inplace=True)
        bedpe.loc[bedpe['overflow'], ['chrom', 'strand']] = bedpe.loc[bedpe['overflow'], ['chrom2', 'strand2']].values
        return bedpe[['chrom', 'strand', 'start', 'end']]

    def region2genome(coord):
        contig_regex = r'(?P<chrom>chr[0-9XYM]?):(?P<strand>[-+]):(?P<start>[0-9]+)-(?P<end>[0-9]+)'
        bed = coord['chrom'].str.extract(contig_regex)
        bed['start'] = pd.to_numeric(bed['start'])
        bed['end'] = pd.to_numeric(bed['end'])
        bed['strand'] = bed['strand'].map({'-': -1, '+': 1})
        bed['end'] = (coord['end'] * bed['strand'] + bed['start']).values
        bed['start'] = (coord['start'] * bed['strand'] + bed['start']).values
        return bed

    chimeric = features.chrom.str.contains(r':[=x]:')
    bed = pd.concat([chimeric2ref(features.loc[chimeric]),
                     region2genome(features.loc[~chimeric])],
                    axis=0).sort_index()
    ix = bed['end'] < bed['start']
    bed.loc[ix, ['start', 'end']] = bed.loc[ix, ['end', 'start']].values
    return bed


def features_to_bins(bed, binsize, chrom_sizes):
    # get unique bins & keep 1 strand
    bins = bed[['chrom', 'start', 'end']].sort_values(by=['chrom', 'start']).drop_duplicates().reset_index(drop=True)
    ix = bins['end'] < bins['start']
    bins.loc[ix, ['start', 'end']] = bins.loc[ix, ['end', 'start']].values

    # Complete affected chromosomes.
    # Make sure SV are unaffected, 1st and last bin may be smaller
    def extend_rows(group):
        chrom = group['chrom'].values[0]
        chrom_max = chrom_sizes[chrom]

        # Extend towards 0
        start = np.arange(group['start'].min() - binsize, 0, -binsize)
        end = start + binsize
        start, end = np.append(start, [0]), np.append(end, start[[-1]])
        towards_zero = pd.DataFrame({'chrom': chrom, 'start': start, 'end': end})

        # Extend towards chrom_sizes
        start = np.arange(group['end'].max(), chrom_max-binsize, binsize)
        end = start + binsize
        start, end = np.append(start, end[[-1]]), np.append(end, [chrom_max])
        towards_max = pd.DataFrame({'chrom': chrom, 'start': start, 'end': end})

        res = pd.concat([towards_zero, group[['chrom', 'start', 'end']], towards_max])
        return res

    bins = (bins
            .groupby('chrom')
            .apply(extend_rows)
            .reset_index(drop=True))

    # Add missing chromosomes
    skip = set(bins['chrom'].unique())
    bins = [bins]
    for chrom, size in chrom_sizes.items():
        if chrom in skip: continue
        starts = np.arange(0, size, binsize)
        ends = np.append(starts[1:], [size])
        bins.append(pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends}))
    bins = pd.concat(bins, axis=0).sort_values(by=['chrom', 'start']).reset_index(drop=True)

    return bins


def normalize_interactions(interactions):
    "Combine strand and other multimaps"
    ix = interactions['bin1_id'] > interactions['bin2_id']
    interactions.loc[ix, ['bin1_id', 'bin2_id']] = interactions.loc[ix, ['bin2_id', 'bin1_id']].values
    interactions = (interactions
                    .groupby(['bin1_id', 'bin2_id'])['count']
                    .mean()
                    .reset_index()
                    .sort_values(['bin1_id', 'bin2_id']))
    return interactions


def main():

    # Parse arguments
    argv = parse_arguments()
    outdir = Path(argv.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_results = not argv.skip_figures
    write_cool = not argv.skip_cool
    write_csv = not argv.skip_csv

    model_weights = argv.model
    seq_prefix = argv.seq_dir
    atac_path = argv.atac
    atac_norm = argv.atac_norm
    ctcf_path = argv.ctcf
    ctcf_norm = argv.ctcf_norm

    window = int(argv.window)
    n_steps = int(argv.nsteps)
    step = window // n_steps
    resolution = int(argv.resolution)
    chrom_sizes = pd.read_csv(argv.chrom_sizes, sep='\t', header=None, index_col=0)
    chrom_sizes = chrom_sizes.squeeze("columns").to_dict()

    # Parse SVs into concatanation of contigs
    SVs = argv.SV
    contigs = [parse_sv(sv, window, chrom_sizes) for sv in SVs]

    # prepare data
    chromosomes = [(interval.chrom for interval in contig[0].interval) for contig in contigs]
    chromosomes = set(chain(*chromosomes))
    genome_data = GenomeData(seq_prefix=seq_prefix, chromosomes=chromosomes,
                             ctcf=ctcf_path, ctcf_norm=ctcf_norm,
                             atac=atac_path, atac_norm=atac_norm)

    # prepare model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = load_default(model_weights).eval()

    # Run analysis
    results = []
    for contig in contigs:
        res = [run_contig_inference(genome_data, c, resolution, step, window, model, device) for c in contig]
        results.extend(res)
        if plot_results:
            fname = re.sub(':[+-]:', '_', str(contig[0]))
            fname = outdir / (fname.replace(':', '_') + '.png')
            plot_sv(res, fname)
    features, interactions = [pd.concat(df_list) for df_list in zip(*results)]

    # write results
    if write_csv:
        features.to_csv(outdir / 'bins.csv.gz', index=False, compression='gzip')
        interactions.to_csv(outdir / 'interactions.csv.gz', index=False, compression='gzip')

    if write_cool:
        metadata = {'window': window, 'resolution': resolution, 'step': step,
                    'model': str(Path(model_weights).resolve()),
                    'atac': str(Path(atac_path).resolve()), 'atac_norm': str(atac_norm),
                    'ctcf': str(Path(ctcf_path).resolve()), 'ctcf_norm': str(ctcf_norm)}
        fname = outdir / 'results.cool'
        write_cooler(features, interactions, chrom_sizes, fname, argv.assembly, metadata)

    return features, interactions


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Predict simple translocations"
    )

    # Positional arguments
    parser.add_argument(
        "chrom_sizes",
        help="Path to tab-separated chromosome size file"
    )
    parser.add_argument(
        "model",
        help="Path to model weights to instantiate Corigami model"
    )
    parser.add_argument(
        "seq_dir",
        help="Path to directory with DNA sequences for Corigami model"
    )
    parser.add_argument(
        "atac",
        help="Path to ATAC signal bigwig file"
    )
    parser.add_argument(
        "ctcf",
        help="Path to CTCF signal bigwig file"
    )
    parser.add_argument(
        "SV",
        nargs="+",
        help="String defining the translocation, see README for details."
    )

    # Optional arguments
    parser.add_argument(
        "-a", "--assembly",
        default="hg38",
        help="Genome assembly code used to generate data [default: hg38]"
    )
    parser.add_argument(
        "-o", "--outdir",
        default="./",
        help="Directory to write results [default: ./]"
    )
    parser.add_argument(
        "-r", "--resolution",
        type=int,
        default=8192,
        help="Bin size measured in base pairs [default: 8192]"
    )
    parser.add_argument(
        "-w", "--window",
        type=int,
        default=2097152,
        help="Size of window that Corigami predicted [default: 2097152]"
    )
    parser.add_argument(
        "-n", "--nsteps",
        type=int,
        default=4,
        help="Number of steps per window used to shift the prediction region [default: 4]"
    )
    parser.add_argument(
        "--atac_norm",
        default="log",
        help="Function used to normalize ATAC signal [default: log]"
    )
    parser.add_argument(
        "--ctcf_norm",
        help="Function used to normalize CTCF signal"
    )
    parser.add_argument(
        "--skip_cool",
        action="store_true",
        help="Do not write results in cool format"
    )
    parser.add_argument(
        "--skip_figures",
        action="store_true",
        help="Do not plot figures of translocation"
    )
    parser.add_argument(
        "--skip_csv",
        action="store_true",
        help="Do not write results in csv format"
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()

