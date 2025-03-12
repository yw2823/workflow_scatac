import os
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

import torch
from torch.utils.data import Dataset, DataLoader

import numpy as np
import pandas as pd
import cooler

from corigami.inference.utils.model_utils import load_default
from corigami.data.data_feature import GenomicFeature, SequenceFeature

from snakemake.script import snakemake


class InputData(Dataset):
    def __init__(self, regions, seq_prefix, ctcf, atac, ctcf_norm=None, atac_norm='log'):
        # cuda = 'cuda:{}'.format(os.environ['CUDA_VISIBLE_DEVICES'])
        # self.device = torch.device(cuda if torch.cuda.is_available() else 'cpu')
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.regions = regions
        chromosomes = list(set(reg.chrom for reg in self.regions))
        chromosomes.sort()
        self.seq = {chrom: SequenceFeature(path=Path(seq_prefix) / f'{chrom}.fa.gz') for chrom in chromosomes}
        self.ctcf = GenomicFeature(path=ctcf, norm=ctcf_norm)
        self.atac = GenomicFeature(path=atac, norm=atac_norm)

    def __len__(self):
        return len(self.regions)

    def _process_input(self, region):
        seq = self.seq[region.chrom].get(region.start, region.end)
        seq = torch.tensor(seq).unsqueeze(0)
        ctcf = self.ctcf.get(region.chrom, region.start, region.end)
        ctcf = torch.tensor(np.nan_to_num(ctcf, 0))
        atac = self.atac.get(region.chrom, region.start, region.end)
        atac = torch.tensor(atac)
        features = [ctcf, atac]
        features = torch.cat([feat.unsqueeze(0).unsqueeze(2) for feat in features], dim = 2)
        inputs = torch.cat([seq, features], dim = 2)
        return inputs.to(self.device)

    def __getitem__(self, idx):
        region = self.regions[idx]
        return region, self._process_input(region)

def read_bed(path, chroms):
    print("Loading windows...")
    windows = pd.read_table(path, sep='\t', index_col=False, header=None, names=['chrom', 'start', 'end'])
    windows = windows[windows['chrom'].isin(chroms)]
    return windows.reset_index(drop=True)


def mat2df(bin0, mat, k=0):
    sym = (mat + mat.T) * .5  # symmetrize matrix
    rows, cols = np.triu_indices_from(sym, k=k)
    df = pd.DataFrame({
        'row': rows + bin0,
        'col': cols + bin0,
        'val': sym[rows, cols]
    })
    return df


def run_inference(windows, snakemake):

    print("Running Inference...")
    model = load_default(snakemake.input['model']).eval()
    # cuda = 'cuda:{}'.format(os.environ['CUDA_VISIBLE_DEVICES'])
    # device = torch.device(cuda if torch.cuda.is_available() else 'cpu')
    # model = model.to(device).eval()
    regions = list(windows.itertuples(index=False, name='Region'))
    dataset = InputData(regions,
                        seq_prefix=snakemake.input['seq'],
                        ctcf=snakemake.input['ctcf'],
                        atac=snakemake.input['atac'])
    # dataloader = DataLoader(dataset, batch_size=1, shuffle=False)

    pred = {}
    res = snakemake.params['resolution']
    with torch.no_grad():
        for reg, x in dataset:
            mat = model(x)[0].detach().cpu().numpy()
            bin0 = reg.start // res
            pred[tuple(reg)] = mat2df(bin0, mat)

    lvl_names = {f'level_{i}': x for i, x in enumerate(('chrom', 'start', 'end'))}
    pred = (pd.concat(pred)
            .reset_index(level=[0,1,2])
            .reset_index(drop=True)
            .rename(columns=lvl_names))
    return pred

def combine_overlaps(pred, res):
    # combine prediction from different runs
    pred['MID'] = 0.5*(pred['start'] + pred['end']) / res  # mid-point of tile
    pred['DIST'] = abs(0.5*(pred['row'] + pred['col']) - pred['MID'])  # dist of pair from mid-tile
    pred['WEIGHT'] = 1. / (pred['DIST'] + 1.)  # harmonic weights

    def weighted_mean(group):
        return (group['val'] * group['WEIGHT']).sum() / group['WEIGHT'].sum()
    pred = pred.groupby(['chrom', 'row', 'col']).apply(weighted_mean, include_groups=True).reset_index(name='count')

    return pred

def tile_chrom(chrom, N, res):
    start = np.arange(0, N, res, dtype=int)
    return pd.DataFrame(data={'chrom': chrom,
                              'start': start,
                              'end': np.minimum(start + res, N),
                              'ix': np.arange(len(start), dtype=int)})

def tile_genome(chrom_sizes, res, chromosomes=None):
    sizes = pd.read_csv(chrom_sizes, sep='\t', header=None, names=['chrom', 'length'])
    if chromosomes is not None:
        sizes = sizes.loc[sizes.chrom.isin(chromosomes)]
    tiles = pd.concat([tile_chrom(chrom, N, res) for chrom, N in sizes.itertuples(index=False)], ignore_index=True)
    return tiles.reset_index()


def index_pair(pixels, bins, anchor, bin_id, chrom_col='chrom', chrom_id='ix', genome_id='index'):
    df = pd.merge(pixels, bins[[chrom_col, chrom_id, genome_id]], how='left',
                  left_on=[chrom_col, anchor], right_on=[chrom_col, chrom_id])
    return df.rename(columns={'index': bin_id}).drop('ix', axis=1)

def write_cooler(snakemake, tiles, hic):
    print("Writing cooler...")
    hic = index_pair(hic, tiles, 'row', 'bin1_id')
    hic = index_pair(hic, tiles, 'col', 'bin2_id')

    # keep only relevant columns
    tiles = tiles[['chrom', 'start', 'end']]
    hic = hic[['bin1_id', 'bin2_id', 'count']]
    metadata = {'resolution': snakemake.params['resolution'],
                'window': snakemake.params['window'],
                'step': snakemake.params['step'],
                'model': str(Path(snakemake.input['model']).resolve()),
                'atac': str(Path(snakemake.input['atac']).resolve()),
                'ctcf': str(Path(snakemake.input['ctcf']).resolve())}
    cooler.create_cooler(snakemake.output[0], tiles, hic, dtypes={'count': np.float64},
                         assembly=snakemake.params['genome'], ordered=True,
                         metadata=metadata)


def main(snakemake):
    windows = read_bed(snakemake.input['windows'], snakemake.params['chromosomes'])
    hic = run_inference(windows, snakemake)
    hic = combine_overlaps(hic, snakemake.params['resolution'])
    tiles = tile_genome(snakemake.params['chrom_sizes'], snakemake.params['resolution'], snakemake.params['chromosomes'])
    write_cooler(snakemake, tiles, hic)


with open(snakemake.log[0], "w") as logfile:
    with redirect_stdout(logfile), redirect_stderr(logfile):
        main(snakemake)
