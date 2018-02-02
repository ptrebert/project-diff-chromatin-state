#!/usr/bin/env python3

import os as os
import sys as sys
import gzip as gz
import argparse as argp
import operator as op
import traceback as trb

import numpy as np
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--de-genes', '-deg', type=str, dest='degenes')
    parser.add_argument('--annotation', '-ann', type=str, dest='annotation')
    parser.add_argument('--diff-suffix', '-do', type=str, dest='diffout')
    parser.add_argument('--norm-suffix', '-no', type=str, dest='normout')
    parser.add_argument('--out-root', '-or', type=str, dest='outroot')
    parser.add_argument('--fc-threshold', '-fc', type=float, default=2., dest='fcthreshold',
                        help='Specify (log2) fold-change threshold to select'
                             ' differentially expressed genes. Default: 2')
    parser.add_argument('--pv-threshold', '-pv', type=float, default=0.01, dest='pvthreshold',
                        help='Specify p-value threshold (corrected for multiple testing)'
                             ' to select differentially expressed genes. Default: 0.01')
    parser.add_argument('--select-chroms', '-slc', type=str, default='chr[0-9X]+', dest='selectchroms',
                        help='Specify regular expression to select chromosomes to dump to file.'
                             ' Default: chr[0-9X]+ (autosomes and chrX)')
    args = parser.parse_args()
    return args


def read_deseq_table(fpath):
    """
    :param fpath:
    :return:
    """
    header = ['geneid', 'mean', 'lfc', 'lfcse', 'stat', 'pv', 'padj']
    with open(fpath, 'r') as table:
        _ = table.readline()
        data = pd.read_csv(table, delimiter='\t', header=None, names=header)
    num_records = data.shape[0]
    # Deseq2 marks identified outliers with NA pvalues
    # it is not recommended to just throw these away
    # when dealing with small sample numbers - keep separate
    outlier = data.loc[pd.isnull(data).any(axis=1), :].copy()
    data.dropna(axis=0, how='any', inplace=True)
    num_outliers = outlier.shape[0]
    #print('Outlier in dataset {} ({}%)'.format(num_outliers, round((num_outliers / num_records) * 100, 1)))
    return data, outlier


def read_annotation(fpath):
    """
    :param fpath:
    :return:
    """
    get_loc = op.itemgetter(*(0, 3, 4, 6))
    gene_bodies = []
    gene_promoters = []
    with gz.open(fpath, 'rt') as gtf:
        for line in gtf:
            cols = line.split()
            try:
                if cols[2] == 'gene':
                    chrom, start, end, strand = get_loc(cols)
                    start, end = int(start), int(end)
                    if strand == '+':
                        prom_start, prom_end = start - 2500, start + 500
                    elif strand == '-':
                        prom_start, prom_end = end - 500, end + 2500
                    else:
                        raise ValueError('Invalid value for strand in line {}'.format(line.strip()))
                    geneid = 'noid'
                    symbol = 'nosymbol'
                    for idx, attr in enumerate(cols[8:], start=8):
                        if attr == 'gene_id':
                            geneid = cols[idx+1].strip('";')
                        elif attr == 'gene_name':
                            symbol = cols[idx+1].strip('";')
                        else:
                            continue
                    gene_bodies.append((chrom, start, end, geneid, 0, strand, symbol))
                    gene_promoters.append((chrom, prom_start, prom_end, geneid, 0, strand, symbol))
            except IndexError:
                continue
    header = ['chrom', 'start', 'end', 'geneid', 'score', 'strand', 'symbol']
    body_ann = pd.DataFrame(gene_bodies, columns=header)
    prom_ann = pd.DataFrame(gene_promoters, columns=header)
    return body_ann, prom_ann


def dump_genes_to_file(fpath, genes, locations, chrom_re):
    """
    :param fpath:
    :param genes:
    :param locations:
    :return:
    """
    dump = genes.merge(locations, how='outer', on='geneid')
    # drop genes that are annotated but not in the dataset "genes"
    dump.dropna(how='any', axis=0, inplace=True)
    dump.sort_values(by=['chrom', 'start', 'end'], inplace=True)
    dump = dump[['chrom', 'start', 'end', 'geneid', 'lfc', 'strand', 'symbol', 'padj']]
    dump.columns = ['chrom', 'start', 'end', 'name', 'log2fc', 'strand', 'symbol', 'pv_adj']
    select_chrom = dump['chrom'].str.match(chrom_re, as_indexer=True)
    dump = dump.loc[select_chrom, :]
    with open(fpath, 'w') as bedfile:
        _ = bedfile.write('#')
        dump.to_csv(bedfile, sep='\t', header=True, index=False)
    return


def main():
    """
    :return:
    """
    cargs = parse_command_line()
    de_table, outlier = read_deseq_table(cargs.degenes)
    outlier.fillna(-1, inplace=True)
    above_fc = np.array(de_table['lfc'].abs() > cargs.fcthreshold, dtype=np.bool)
    below_pv = np.array(de_table['padj'] < cargs.pvthreshold, dtype=np.bool)
    select_rows = np.logical_and(above_fc, below_pv)
    de_genes = de_table.loc[select_rows, :].copy()
    nm_genes = de_table.loc[~select_rows, :].copy()
    nm_genes = pd.concat([nm_genes, outlier], axis=0, ignore_index=False)
    bodies, promoters = read_annotation(cargs.annotation)

    base_out = os.path.splitext(os.path.basename(cargs.degenes))
    base_out = os.path.join(cargs.outroot, base_out[0])

    dump_de_bodies = base_out + cargs.diffout + '_body.bed'
    dump_genes_to_file(dump_de_bodies, de_genes, bodies, cargs.selectchroms)
    dump_de_promoters = base_out + cargs.diffout + '_promoter.bed'
    dump_genes_to_file(dump_de_promoters, de_genes, promoters, cargs.selectchroms)

    dump_nm_bodies = base_out + cargs.normout + '_body.bed'
    dump_genes_to_file(dump_nm_bodies, nm_genes, bodies, cargs.selectchroms)
    dump_nm_promoters = base_out + cargs.normout + '_promoter.bed'
    dump_genes_to_file(dump_nm_promoters, nm_genes, promoters, cargs.selectchroms)

    return


if __name__ == '__main__':
    try:
        main()
    except Exception:
        trb.print_exc()
        sys.exit(1)
    else:
        sys.exit(0)
