#!/usr/bin/env python3

import os as os
import sys as sys
import gzip as gz
import argparse as argp
import operator as op
import traceback as trb
import csv as csv

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


def read_deseq_table(fpath, t_padj, t_fc):
    """
    :param fpath:
    :return:
    """
    assert 0 < t_padj < 1, 'Invalid argument for p-value threshold: {}'.format(t_padj)
    assert t_fc > 0, 'Invalid argument for fold-change threshold: {}'.format(t_fc)
    header = ['name', 'baseMean', 'log2fc', 'fcSE', 'statistic', 'pvalue', 'padj']
    collector = []
    with open(fpath, 'r', newline='') as table:
        _ = table.readline()
        rows = csv.DictReader(table, fieldnames=header, delimiter='\t',
                              dialect=csv.unix_dialect)
        for row in rows:
            row['name'] = row['name'].split('.')[0]
            if row['log2fc'] == 'NA':
                # this can happen for genes that
                # do not get any read counts across
                # all samples
                row['is_de'] = -1
                row['padj'] = 1
                row['pvalue'] = 1
                row['statistic'] = 0
                row['fcSE'] = 0
                row['log2fc'] = 0
                assert int(row['baseMean']) == 0, 'Weird gene in DESeq output: {}'.format(row)
            elif row['padj'] == 'NA':
                row['is_de'] = 0
                row['padj'] = -1
            else:
                fc = abs(float(row['log2fc']))
                padj = float(row['padj'])
                if fc > t_fc and padj < t_padj:
                    row['is_de'] = 1
                else:
                    row['is_de'] = -1
            collector.append(row)
    df = pd.DataFrame.from_records(collector, columns=header + ['is_de'])
    return df


def read_gtf_annotation(fpath):
    """
    :param fpath:
    :return:
    """
    get_loc = op.itemgetter(*(0, 3, 4, 6))
    gene_bodies = []
    gene_promoters = []
    gene_loci = []
    with gz.open(fpath, 'rt') as gtf:
        for line in gtf:
            cols = line.split()
            try:
                if cols[2] == 'gene':
                    chrom, start, end, strand = get_loc(cols)
                    start, end = int(start), int(end)
                    if strand == '+':
                        prom_start, prom_end = start - 2500, start + 500
                        loc_start, loc_end = start - 2500, end
                    elif strand == '-':
                        prom_start, prom_end = end - 500, end + 2500
                        loc_start, loc_end = start, end + 2500
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
                    gene_loci.append((chrom, loc_start, loc_end, geneid, 0, strand, symbol))
            except IndexError:
                continue
    header = ['chrom', 'start', 'end', 'geneid', 'score', 'strand', 'symbol']
    body_ann = pd.DataFrame(gene_bodies, columns=header)
    prom_ann = pd.DataFrame(gene_promoters, columns=header)
    locus_ann = pd.DataFrame(gene_loci, columns=header)
    return body_ann, prom_ann, locus_ann


def read_bed_annotation(fpath):
    """
    :param fpath:
    :return:
    """
    get_loc = op.itemgetter(*('chrom', 'start', 'end', 'strand'))
    get_desc = op.itemgetter(*('name', 'symbol', 'score'))
    gene_bodies = []
    gene_promoters = []
    gene_loci = []
    with open(fpath, 'r', newline='') as bed:
        header = bed.readline().strip('#\n ').split('\t')
        rows = csv.DictReader(bed, fieldnames=header, delimiter='\t',
                              dialect=csv.unix_dialect)
        for row in rows:
            c, s, e, strand = get_loc(row)
            s, e = int(s), int(e)
            name, symbol, score = get_desc(row)
            name = name.split('.')[0]  # better safe than sorry...
            if strand == '+':
                prom_start, prom_end = s - 2500, s + 500
                loc_start, loc_end = s - 2500, e
            elif strand == '-':
                prom_start, prom_end = e - 500, e + 2500
                loc_start, loc_end = s, e + 2500
            else:
                raise ValueError('Invalid value for strand in line {}'.format(strand))
            gene_bodies.append((c, s, e, name, score, strand, symbol))
            gene_promoters.append((c, prom_start, prom_end, name, score, strand, symbol))
            gene_loci.append((c, loc_start, loc_end, name, score, strand, symbol))

    header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'symbol']
    body_ann = pd.DataFrame(gene_bodies, columns=header)
    prom_ann = pd.DataFrame(gene_promoters, columns=header)
    locus_ann = pd.DataFrame(gene_loci, columns=header)
    assert body_ann.shape[0] == prom_ann.shape[0] == locus_ann.shape[0], \
        'Shape mismatch in annotation DataFrames'
    return body_ann, prom_ann, locus_ann


def dump_genes_to_file(fpath, genes, locations, chrom_re):
    """
    :param fpath:
    :param genes:
    :param locations:
    :return:
    """
    shared_cols = sorted(set(genes.columns).intersection(set(locations.columns)))
    dump = genes.merge(locations, how='outer', on=shared_cols)
    # drop genes that are annotated but not in the dataset "genes"
    dump.dropna(how='any', axis=0, inplace=True)
    dump.sort_values(by=['chrom', 'start', 'end'], inplace=True)
    dump_columns = ['chrom', 'start', 'end', 'name', 'log2fc', 'strand', 'symbol', 'padj', 'baseMean']
    if '_all_' in fpath:
        dump_columns.extend(['is_de'])
    dump = dump[dump_columns]
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

    if cargs.annotation.endswith('gtf.gz'):
        bodies, promoters, loci = read_gtf_annotation(cargs.annotation)
    elif cargs.annotation.endswith('.bed'):
        bodies, promoters, loci = read_bed_annotation(cargs.annotation)
    else:
        raise ValueError('Unexpected file format for annotation: {}'.format(cargs.annotation))

    genes = read_deseq_table(cargs.degenes, cargs.pvthreshold, cargs.fcthreshold)
    genes = genes.loc[genes['name'].isin(loci['name']), :].copy()

    assert genes.shape[0] == loci.shape[0], \
        'Too many genes in annotation: {} v {}'.format(genes.shape[0], loci.shape[0])

    de_genes = genes.loc[genes['is_de'] == 1, :].copy()
    nm_genes = genes.loc[genes['is_de'] < 1, :].copy()

    base_out = os.path.splitext(os.path.basename(cargs.degenes))
    base_out = os.path.join(cargs.outroot, base_out[0])

    dump_de_bodies = base_out + cargs.diffout + '_body.bed'
    dump_genes_to_file(dump_de_bodies, de_genes, bodies, cargs.selectchroms)
    dump_de_promoters = base_out + cargs.diffout + '_promoter.bed'
    dump_genes_to_file(dump_de_promoters, de_genes, promoters, cargs.selectchroms)
    dump_de_loci = base_out + cargs.diffout + '_locus.bed'
    dump_genes_to_file(dump_de_loci, de_genes, loci, cargs.selectchroms)

    dump_nm_bodies = base_out + cargs.normout + '_body.bed'
    dump_genes_to_file(dump_nm_bodies, nm_genes, bodies, cargs.selectchroms)
    dump_nm_promoters = base_out + cargs.normout + '_promoter.bed'
    dump_genes_to_file(dump_nm_promoters, nm_genes, promoters, cargs.selectchroms)
    dump_nm_loci = base_out + cargs.normout + '_locus.bed'
    dump_genes_to_file(dump_nm_loci, nm_genes, loci, cargs.selectchroms)

    dump_all_bodies = base_out + '_all' + '_body.bed'
    dump_genes_to_file(dump_all_bodies, genes, bodies, cargs.selectchroms)
    dump_all_promoters = base_out + '_all' + '_promoter.bed'
    dump_genes_to_file(dump_all_promoters, genes, promoters, cargs.selectchroms)
    dump_all_loci = base_out + '_all' + '_locus.bed'
    dump_genes_to_file(dump_all_loci, genes, loci, cargs.selectchroms)

    return


if __name__ == '__main__':
    try:
        main()
    except Exception:
        trb.print_exc()
        sys.exit(1)
    else:
        sys.exit(0)
