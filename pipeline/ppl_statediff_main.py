# coding=utf-8

import os as os
import re as re
import csv as csv
import operator as op
import collections as col
import itertools as itt
import time as ti

from ruffus import *


def link_deep_samples(basedir, sample_ids, libraries, outdir, input_complete):
    """
    :param basedir:
    :param sample_ids:
    :return:
    """
    replace_ids = {'01_HepG2_LiHG_Ct1': '01_Hc01_LiHG_Ct',
                   '01_HepG2_LiHG_Ct2': '01_Hc02_LiHG_Ct'}

    if input_complete:
        all_links = os.listdir(outdir)
        all_links = [os.path.join(outdir, fn) for fn in all_links]
        return all_links

    sample_complete = col.defaultdict(list)
    linked_files = []
    os.makedirs(outdir, exist_ok=True)
    for root, dirs, data_files in os.walk(basedir):
        if root.endswith('qualitycontrol'):
            continue
        if data_files:
            for df in data_files:
                if df.endswith('.bam') or df.endswith('.bai'):
                    this_sample, library, _, rep = df.rsplit('_', 3)
                    if library not in libraries:
                        continue
                    if 'bwamem' not in df and this_sample != '01_HepG2_LiHG_Ct1':
                        print('Skipping {}'.format(df))
                        continue
                    if '_S_2' in df:  # ignore HepG2 Ct1 duplicates
                        continue
                    if 'GALv2' in df:
                        if library in ['H3K27me3', 'H3K4me1', 'H3K4me3', 'Input']:
                            continue
                    if this_sample in sample_ids:
                        this_analyte = df.split('.')[0]
                        this_ext = df.split('.')[-1]
                        if this_sample in replace_ids:
                            this_analyte = this_analyte.replace(this_sample, replace_ids[this_sample])
                        out_path = os.path.join(outdir, this_analyte + '.' + this_ext)
                        if not os.path.islink(out_path):
                            os.symlink(os.path.join(root, df), out_path)
                        if df.endswith('.bam'):
                            sample_complete[this_sample].append(library)
                        linked_files.append(out_path)
    for k, v in sample_complete.items():
        assert len(v) == 7, 'Incomplete sample: {} - {}'.format(k, sorted(v))
    return linked_files


def prep_chromhmm_table(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    entries = dict()
    for bamfile in inputfiles:
        fn = os.path.basename(bamfile)
        cell = fn.rsplit('_', 3)[0]
        mark = fn.split('_')[4]
        entries[(cell, mark)] = fn
    with open(outputfile, 'w') as table:
        writer = csv.writer(table, delimiter='\t')
        for (c, m) in sorted(entries.keys()):
            if m == 'Input':
                continue
            signal = entries[(c, m)]
            control = entries[(c, 'Input')]
            row = [c, m, signal, control]
            writer.writerow(row)
    return outputfile


def check_run_binbam(sample_table, binary_dir):
    """
    :param sample_table:
    :param binary_dir:
    :return:
    """
    samples = set()
    if not os.path.isfile(sample_table):
        return True
    with open(sample_table) as tab:
        for line in tab:
            if line.strip():
                smp, _ = line.split('\t', 1)
                samples.add(smp)
    bin_files = os.listdir(binary_dir)
    chrom_count = col.Counter()
    for bf in bin_files:
        smp, _, _ = bf.rsplit('_', 2)
        chrom_count[smp] += 1
    chrom_missing = len(set(chrom_count.values())) > 1
    sample_missing = len(samples) != len(chrom_count.keys())
    incomplete = chrom_missing or sample_missing
    return incomplete


def make_ecs_count_params(bamfolder, samples, outdir, cmd, jobcall):
    """
    :param bamfolder:
    :param samples:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    sort_samples = col.defaultdict(list)
    samples = [s.strip() for s in samples.split()]
    for bf in os.listdir(bamfolder):
        if bf.endswith('.bam'):
            if any(bf.startswith(x) for x in samples):
                smp, mark, _, _ = bf.rsplit('_', 3)
                sort_samples[smp].append((os.path.join(bamfolder, bf), mark))
    args = []
    for smp, datafiles in sort_samples.items():
        lab_input = ''
        inputfiles = []
        datafiles = sorted(datafiles)
        for fp, mark in datafiles:
            lab_input += '-m {}:{} '.format(mark, fp)
            inputfiles.append(fp)
        outfile = smp + '.ecs-cnt.txt'
        outpath = os.path.join(outdir, outfile)
        tmp = cmd.format(**{'labeled_inputs': lab_input})
        args.append([inputfiles, outpath, tmp, jobcall])
    return args


def make_ecs_norm_params(matfolder, cmd, jobcall):
    """
    :param matfolder:
    :param cmd:
    :param jobcall:
    :return:
    """
    file_args = []
    cli_args = ''
    for fn in os.listdir(matfolder):
        if fn.endswith('.ecs-cnt.txt'):
            fp = os.path.join(matfolder, fn)
            cli_args += '-c {} '.format(fp)
            file_args.append(fp)
    tmp = cmd.format(**{'matrices': cli_args})
    args = [[file_args, os.path.join(matfolder, 'ecs_norm_mat.chk'), tmp, jobcall]]
    return args


def make_ecs_seg_params(matfolder, out_base, nstates, cmd, jobcall):
    """
    :param matfolder:
    :param out_base:
    :param nstates:
    :param cmd:
    :param jobcall:
    :return:
    """
    file_args = []
    cli_args = ''
    for fn in sorted(os.listdir(matfolder)):
        if fn.endswith('_norm.txt'):
            fp = os.path.join(matfolder, fn)
            file_args.append(fp)
            smp, _ = fn.split('.', 1)
            cli_args += '-c {}:{} '.format(smp, fp)
    args = []
    for n in nstates:
        outdir = os.path.join(out_base, 'seg_{}'.format(n))
        os.makedirs(outdir, exist_ok=True)
        tmp = cmd.format(**{'labeled_input': cli_args,
                            'num_states': n, 'outdir': outdir})
        outfile = os.path.join(outdir, 'ecs_deep_seg_{}.chk'.format(n))
        args.append([file_args, outfile, tmp, jobcall])
    return args


def make_salmon_calls(input_folder, base_out, cmd, jobcall):
    """
    :param input_folder:
    :param base_out:
    :param cmd:
    :param jobcall:
    :return:
    """
    expfiles = os.listdir(input_folder)
    if not expfiles:
        return []
    # put files in appropriate buckets
    samples = dict()
    for fname in expfiles:
        fpath = os.path.join(input_folder, fname)
        parts = fname.split('_', 4)
        sid = '_'.join(parts[0:4])
        if sid not in samples:
            samples[sid] = col.defaultdict(list)
        if '_R1_' in parts[4]:
            samples[sid]['R1'].append(fpath)
        else:
            samples[sid]['R2'].append(fpath)
    args = []
    for sid, reads in samples.items():
        outpath = os.path.join(base_out, sid)
        reads1 = ' '.join(sorted(reads['R1']))
        reads2 = ' '.join(sorted(reads['R2']))
        inputfiles = reads['R1'] + reads['R2']
        outputfile = os.path.join(outpath, 'quant.sf')
        tmp = cmd.format(**{'outpath': outpath, 'reads1': reads1, 'reads2': reads2})
        args.append([inputfiles, outputfile, tmp, jobcall])
    return args


def make_zerone_deep_params(inputfolder, outputfolder, cmd, jobcall):
    """
    :param inputfolder:
    :param outputfolder:
    :param cmd:
    :param jobcall:
    :return:
    """
    all_files = os.listdir(inputfolder)
    if not all_files:
        return []
    replicates = col.defaultdict(list)
    args = []
    for fn in all_files:
        if fn.endswith('.bam'):
            smp, lib, _, _ = fn.rsplit('_', 3)
            smp_key = smp[:5] + 'NN' + smp[7:]
            replicates[(smp_key, lib)].append(fn)
    for k, v in replicates.items():
        if k[1] == 'Input':
            continue
        input_bams = replicates[(k[0], 'Input')]
        control_files = ','.join([os.path.join(inputfolder, fn) for fn in input_bams])
        if len(v) > 1:
            out_base = k[0] + '_' + k[1] + '.zerone-pk.txt'
            signal_files = ','.join([os.path.join(inputfolder, fn) for fn in v])
        else:
            out_base = v[0].rsplit('_', 2)[0] + '.zerone-pk.txt'
            signal_files = os.path.join(inputfolder, v[0])
        tmp = cmd.format(**{'signal_files': signal_files, 'control_files': control_files})
        inputfiles = [os.path.join(inputfolder, fn) for fn in v]
        inputfiles.extend([os.path.join(inputfolder, fn) for fn in input_bams])
        outpath = os.path.join(outputfolder, out_base)
        args.append([inputfiles, outpath, tmp, jobcall])
    return args


def make_pepr_diff_params(inputfolder, outputfolder, subset, cmd, jobcall):
    """
    :param inputfolder:
    :param outputfolder:
    :param subset:
    :param cmd:
    :param jobcall:
    :return:
    """
    peaktypes = {'H3K27me3': 'broad', 'H3K9me3': 'broad', 'H3K36me3': 'broad'}
    comparisons = [('HG', 'He'), ('HG', 'Ma'), ('HG', 'Mo'),
                   ('He', 'Ma'), ('He', 'Mo'), ('Mo', 'Ma')]
    marks = dict((t[0], col.defaultdict(list)) for t in comparisons)
    marks['Ma'] = col.defaultdict(list)
    bam_input = col.defaultdict(list)
    for fn in os.listdir(inputfolder):
        if fn.endswith('.bam'):
            sid, lib, _, _ = fn.rsplit('_', 3)
            if sid in subset:
                celltype = sid[10:12]
                if lib == 'Input':
                    bam_input[celltype].append(fn)
                else:
                    marks[celltype][lib].append(fn)
    args = []
    for g1, g2 in comparisons:
        g1_marks = marks[g1]
        g2_marks = marks[g2]
        inputs1 = ','.join(sorted(bam_input[g1]))
        inputs2 = ','.join(sorted(bam_input[g2]))
        out_dir = os.path.join(outputfolder, '{}_vs_{}'.format(g1, g2))
        os.makedirs(out_dir, exist_ok=True)
        for lib, files1 in g1_marks.items():
            files2 = g2_marks[lib]
            infos = {'chip1': ','.join(sorted(files1)), 'chip2': ','.join(sorted(files2)),
                     'input1': inputs1, 'input2': inputs2, 'name': '{}_vs_{}_{}'.format(g1, g2, lib),
                     'inputdir': inputfolder, 'outputdir': out_dir, 'peaktype': peaktypes.get(lib, 'sharp')}
            tmp = cmd.format(**infos)
            inputfiles = [os.path.join(inputfolder, fn) for fn in files1]
            inputfiles.extend([os.path.join(inputfolder, fn) for fn in files2])
            outputfile = '{}_vs_{}_{}__PePr_parameters.txt'.format(g1, g2, lib)
            outpath = os.path.join(out_dir, outputfile)
            args.append([inputfiles, outpath, tmp, jobcall])
    return args


def make_thor_diff_params(inputfolder, outputfolder, subset, cmd, jobcall):
    """
    :param inputfolder:
    :param outputfolder:
    :param subset:
    :param cmd:
    :param jobcall:
    :return:
    """

    raw_cfg = """#rep1
            {sample1_mark}
            #rep2
            {sample2_mark}
            #genome
            /TL/deep/fhgfs/projects/pebert/thesis/refdata/rgtdata/hg38/genome_hg38.fa
            #chrom_sizes
            /TL/deep/fhgfs/projects/pebert/thesis/refdata/rgtdata/hg38/chrom.sizes.hg38
            #inputs1
            {sample1_input}
            #inputs2
            {sample2_input}
            """
    raw_cfg = raw_cfg.replace(' ', '')
    config_folder = os.path.join(os.path.split(outputfolder)[0], 'run_configs')
    os.makedirs(config_folder, exist_ok=True)
    comparisons = [('HG', 'He'), ('HG', 'Ma'), ('HG', 'Mo'),
                   ('He', 'Ma'), ('He', 'Mo'), ('Ma', 'Mo')]
    marks = dict((t[0], col.defaultdict(list)) for t in comparisons)
    marks['Mo'] = col.defaultdict(list)
    bam_input = col.defaultdict(list)
    for fn in os.listdir(inputfolder):
        if fn.endswith('.bam'):
            sid, lib, _, _ = fn.rsplit('_', 3)
            if sid in subset:
                celltype = sid[10:12]
                if lib == 'Input':
                    bam_input[celltype].append(os.path.join(inputfolder, fn))
                else:
                    marks[celltype][lib].append(os.path.join(inputfolder, fn))
    args = []
    for g1, g2 in comparisons:
        g1_marks = marks[g1]
        g2_marks = marks[g2]
        inputs1 = '\n'.join(sorted(bam_input[g1]))
        inputs2 = '\n'.join(sorted(bam_input[g2]))
        out_dir = os.path.join(outputfolder, '{}_vs_{}'.format(g1, g2))
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        for lib, files1 in g1_marks.items():
            files2 = g2_marks[lib]
            run_name = '{}_vs_{}_{}'.format(g1, g2, lib)
            cfg = {'sample1_mark': '\n'.join(sorted(files1)), 'sample2_mark': '\n'.join(sorted(files2)),
                   'sample1_input': inputs1, 'sample2_input': inputs2}
            tmp_cfg = raw_cfg.format(**cfg)
            cfg_path = os.path.join(config_folder, run_name + '.thor.cfg')

            with open(cfg_path, 'w') as dump:
                _ = dump.write(tmp_cfg)

            cmd_opt = {'outputdir': out_dir, 'name': run_name}
            tmp = cmd.format(**cmd_opt)
            outputfile = run_name + '-setup.info'
            outpath = os.path.join(out_dir, outputfile)
            args.append([cfg_path, outpath, tmp, jobcall])
    return args


def make_pepr_post_params(inputroot, cmd, jobcall):
    """
    :param inputroot:
    :param cmd:
    :param jobcall:
    :return:
    """
    args = []
    for root, dirs, datafiles in os.walk(inputroot):
        if datafiles:
            datafiles = [df for df in datafiles if df.endswith('__PePr_parameters.txt')]
            for df in datafiles:
                param_file = os.path.join(root, df)
                records = col.defaultdict(list)
                with open(param_file, 'r') as tsv:
                    _ = tsv.readline()
                    for line in tsv:
                        k, v = line.strip().split()[:2]
                        records[k].append(v)
                output_file = param_file.replace('parameters.txt', 'post-peak1.chk')
                data_dir = records['input-directory'][0]
                peak1 = param_file.replace('parameters.txt', 'chip1_peaks.bed')
                chip1 = ','.join(sorted([os.path.join(data_dir, fn) for fn in records['chip1']]))
                input1 = ','.join(sorted([os.path.join(data_dir, fn) for fn in records['input1']]))
                tmp1 = cmd.format(**{'peak': peak1, 'chip': chip1, 'input': input1})
                args.append([param_file, output_file, tmp1, jobcall])

                output_file = param_file.replace('parameters.txt', 'post-peak2.chk')
                peak2 = param_file.replace('parameters.txt', 'chip2_peaks.bed')
                chip2 = ','.join(sorted([os.path.join(data_dir, fn) for fn in records['chip2']]))
                input2 = ','.join(sorted([os.path.join(data_dir, fn) for fn in records['input2']]))
                tmp2 = cmd.format(**{'peak': peak2, 'chip': chip2, 'input': input2})
                args.append([param_file, output_file, tmp2, jobcall])
    return args


def filter_sort_pepr_peaks(inputfile, outputfile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    if inputfile.endswith('peak1.chk'):
        peakfile = inputfile.replace('post-peak1.chk', 'chip1_peaks.bed.passed')
    else:
        peakfile = inputfile.replace('post-peak2.chk', 'chip2_peaks.bed.passed')
    header = ['chrom', 'start', 'end', 'name', 'foldChange',
              'strand', 'signalValue', 'pvalue', 'qvalue']
    recs = []
    with open(peakfile, 'r', newline='') as peaks:
        rows = csv.DictReader(peaks, delimiter='\t', fieldnames=header)
        for row in rows:
            if float(row['qvalue']) < 0.01:
                recs.append(row)
    recs = sorted(recs, key=lambda d: (d['chrom'], int(d['start']), int(d['end'])))
    peak_lens = sorted([int(r['end']) - int(r['start']) for r in recs])
    total_peaks = len(peak_lens)
    shortest = peak_lens[0]
    longest = peak_lens[-1]
    pct10 = int(total_peaks * 0.1)
    pct25 = int(total_peaks * 0.25)
    pct10_len = peak_lens[:pct10][-1]
    pct25_len = peak_lens[:pct25][-1]
    with open(outputfile, 'w', encoding='utf-8', newline='') as filtered:
        _ = filtered.write('#')
        writer = csv.DictWriter(filtered, delimiter='\t',
                                fieldnames=header, dialect='unix',
                                quoting=csv.QUOTE_NONE)
        writer.writeheader()
        writer.writerows(recs)
    statfile = outputfile.replace('sig.bed', 'stat.txt')
    with open(statfile, 'w') as stats:
        for k, v in zip(['total_num', 'shortest_peak', 'longest_peak',
                         'percentile_10', 'percentile_25'],
                        [total_peaks, shortest, longest, pct10_len, pct25_len]):
            _ = stats.write(k + '\t' + str(v) + '\n')
    return outputfile


def make_sciddo_scan_params(inputfolder, sub_folder, groups, scorings, cmd, jobcall):
    """
    :param inputfolder:
    :param sub_folder:
    :param groups:
    :param scorings:
    :param cmd:
    :param jobcall:
    :return:
    """
    args = []
    for root, dirs, data in os.walk(inputfolder):
        if root.endswith('_run') or root.endswith('data_dump'):
            continue
        if data:
            for fn in data:
                if fn.endswith('.h5') and fn.startswith('sciddo-data_'):
                    in_path = os.path.join(root, fn)
                    out_dir = os.path.join(root, sub_folder)
                    os.makedirs(out_dir, exist_ok=True)
                    outfile = fn.replace('sciddo-data_', 'sciddo-run_')
                    for g1, g2 in itt.combinations(groups, 2):
                        for s in scorings:
                            c1 = g1.split('_')[1]
                            c2 = g2.split('_')[1]
                            tmp = cmd.format(**{'grp1': g1, 'grp2': g2, 'scoring': s})
                            out_tmp = outfile.replace('.h5', '_' + s + '_' + c1 + '_vs_' + c2 + '.h5')
                            out_path = os.path.join(out_dir, out_tmp)
                            args.append([in_path, out_path, tmp, jobcall])
                    for s in scorings:
                        g1 = 'TISSUE_Li'
                        g2 = 'TISSUE_Bl'
                        tmp = cmd.format(**{'grp1': g1, 'grp2': g2, 'scoring': s})
                        out_tmp = outfile.replace('.h5', '_' + s + '_Liver_vs_Blood.h5')
                        out_path = os.path.join(out_dir, out_tmp)
                        args.append([in_path, out_path, tmp, jobcall])
    return args


def make_sciddo_dump_params(inputroot, sub_output, scorings, hsp_path, cmd, jobcall):
    """
    :param inputroot:
    :param sub_output:
    :param scorings:
    :param cmd:
    :param jobcall:
    :return:
    """
    args = []
    extract_infos = op.itemgetter(*tuple([0, 1, 2, 3, 5, 6, 8]))
    for root, dirs, sciddo_files in os.walk(inputroot):
        if files and root.endswith('/hsp_run'):
            for f in sciddo_files:
                fn_parts = f.split('.')[0].split('_')
                outbase = '_'.join(extract_infos(fn_parts))
                for s in scorings:
                    tmp_out = outbase + '_' + s + '.bed'
                    tmp_out = os.path.join(root, sub_output, tmp_out)
                    tmp_cmd = cmd.format(**{'dump_path': os.path.join(hsp_path, s)})
                    args.append([os.path.join(root, f), tmp_out, tmp_cmd, jobcall])
    return args


def make_btools_all_intersect(inputroot, outroot, datatype, cmd, jobcall):
    """
    :param inputroot:
    :param outroot:
    :param cmd:
    :param jobcall:
    :return:
    """
    # subfolder hierarchy:
    # 1) identical/different segmentation
    # 2) identical/different sample comparison
    # 3) identical/different scoring scheme
    subfolders = {(1, 1, 1): 'ident_ident_ident', (1, 0, 1): 'ident_diff_ident',
                  (0, 0, 0): 'diff_diff_diff', (0, 0, 1): 'diff_diff_ident',
                  (0, 1, 1): 'diff_ident_ident', (1, 0, 0): 'ident_diff_diff',
                  (1, 1, 0): 'ident_ident_diff', (0, 1, 0): 'diff_ident_diff'}
    for folder in subfolders.values():
        os.makedirs(os.path.join(outroot, folder), exist_ok=True)
    # cmm18_hg38_He_vs_Ma_ordem_raw-t1.bed
    match_info = re.compile("(?P<SEGMENT>[a-z0-9]+)_hg38_(?P<COMP>\w+)_(?P<SCORING>[pordemn]+)_(?P<DATA>(hsp|raw)\-t[10]+)\.bed")
    inputfiles = []
    for root, dirs, datafiles in os.walk(inputroot):
        if datafiles and 'data_dump' in root:
            for bf in datafiles:
                mobj = match_info.match(bf)
                if mobj is not None:
                    if mobj.group('DATA') == datatype:
                        inputfiles.append((os.path.join(root, bf), mobj))
    args = []
    for b1, b2 in itt.permutations(sorted(inputfiles), 2):
        b1_path = b1[0]
        b2_path = b2[0]
        b1_info = b1[1]
        b2_info = b2[1]
        outname = '{}-{}-{}-ovl-{}-{}-{}_{}.tsv'.format(b1_info.group('SEGMENT'),
                                                        b1_info.group('COMP'),
                                                        b1_info.group('SCORING'),
                                                        b2_info.group('SEGMENT'),
                                                        b2_info.group('COMP'),
                                                        b2_info.group('SCORING'),
                                                        datatype)

        a = int(b1_info.group('SEGMENT') == b2_info.group('SEGMENT'))
        b = int(b1_info.group('COMP') == b2_info.group('COMP'))
        c = int(b1_info.group('SCORING') == b2_info.group('SCORING'))
        if c == 0:
            # skip over different scorings
            continue
        seg1, seg2 = b1_info.group('SEGMENT'), b2_info.group('SEGMENT')
        if (seg1, seg2) in [('cmm18', 'ecs10'), ('ecs10', 'cmm18')]:
            continue
        subfolder = subfolders[(a, b, c)]
        outpath = os.path.join(outroot, subfolder, outname)
        args.append([[b1_path, b2_path], outpath, cmd, jobcall])
    return args


def make_btools_hsp_pepr_isect(hsp_root, pepr_root, deseq_root, regtype, condition, outroot, cmd, jobcall):
    """
    :return:
    """
    args = []
    genes_raw = os.path.join(deseq_root, 'deseq2_{}_{}_{}.bed')
    for root, dirs, dumpfile in os.walk(hsp_root):
        if root.endswith('data_dump'):
            hspfiles = [f for f in dumpfile if '_hsp-' in f and f.endswith('.bed')]
            for hsp in hspfiles:
                peakfiles = []
                # deep_hg38_cmm18_HG_vs_Mo_hsp-emission.bed
                parts = hsp.split('.')[0].split('_')
                seg = parts[2]
                comp = parts[3] + '_vs_' + parts[5]
                score = parts[-1].split('-')[-1]
                chip_lut = {'chip1': parts[3], 'chip2': parts[5]}
                pepr_folder = os.path.join(pepr_root, comp)
                pepr_files = os.listdir(pepr_folder)
                pepr_files = [f for f in pepr_files if f.endswith('peak.sig.bed')]
                for pf in pepr_files:
                    id_parts = pf.split('_')
                    mark = id_parts[3]
                    cell = chip_lut[id_parts[6]]
                    file_label = cell + '-' + mark
                    file_path = os.path.join(pepr_folder, pf)
                    peakfiles.append((file_label, file_path))
                peakfiles = sorted(peakfiles)
                isect_labels = ['HSP'] + [t[0] for t in peakfiles]
                isect_paths = [os.path.join(root, hsp)] + [t[1] for t in peakfiles]
                genefile = genes_raw.format(comp, condition, regtype)
                outfile = '{}-genes_{}_hsp_pepr_{}_{}_{}.tsv'.format(condition, regtype, seg, comp, score)
                outpath = os.path.join(outroot, outfile)
                tmp = cmd.format(**{'isect_files': ' '.join(isect_paths),
                                    'isect_labels': ' '.join(isect_labels)})
                args.append([genefile, outpath, tmp, jobcall])
    return args


def make_invert_intersect_params(gene_root, enhancer_file, sciddo_root, out_folder, cmd, jobcall):
    """
    Inverse intersection - intersect genes/enhancers with HSPs
    to see how HSP hits alter gene expression patterns

    :param gene_root: folder with DeSeq2 output files
    :param enhancer_file:
    :param sciddo_root: root for HSP dumps
    :param out_folder:
    :param cmd:
    :param jobcall:
    :return:
    """
    args = []
    for gene_file in os.listdir(gene_root):
        if not gene_file.endswith('_body.bed'):
            continue
        gene_path = os.path.join(gene_root, gene_file)
        parts = gene_file.split('_')
        comp = '_vs_'.join([parts[1], parts[3]])
        gene_type = parts[4]

        for root, dirs, data_files in os.walk(sciddo_root):
            if root.endswith('data_dump'):
                hsp_files = [f for f in data_files if '_hsp-' in f and f.endswith('.bed') and comp in f]
                for f in hsp_files:
                    hsp_path = os.path.join(root, f)
                    seg = f.split('_')[2]
                    scoring = f.split('.')[0].split('-')[-1]
                    gene_prefix = 'st'
                    if gene_type == 'diff':
                        gene_prefix = 'de'
                        infile1 = enhancer_file
                        infile2 = hsp_path
                        outfile = 'deep_enh_ovl_hsp_{}_{}_{}.tsv'.format(seg, comp, scoring)
                        outpath = os.path.join(out_folder, outfile)
                        tmp = str(cmd)
                        args.append([(infile1, infile2), outpath, tmp, jobcall])
                    infile1 = gene_path
                    infile2 = hsp_path
                    outfile = 'deep_{}genes_ovl_hsp_{}_{}_{}.tsv'.format(gene_prefix, seg, comp, scoring)
                    outpath = os.path.join(out_folder, outfile)
                    tmp = str(cmd)
                    args.append([(infile1, infile2), outpath, tmp, jobcall])

    return args


def make_pepr_any_intersect_params(gene_root, enhancer_file, pepr_root, out_folder, cmd, jobcall):
    """
    Inverse intersection - intersect genes/enhancers with HSPs
    to see how HSP hits alter gene expression patterns

    :param gene_root: folder with DeSeq2 output files
    :param enhancer_file:
    :param pepr_root: root for PePr peak files
    :param out_folder:
    :param cmd:
    :param jobcall:
    :return:
    """
    args = []
    use_marks = ['H3K27ac', 'H3K4me3', 'H3K36me3']
    for gene_file in os.listdir(gene_root):
        if gene_file.endswith('locus.bed'):
            continue
        if '_all_' not in gene_file:
            continue
        gene_path = os.path.join(gene_root, gene_file)
        region = gene_file.split('.')[0].split('_')[-1]

        parts = gene_file.split('_')
        s1, s2 = parts[1], parts[3]
        sample_lut = {'chip1': s1, 'chip2': s2}
        comp = '_vs_'.join([s1, s2])

        done_peaks = set()
        for root, dirs, data_files in os.walk(pepr_root):
            if root.endswith(comp):
                peak_files = [f for f in data_files if f.endswith('.sig.bed')]
                for p in peak_files:
                    mark = p.split('_')[3]
                    if mark not in use_marks:
                        continue
                    peak_path = os.path.join(root, p)

                    condition = p.split('_')[5]
                    if mark == 'H3K27ac':
                        infile1 = enhancer_file
                        infile2 = peak_path
                        outfile = 'enh_ovl_{}-{}_{}.tsv'.format(sample_lut[condition], mark, comp)
                        if outfile not in done_peaks:
                            outpath = os.path.join(out_folder, outfile)
                            tmp = str(cmd)
                            args.append([(infile1, infile2), outpath, tmp, jobcall])
                    elif mark == 'H3K4me3' and region == 'promoter':
                        infile1 = gene_path
                        infile2 = peak_path
                        outfile = 'prom_ovl_{}-{}_{}.tsv'.format(sample_lut[condition], mark, comp)
                        if outfile not in done_peaks:
                            outpath = os.path.join(out_folder, outfile)
                            tmp = str(cmd)
                            args.append([(infile1, infile2), outpath, tmp, jobcall])
                    elif mark == 'H3K36me3' and region == 'body':
                        infile1 = gene_path
                        infile2 = peak_path
                        outfile = 'gene_ovl_{}-{}_{}.tsv'.format(sample_lut[condition], mark, comp)
                        if outfile not in done_peaks:
                            outpath = os.path.join(out_folder, outfile)
                            tmp = str(cmd)
                            args.append([(infile1, infile2), outpath, tmp, jobcall])
                    else:
                        pass
    return args


def collect_random_set_files(workdir):

    folders = [
        'sciddo/deep/cmm18/data_dump/rand_t1',
        'sciddo/deep/ecs18/data_dump/rand_t1',
        'sciddo/deep/ecs10/data_dump/rand_t1'
    ]

    random_sets = []
    for f in folders:
        bed_files = os.listdir(os.path.join(workdir, f))
        bed_files = list(filter(lambda x: x.endswith('.bed'), bed_files))
        bed_files = [os.path.join(workdir, f, bf) for bf in bed_files]
        random_sets.extend(bed_files)
    return random_sets


def touch_checkfile(inputfiles, outputfile):

    with open(outputfile, 'w') as checkfile:
        _ = checkfile.write(ti.ctime() + '\n\n')
        _ = checkfile.write('\n'.join(sorted(inputfiles)) + '\n')
    return outputfile


def build_pipeline(args, config, sci_obj, pipe):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """
    if pipe is None:
        pipe = Pipeline(name=config.get('Pipeline', 'name'))
    else:
        # remove all previous tasks from pipeline
        pipe.clear()
        # turns out clear() seems NOT to have the effect
        # of really clearing the pipeline object, do it manually...
        pipe.task_names = set()
        pipe.tasks = set()

    workdir = config.get('EnvPaths', 'workdir')

    # ============================================================
    # Task configuration

    linked_input = os.path.join(workdir, 'linked_input')
    loaded_input = os.path.join(workdir, 'loaded_input')

    deep_input = os.path.join(linked_input, 'deep')
    deep_folder = config.get('Samples', 'deep_folder')
    deep_samples = config.get('Samples', 'deep_samples').split()
    histone_libs = config.get('Samples', 'libraries').split()
    assume_complete = config.getboolean('Pipeline', 'assume_complete')
    init_deep_samples = pipe.originate(lambda x: x,
                                       link_deep_samples(deep_folder, deep_samples,
                                                         histone_libs, deep_input,
                                                         assume_complete),
                                       name='init_deep_samples')
    init_deep_samples = init_deep_samples.mkdir(linked_input)
    init_deep_samples = init_deep_samples.mkdir(deep_input)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    deep_filtered_output = os.path.join(workdir, 'filtered_input', 'deep')
    cmd = config.get('Pipeline', 'bam_qfilter')
    bam_qfilter_deep = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='bam_qfilter_deep',
                                      input=output_from(init_deep_samples),
                                      filter=suffix('.bam'),
                                      output='.filt.bam',
                                      output_dir=deep_filtered_output,
                                      extras=[cmd, jobcall])
    bam_qfilter_deep = bam_qfilter_deep.mkdir(deep_filtered_output)

    # sort BAM files by read name ("qname"): expected by PePr

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    deep_sorted_output = os.path.join(workdir, 'sorted_input', 'deep', 'qname_sort')
    cmd = config.get('Pipeline', 'bam_qname_sort')
    subset_re = '(?P<SAMPLE>(01_Hc|41_Hf|43_Hm)[0-9]+_(Bl|Li)(HG|He|Mo|Ma)_Ct\w+)\.filt\.bam$'
    bam_qname_sort_deep = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                         name='bam_qname_sort_deep',
                                         input=output_from(bam_qfilter_deep),
                                         filter=formatter(subset_re),
                                         output=os.path.join(deep_sorted_output, '{SAMPLE[0]}.qnsrt.bam'),
                                         extras=[cmd, jobcall])
    bam_qname_sort_deep = bam_qname_sort_deep.mkdir(deep_sorted_output)

    # =====================================================================
    # Chromatine state segmentation with ChromHMM
    # Set folders for DEEP samples
    cmm_deep_folder = os.path.join(workdir, 'chromhmm', 'deep')
    cmm_deep_binarized = os.path.join(cmm_deep_folder, 'binarized')
    cmm_deep_segmentation = os.path.join(cmm_deep_folder, 'segmentation')

    # Create ChromHMM DEEP sample sheet
    cmm_deep_sheet = pipe.merge(task_func=prep_chromhmm_table,
                                name='cmm_deep_sheet',
                                input=output_from(bam_qfilter_deep),
                                output=os.path.join(cmm_deep_folder, 'cmm_deep_samples.tsv'))
    cmm_deep_sheet = cmm_deep_sheet.mkdir(cmm_deep_folder)

    # Binarize all BAM files
    cmd_binbam_deep = config.get('Pipeline', 'cmm_binbam')
    tmp_binbam_deep = cmd_binbam_deep.format(**{'inputdir': deep_filtered_output,
                                                'outputdir': cmm_deep_binarized})
    cmm_deep_binbam = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                                     name='cmm_deep_binbam',
                                     input=output_from(cmm_deep_sheet),
                                     filter=formatter('cmm_deep_samples\.tsv'),
                                     output=os.path.join(cmm_deep_binarized, '*.txt'),
                                     extras=[cmm_deep_binarized, '*.txt', tmp_binbam_deep, jobcall])
    cmm_deep_binbam = cmm_deep_binbam.mkdir(cmm_deep_binarized)
    cmm_deep_binbam = cmm_deep_binbam.active_if(check_run_binbam(os.path.join(cmm_deep_folder, 'cmm_deep_samples.tsv'),
                                                                 cmm_deep_binarized))

    cmd_mkseg_deep = config.get('Pipeline', 'cmm_mkseg')
    tmp_mkseg_deep = cmd_mkseg_deep.format(**{'inputdir': cmm_deep_binarized,
                                              'outputdir': cmm_deep_segmentation})
    cmm_deep_mkseg = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                                name='cmm_deep_mkseg',
                                input=output_from(cmm_deep_binbam),
                                output=os.path.join(cmm_deep_folder, 'cmm_deep_seg.chk'),
                                extras=[tmp_mkseg_deep, jobcall])
    cmm_deep_mkseg = cmm_deep_mkseg.mkdir(cmm_deep_segmentation)

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # learn a new ChromHMM 18-state model just on the DEEP data to potentially
    # check for differences in the segmentations and state emissions

    cmd = config.get('Pipeline', 'cmm_learn')
    cmm_deep_learn = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                    name='cmm_deep_learn',
                                    input=output_from(cmm_deep_sheet),
                                    filter=formatter(),
                                    output='{path[0]}/cmm_deep_learn.chk',
                                    extras=[cmd, jobcall])
    cmm_deep_learn = cmm_deep_learn.follows(cmm_deep_binbam)
    cmm_deep_learn = cmm_deep_learn.mkdir(cmm_deep_folder, 'test_new18')

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # With EpiCSeg - 18 and 10 states
    ecs_deep_folder = os.path.join(workdir, 'epicseg', 'deep')
    ecs_deep_count_sub = os.path.join(ecs_deep_folder, 'count_mat')
    cmd = config.get('Pipeline', 'ecs_count')
    ecs_deep_count_params = make_ecs_count_params(deep_filtered_output, config.get('Pipeline', 'deep_subset'),
                                                  ecs_deep_count_sub, cmd, jobcall)
    ecs_deep_count = pipe.files(sci_obj.get_jobf('ins_out'),
                                ecs_deep_count_params,
                                name='ecs_deep_count')
    ecs_deep_count = ecs_deep_count.mkdir(ecs_deep_count_sub)
    ecs_deep_count = ecs_deep_count.follows(bam_qfilter_deep)

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'ecs_norm')
    ecs_deep_norm_params = make_ecs_norm_params(ecs_deep_count_sub, cmd, jobcall)
    ecs_deep_norm = pipe.files(sci_obj.get_jobf('ins_out'),
                               ecs_deep_norm_params,
                               name='ecs_deep_norm')
    ecs_deep_norm = ecs_deep_norm.follows(ecs_deep_count)

    cmd = config.get('Pipeline', 'ecs_seg')
    ecs_deep_seg_params = make_ecs_seg_params(ecs_deep_count_sub, ecs_deep_folder, ['18', '10'],
                                              cmd, jobcall)
    ecs_deep_seg = pipe.files(sci_obj.get_jobf('ins_out'),
                              ecs_deep_seg_params,
                              name='ecs_deep_seg')
    ecs_deep_seg = ecs_deep_seg.follows(ecs_deep_norm)

    # =====================================================================
    # Peak calling with replicates
    # Set folders for DEEP samples

    zerone_deep_folder = os.path.join(workdir, 'zerone', 'deep')
    zerone_deep_rawout = os.path.join(zerone_deep_folder, 'raw_output')

    cmd = config.get('Pipeline', 'pk_zerone')
    zerone_deep_params = make_zerone_deep_params(deep_filtered_output,
                                                 zerone_deep_rawout,
                                                 cmd, jobcall)

    zerone_deep_peaks = pipe.files(sci_obj.get_jobf('ins_out'),
                                   zerone_deep_params,
                                   name='zerone_deep_peaks')
    zerone_deep_peaks = zerone_deep_peaks.mkdir(zerone_deep_rawout)
    zerone_deep_peaks = zerone_deep_peaks.follows(bam_qfilter_deep)

    # =====================================================================
    # Differential peak calling with replicates (PePr)
    # Set folder for DEEP samples
    pepr_deep_folder = os.path.join(workdir, 'pepr', 'deep')
    pepr_deep_output = os.path.join(pepr_deep_folder, 'diff_out')

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'pkdiff_pepr')
    pepr_deep_diff_params = make_pepr_diff_params(deep_sorted_output, pepr_deep_output,
                                                  config.get('Pipeline', 'deep_subset'), cmd, jobcall)
    pepr_deep_diff = pipe.files(sci_obj.get_jobf('ins_out'),
                                pepr_deep_diff_params,
                                name='pepr_deep_diff')
    pepr_deep_diff = pepr_deep_diff.mkdir(pepr_deep_folder)
    pepr_deep_diff = pepr_deep_diff.follows(bam_qname_sort_deep)

    cmd = config.get('Pipeline', 'pkpost_pepr')
    pepr_deep_post_params = make_pepr_post_params(pepr_deep_output, cmd, jobcall)

    pepr_deep_post = pipe.files(sci_obj.get_jobf('in_out'),
                                pepr_deep_post_params,
                                name='pepr_deep_post')
    pepr_deep_post = pepr_deep_post.follows(pepr_deep_diff)

    cmd = config.get('Pipeline', 'pkuniq12_pepr')
    pepr_peak12_uniq = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='pepr_peak12_uniq',
                                      input=output_from(pepr_deep_post),
                                      filter=formatter('(?P<SAMPLEID>\w+)__PePr_post\-peak1\.chk'),
                                      output=os.path.join('{path[0]}', '{SAMPLEID[0]}_PePr_chip1_peaks.uniq.bed'),
                                      extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'pkuniq21_pepr')
    pepr_peak21_uniq = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='pepr_peak21_uniq',
                                      input=output_from(pepr_deep_post),
                                      filter=formatter('(?P<SAMPLEID>\w+)__PePr_post\-peak2\.chk'),
                                      output=os.path.join('{path[0]}', '{SAMPLEID[0]}_PePr_chip2_peaks.uniq.bed'),
                                      extras=[cmd, jobcall])

    # filter "passed" BED files for q-value
    # Mo_vs_Ma_H3K4me1__PePr_post-peak1.chk
    pepr_deep_filter = pipe.transform(task_func=filter_sort_pepr_peaks,
                                      name='pepr_deep_filter',
                                      input=output_from(pepr_peak12_uniq, pepr_peak21_uniq),
                                      filter=formatter('(?P<SAMPLEID>\w+)_PePr_chip(?P<REP>(1|2))_peaks\.uniq\.bed'),
                                      output=os.path.join('{path[0]}', '{SAMPLEID[0]}_PePr_chip{REP[0]}_peak.sig.bed'))
    pepr_deep_filter = pepr_deep_filter.jobs_limit(4)
    pepr_deep_filter = pepr_deep_filter.follows(pepr_peak12_uniq)
    pepr_deep_filter = pepr_deep_filter.follows(pepr_peak21_uniq)
    
    # =====================================================================
    # Differential peak calling with replicates (THOR)
    # Set folder for DEEP samples
    thor_deep_folder = os.path.join(workdir, 'thor', 'deep')
    thor_deep_output = os.path.join(thor_deep_folder, 'diff_out')

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('Py2EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'pkdiff_thor')
    thor_deep_diff_params = make_thor_diff_params(deep_filtered_output, thor_deep_output,
                                                  config.get('Pipeline', 'deep_subset'), cmd, jobcall)
    thor_deep_diff = pipe.files(sci_obj.get_jobf('in_out'),
                                thor_deep_diff_params,
                                name='thor_deep_diff')
    thor_deep_diff = thor_deep_diff.mkdir(thor_deep_folder)
    thor_deep_diff = thor_deep_diff.follows(bam_qfilter_deep)
    # there is still a bug in THOR that prevents execution
    thor_deep_diff = thor_deep_diff.active_if(False)

    # =====================================================================
    # Gene expression quantification
    # Set folders for DEEP samples

    salmon_deep_folder = os.path.join(workdir, 'salmon', 'deep')
    salmon_deep_quant_sub = os.path.join(salmon_deep_folder, 'quant')

    deep_exp_fastq = os.path.join(loaded_input, 'deep', 'rna_data', 'sequence')

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'salmon_deep_quant')
    salmon_deep_quant_params = make_salmon_calls(deep_exp_fastq, salmon_deep_quant_sub, cmd, jobcall)

    salmon_deep_quant = pipe.files(sci_obj.get_jobf('ins_out'),
                                   salmon_deep_quant_params,
                                   name='salmon_deep_quant')
    salmon_deep_quant = salmon_deep_quant.mkdir(salmon_deep_quant_sub)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    deseq_deep_folder = os.path.join(workdir, 'deseq', 'deep')

    cmd = config.get('Pipeline', 'deseq_deep_diff')
    tmp = cmd.format(**{'outpath': deseq_deep_folder,
                        'quantpath': salmon_deep_quant_sub})
    deseq_deep_diff = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                                 name='deseq_deep_diff',
                                 input=output_from(salmon_deep_quant),
                                 output=os.path.join(workdir, 'deseq', 'deseq_deep_diff.chk'),
                                 extras=[tmp, jobcall])
    deseq_deep_diff = deseq_deep_diff.mkdir(deseq_deep_folder)
    deseq_deep_diff = deseq_deep_diff.follows(salmon_deep_quant)

    deep_to_bed = os.path.join(workdir, 'deseq', 'bed_out')
    cmd = config.get('Pipeline', 'conv_deep_diff')
    conv_deep_diff = pipe.subdivide(task_func=sci_obj.get_jobf('in_outpair'),
                                    name='conv_deep_diff',
                                    input=[os.path.join(deseq_deep_folder, fn) for fn in os.listdir(deseq_deep_folder)],
                                    filter=formatter(),
                                    output=[os.path.join(deep_to_bed, '{basename[0]}_diff_body.bed'),
                                            os.path.join(deep_to_bed, '{basename[0]}_diff_promoter.bed')],
                                            # to keep output as pair, ignore that gene loci are also produced
                                            #os.path.join(deep_to_bed, '{basename[0]}_diff_locus.bed')],
                                    extras=[cmd, jobcall])
    conv_deep_diff = conv_deep_diff.mkdir(deep_to_bed)
    conv_deep_diff = conv_deep_diff.follows(deseq_deep_diff)

    # =======================================================
    # Process state segmentation of DEEP samples with SCIDDO
    # Set folders for DEEP samples
    # =======================================================

    # NB: create sub-folders for each segmentation
    sciddo_deep_folder = os.path.join(workdir, 'sciddo', 'deep')

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'sciddo_deep_conv_cmm18')
    sciddo_deep_conv_cmm18 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='sciddo_deep_conv_cmm18',
                                            input=output_from(cmm_deep_mkseg),
                                            filter=formatter('chk$'),
                                            output=os.path.join(sciddo_deep_folder, 'cmm18',
                                                                'sciddo-data_hg38_cmm18.conv-chk'),
                                            extras=[cmd, jobcall])
    sciddo_deep_conv_cmm18 = sciddo_deep_conv_cmm18.mkdir(sciddo_deep_folder)

    # convert EpiCSeg
    cmd = config.get('Pipeline', 'sciddo_deep_conv_ecs')
    sciddo_deep_conv_ecs = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='sciddo_deep_conv_ecs',
                                          input=output_from(ecs_deep_seg),
                                          filter=formatter('ecs_deep_seg_(?P<STNUM>[0-9]+)\.chk'),
                                          output=os.path.join(sciddo_deep_folder, 'ecs{STNUM[0]}',
                                                              'sciddo-data_hg38_ecs{STNUM[0]}.conv-chk'),
                                          extras=[cmd, jobcall])
    sciddo_deep_conv_ecs = sciddo_deep_conv_ecs.mkdir(sciddo_deep_folder)

    # compute counting statistics - identical for all
    cmd = config.get('Pipeline', 'sciddo_deep_stats')
    sciddo_deep_stats = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='sciddo_deep_stats',
                                       input=output_from(sciddo_deep_conv_cmm18, sciddo_deep_conv_ecs),
                                       filter=formatter('conv\-chk$'),
                                       output='{path[0]}/{basename[0]}.stat-chk',
                                       extras=[cmd, jobcall])

    # add scoring matrices - identical for all
    cmd = config.get('Pipeline', 'sciddo_deep_score')
    sciddo_deep_score = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='sciddo_deep_score',
                                       input=output_from(sciddo_deep_stats),
                                       filter=formatter('stat\-chk$'),
                                       output='{path[0]}/{basename[0]}.score-chk',
                                       extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # ==============================
    # make a baseline run
    # that compares all replicates
    # in a 1-vs-1 manner against
    # each other

    cmd = config.get('Pipeline', 'sciddo_deep_base')
    sciddo_deep_base = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='sciddo_deep_base',
                                      input=output_from(sciddo_deep_score),
                                      filter=formatter('score\-chk$'),
                                      output=os.path.join('{path[0]}', 'baseline_run', '{basename[0]}_baseline.h5'),
                                      extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_smprand')
    sciddo_deep_smprand = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                         name='sciddo_deep_smprand',
                                         input=output_from(sciddo_deep_score),
                                         filter=formatter('score\-chk$'),
                                         output=os.path.join('{path[0]}', 'baseline_run', '{basename[0]}_smprand.h5'),
                                         extras=[cmd, jobcall])

    # ==============================
    # Scan for differential regions
    # based on cell-type grouping
    cmd = config.get('Pipeline', 'sciddo_deep_scan')
    groups = config.get('Pipeline', 'sciddo_deep_groups').split()
    scorings = config.get('Pipeline', 'use_scorings').split()
    sciddo_deep_scan_params = make_sciddo_scan_params(sciddo_deep_folder, 'hsp_run',
                                                      groups, scorings, cmd, jobcall)
    sciddo_deep_scan = pipe.files(sci_obj.get_jobf('in_out'),
                                  sciddo_deep_scan_params,
                                  name='sciddo_deep_scan')
    sciddo_deep_scan = sciddo_deep_scan.follows(sciddo_deep_score)
    # ================================================================================
    # Dump data for visualization
    # and HSPs

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # dump fully annotated state segmentations
    deep_dset_re = 'sciddo-data_hg38_(?P<SEGMENT>[cmes180]{5})\.conv-chk'
    cmd = config.get('Pipeline', 'sciddo_deep_dump_seg')
    sciddo_deep_dump_seg = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                                          name='sciddo_deep_dump_seg',
                                          input=output_from(sciddo_deep_conv_cmm18, sciddo_deep_conv_ecs),
                                          filter=formatter(deep_dset_re),
                                          output='{path[0]}/data_dump/state_maps/{SEGMENT[0]}_*.bed',
                                          extras=['{path[0]}/data_dump/state_maps', '{SEGMENT[0]}_*.bed', cmd, jobcall])

    # Regular expression to match SCIDDO run files
    deep_scan_re = 'sciddo-run_hg38_(?P<SEGMENT>[cmes180]{5})_(?P<SCORING>[a-z]+)_(?P<C1>[HGeMoaLivr]{2,5})_vs_(?P<C2>[HGeMoaBld]{2,5})\.h5'

    cmd = config.get('Pipeline', 'sciddo_deep_dump_scores')
    sciddo_deep_dump_scores = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                             name='sciddo_deep_dump_scores',
                                             input=output_from(sciddo_deep_scan),
                                             filter=formatter(deep_scan_re),
                                             output='{subpath[0][1]}/data_dump/score_tracks/{SEGMENT[0]}_hg38_{C1[0]}_vs_{C2[0]}_{SCORING[0]}_scores.bg',
                                             extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_trans')
    sciddo_deep_dump_trans = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='sciddo_deep_dump_trans',
                                            input=output_from(sciddo_deep_scan),
                                            filter=formatter(deep_scan_re),
                                            output='{subpath[0][1]}/data_dump/transition_counts/{SEGMENT[0]}_hg38_{C1[0]}_vs_{C2[0]}_{SCORING[0]}_transitions.tsv',
                                            extras=[cmd, jobcall])

    deep_scan_re_hema = 'sciddo-run_hg38_(?P<SEGMENT>[cm18]{5})_(?P<SCORING>[a-z]+)_(?P<C1>[He]{2,5})_vs_(?P<C2>[Ma]{2,5})\.h5'
    cmd = config.get('Pipeline', 'sciddo_deep_dump_scores_hema')
    sciddo_deep_dump_scores_hema = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                                                  name='sciddo_deep_dump_scores_hema',
                                                  input=output_from(sciddo_deep_scan),
                                                  filter=formatter(deep_scan_re_hema),
                                                  output='{subpath[0][1]}/data_dump/pw_score_tracks/*-scores.bg',
                                                  extras=['{subpath[0][1]}/data_dump/pw_score_tracks', '*.bg', cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_hsp_t1')
    sciddo_deep_dump_hsp_t1 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                             name='sciddo_deep_dump_hsp_t1',
                                             input=output_from(sciddo_deep_scan),
                                             filter=formatter(deep_scan_re),
                                             output='{subpath[0][1]}/data_dump/hsp_t1/{SEGMENT[0]}_hg38_{C1[0]}_vs_{C2[0]}_{SCORING[0]}_hsp-t1.bed',
                                             extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_hsp_t100')
    sciddo_deep_dump_hsp_t100 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                               name='sciddo_deep_dump_hsp_t100',
                                               input=output_from(sciddo_deep_scan),
                                               filter=formatter(deep_scan_re),
                                               output='{subpath[0][1]}/data_dump/hsp_t100/{SEGMENT[0]}_hg38_{C1[0]}_vs_{C2[0]}_{SCORING[0]}_hsp-t100.bed',
                                               extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_raw_t1')
    sciddo_deep_dump_raw_t1 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                             name='sciddo_deep_dump_raw_t1',
                                             input=output_from(sciddo_deep_scan),
                                             filter=formatter(deep_scan_re),
                                             output='{subpath[0][1]}/data_dump/raw_t1/{SEGMENT[0]}_hg38_{C1[0]}_vs_{C2[0]}_{SCORING[0]}_raw-t1.bed',
                                             extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # convert score tracks to bigWig format
    cmd = config.get('Pipeline', 'bg_to_bw')
    bg_to_bw_scores = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                     name='bg_to_bw_scores',
                                     input=output_from(sciddo_deep_dump_scores, sciddo_deep_dump_scores_hema),
                                     filter=suffix('.bg'),
                                     output='.bw',
                                     extras=[cmd, jobcall])

    # root folder for all intersections
    btl_deep_folder = os.path.join(workdir, 'bedtools', 'deep')

    match_hsp_bed = '(?P<SEG>[a-z0-9]+)_hg38_(?P<COMP>\w+)_(?P<SCORING>[a-z]+)_hsp\-(?P<THRESHOLD>t[10]+)\.bed'

    # ======= ORIENTATION =======
    # INTERSECT HSPs WITH FILE X
    # ===========================

    # intersect HSPs with various annotation files
    cmd = config.get('Pipeline', 'btl_isect_hsp_any')
    out_folder = os.path.join(btl_deep_folder, 'isect_hsp_any')
    out_name = 'hsp_ovl_reg_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_isect_hsp_any = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='btl_isect_hsp_any',
                                       input=output_from(sciddo_deep_dump_hsp_t1),
                                       filter=formatter(match_hsp_bed),
                                       output=os.path.join(out_folder, out_name),
                                       extras=[cmd, jobcall])
    btl_isect_hsp_any = btl_isect_hsp_any.mkdir(out_folder)

    # intersect HSPs with Ensembl Regulatory Build for more detailed analysis
    cmd = config.get('Pipeline', 'btl_isect_hsp_rgb')
    out_folder = os.path.join(btl_deep_folder, 'isect_hsp_rgb')
    out_name = 'hsp_ovl_rgb_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_isect_hsp_rgb = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='btl_isect_hsp_rgb',
                                       input=output_from(sciddo_deep_dump_hsp_t1),
                                       filter=formatter(match_hsp_bed),
                                       output=os.path.join(out_folder, out_name),
                                       extras=[cmd, jobcall])
    btl_isect_hsp_rgb = btl_isect_hsp_rgb.mkdir(out_folder)

    # intersect HSPs with gene for more detailed analysis (similar to RGB boxplot)
    cmd = config.get('Pipeline', 'btl_isect_hsp_gene')
    out_folder = os.path.join(btl_deep_folder, 'isect_hsp_gene')
    out_name = 'hsp_ovl_gene_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_isect_hsp_gene = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_isect_hsp_gene',
                                        input=output_from(sciddo_deep_dump_hsp_t1),
                                        filter=formatter(match_hsp_bed),
                                        output=os.path.join(out_folder, out_name),
                                        extras=[cmd, jobcall])
    btl_isect_hsp_gene = btl_isect_hsp_gene.mkdir(out_folder)

    # to check if there is a problem detecting shorter DE genes
    # intersect genes also with EXP=100 thresholded HSPs
    cmd = config.get('Pipeline', 'btl_isect_gene_hsp')
    out_folder = os.path.join(btl_deep_folder, 'isect_gene_hsp100')
    out_name = 'gene_ovl_hsp-t100_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_deep_isect_t100 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                         name='btl_isect_gene_hsp100',
                                         input=output_from(sciddo_deep_dump_hsp_t100),
                                         filter=formatter(match_hsp_bed),
                                         output=os.path.join(out_folder, out_name),
                                         extras=[cmd, jobcall])
    btl_deep_isect_t100 = btl_deep_isect_t100.mkdir(out_folder)

    # intersect raw (unmerged) segments with self to show consistency
    # of results across all replicate comparisons
    cmd = config.get('Pipeline', 'btl_isect_self')
    out_folder = os.path.join(btl_deep_folder, 'isect_raw_self')
    btl_isect_raw_self = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_isect_raw_self',
                                        input=output_from(sciddo_deep_dump_raw_t1),
                                        filter=suffix('.bed'),
                                        output='.tsv',
                                        output_dir=out_folder,
                                        extras=[cmd, jobcall])
    btl_isect_raw_self = btl_isect_raw_self.mkdir(out_folder)

    # intersect raw (unmerged) segments with segments from
    # other segmentations to show robust identification of
    # hotspots of changing chromatin
    cmd = config.get('Pipeline', 'btl_isect_pair')
    btl_isect_raw_other_params = make_btools_all_intersect(sciddo_deep_folder,
                                                           os.path.join(btl_deep_folder, 'isect_raw_other'),
                                                           'raw-t1',
                                                           cmd, jobcall)
    btl_isect_raw_other = pipe.files(sci_obj.get_jobf('inpair_out'),
                                     btl_isect_raw_other_params,
                                     name='btl_isect_raw_other')
    btl_isect_raw_other = btl_isect_raw_other.mkdir(os.path.join(btl_deep_folder, 'isect_raw_other'))
    btl_isect_raw_other = btl_isect_raw_other.follows(sciddo_deep_dump_raw_t1)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Create background set by random sampling
    cmd = config.get('Pipeline', 'sample_random_regions')
    smp_random_genomic_regions = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                                name='smp_random_genomic_regions',
                                                input=output_from(sciddo_deep_dump_hsp_t1),
                                                filter=formatter(match_hsp_bed),
                                                output='{subpath[0][1]}/rand_t1/{SEG[0]}_hg38_{COMP[0]}_{SCORING[0]}_rand-{THRESHOLD[0]}.bed',
                                                extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    collect_random_sets = pipe.originate(lambda x: x,
                                         collect_random_set_files(workdir),
                                         name='collect_random_sets')
    collect_random_sets = collect_random_sets.follows(smp_random_genomic_regions)

    # intersect random regions with various annotation files
    cmd = config.get('Pipeline', 'btl_isect_hsp_any')
    out_folder = os.path.join(btl_deep_folder, 'isect_rand_any')
    match_rand_bed = '(?P<SEG>[a-z0-9]+)_hg38_(?P<COMP>\w+)_(?P<SCORING>[a-z]+)_rand\-(?P<THRESHOLD>t[0-9\.]+)\.bed'
    out_name = 'rand_ovl_reg_{SEG[0]}_{COMP[0]}_{SCORING[0]}.{THRESHOLD[0]}.tsv'
    btl_isect_rand_any = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_isect_rand_any',
                                        input=output_from(collect_random_sets),
                                        filter=formatter(match_rand_bed),
                                        output=os.path.join(out_folder, out_name),
                                        extras=[cmd, jobcall])
    btl_isect_rand_any = btl_isect_rand_any.mkdir(out_folder)

    # === CHANGED ORIENTATION ===
    # INTERSECT FILE X WITH HSPs
    # ===========================

    # intersect genes with HSPs
    cmd = config.get('Pipeline', 'btl_isect_gene_hsp')
    out_folder = os.path.join(btl_deep_folder, 'isect_gene_hsp')
    out_name = 'gene_ovl_hsp_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_isect_gene_hsp = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_isect_gene_hsp',
                                        input=output_from(sciddo_deep_dump_hsp_t1),
                                        filter=formatter(match_hsp_bed),
                                        output=os.path.join(out_folder, out_name),
                                        extras=[cmd, jobcall])
    btl_isect_gene_hsp = btl_isect_gene_hsp.mkdir(out_folder)

    # intersect promoters with HSPs
    cmd = config.get('Pipeline', 'btl_isect_prom_hsp')
    out_folder = os.path.join(btl_deep_folder, 'isect_prom_hsp')
    out_name = 'prom_ovl_hsp_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_isect_prom_hsp = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_isect_prom_hsp',
                                        input=output_from(sciddo_deep_dump_hsp_t1),
                                        filter=formatter(match_hsp_bed),
                                        output=os.path.join(out_folder, out_name),
                                        extras=[cmd, jobcall])
    btl_isect_prom_hsp = btl_isect_prom_hsp.mkdir(out_folder)

    # intersect enhancers with HSPs
    cmd = config.get('Pipeline', 'btl_isect_enh_hsp')
    out_folder = os.path.join(btl_deep_folder, 'isect_enh_hsp')
    out_name = 'enh_ovl_hsp_{SEG[0]}_{COMP[0]}_{SCORING[0]}.tsv'
    btl_isect_enh_hsp = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='btl_isect_enh_hsp',
                                       input=output_from(sciddo_deep_dump_hsp_t1),
                                       filter=formatter(match_hsp_bed),
                                       output=os.path.join(out_folder, out_name),
                                       extras=[cmd, jobcall])
    btl_isect_enh_hsp = btl_isect_enh_hsp.mkdir(out_folder)

    # intersect genes/enhancers with PePr peaks
    cmd = config.get('Pipeline', 'btl_isect_pair')
    out_folder = os.path.join(btl_deep_folder, 'isect_any_pepr')
    enhancer_file = os.path.join(config.get('References', 'project_ref'), 'hg38_genehancer_gencodeV21_complete.bed')
    btl_isect_any_pepr_params = make_pepr_any_intersect_params(deep_to_bed, enhancer_file,
                                                               pepr_deep_output, out_folder,
                                                               cmd, jobcall)
    btl_isect_any_pepr = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    btl_isect_any_pepr_params,
                                    name='btl_isect_any_pepr')
    btl_isect_any_pepr = btl_isect_any_pepr.mkdir(out_folder)


    # cmd = config.get('Pipeline', 'btl_deep_isect_hsp')
    # btl_deep_isect_raw_params = make_btools_all_intersect(sciddo_deep_folder,
    #                                                       os.path.join(btl_deep_folder, 'raw_isect'),
    #                                                       'raw',
    #                                                       cmd, jobcall)
    # btl_deep_isect_raw = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                                 btl_deep_isect_raw_params,
    #                                 name='btl_deep_isect_raw')
    # btl_deep_isect_raw = btl_deep_isect_raw.follows(btl_deep_isect_hsp)
    #
    # # now intersect merged segments with E < 1
    # # with various annotations for later analysis
    #
    # # define regexp for dumped HSP BED files
    # # deep_hg38_cmm18_HG_vs_He_hsp-emission.bed
    # match_hsp_bed = 'deep_hg38_(?P<SEG>[a-z0-9]+)_(?P<COMP>\w+)_hsp\-(?P<SCORE>[a-z]+)\.bed'
    #
    #
    #
    #
    #
    # cmd = config.get('Pipeline', 'btl_deep_isect_enh')
    # out_folder = os.path.join(btl_deep_folder, 'enh_isect')
    # out_name = 'deep_hsp_ovl_enh_{SEG[0]}_{COMP[0]}_{SCORE[0]}.tsv'
    # btl_deep_isect_enh = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                                     name='btl_deep_isect_enh',
    #                                     input=output_from(sciddo_deep_dump_hsp_t1),
    #                                     filter=formatter(match_hsp_bed),
    #                                     output=os.path.join(out_folder, out_name),
    #                                     extras=[cmd, jobcall])
    # btl_deep_isect_enh = btl_deep_isect_enh.mkdir(out_folder)
    #
    # # intersect genes and enhancers with HSPs
    # # to better understand how HSP hits alter
    # # gene expression patterns
    # cmd = config.get('Pipeline', 'btl_deep_isect_inv')
    # out_folder = os.path.join(btl_deep_folder, 'inv_isect')
    # enhancer_file = os.path.join(config.get('References', 'project_ref'), 'hg38_genehancer_gencodeV21_complete.bed')
    # btl_deep_isect_inv_params = make_invert_intersect_params(deep_to_bed, enhancer_file,
    #                                                          sciddo_deep_folder, out_folder,
    #                                                          cmd, jobcall)
    # btl_deep_isect_inv = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                                 btl_deep_isect_inv_params,
    #                                 name='btl_deep_isect_inv')
    # btl_deep_isect_inv = btl_deep_isect_inv.mkdir(out_folder)
    #
    #
    # outfolder = os.path.join(btl_deep_folder, 'pepr_locus_isect')
    # cmd = config.get('Pipeline', 'btl_deep_isect_hsp_deg_pepr')
    # btl_deep_isect_hsp_deg_pepr_params = make_btools_hsp_pepr_isect(sciddo_deep_folder,
    #                                                                 pepr_deep_output,
    #                                                                 deep_to_bed,
    #                                                                 'locus', 'diff',
    #                                                                 outfolder, cmd, jobcall)
    #
    # btl_deep_isect_hsp_deg_pepr = pipe.files(sci_obj.get_jobf('in_out'),
    #                                          btl_deep_isect_hsp_deg_pepr_params,
    #                                          name='btl_deep_isect_hsp_deg_pepr')
    # btl_deep_isect_hsp_deg_pepr = btl_deep_isect_hsp_deg_pepr.mkdir(outfolder)
    # btl_deep_isect_hsp_deg_pepr = btl_deep_isect_hsp_deg_pepr.follows(conv_deep_diff)
    # btl_deep_isect_hsp_deg_pepr = btl_deep_isect_hsp_deg_pepr.follows(sciddo_deep_dump_hsp_t1)
    # btl_deep_isect_hsp_deg_pepr = btl_deep_isect_hsp_deg_pepr.follows(pepr_deep_filter)
    #
    # outfolder = os.path.join(btl_deep_folder, 'pepr_locus_isect')
    # cmd = config.get('Pipeline', 'btl_deep_isect_hsp_deg_pepr')
    # btl_deep_isect_hsp_stg_pepr_params = make_btools_hsp_pepr_isect(sciddo_deep_folder,
    #                                                                 pepr_deep_output,
    #                                                                 deep_to_bed,
    #                                                                 'locus', 'stable',
    #                                                                 outfolder, cmd, jobcall)
    #
    # btl_deep_isect_hsp_stg_pepr = pipe.files(sci_obj.get_jobf('in_out'),
    #                                          btl_deep_isect_hsp_stg_pepr_params,
    #                                          name='btl_deep_isect_hsp_stg_pepr')
    # btl_deep_isect_hsp_stg_pepr = btl_deep_isect_hsp_stg_pepr.mkdir(outfolder)
    # btl_deep_isect_hsp_stg_pepr = btl_deep_isect_hsp_stg_pepr.follows(conv_deep_diff)
    # btl_deep_isect_hsp_stg_pepr = btl_deep_isect_hsp_stg_pepr.follows(sciddo_deep_dump_hsp_t1)
    # btl_deep_isect_hsp_stg_pepr = btl_deep_isect_hsp_stg_pepr.follows(pepr_deep_filter)
    # btl_deep_isect_hsp_stg_pepr = btl_deep_isect_hsp_stg_pepr.follows(btl_deep_isect_hsp_deg_pepr)
    #
    # ======================================
    # SPECIAL TASKS FOR INDIVIDUAL USECASES
    # ======================================

    # use case 1: HepG2 / EP300 enhancer activity
    # cmm18_hg38_HG_vs_Ma_penem_hsp-t1.bed

    usecase_out = os.path.join(workdir, 'usecases')

    uc1_folder = os.path.join(usecase_out, 'uc1_p300')

    cmd = config.get('Pipeline', 'btl_uc1_isect_p300')
    btl_uc1_isect_p300 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_uc1_isect_p300',
                                        input=output_from(sciddo_deep_dump_hsp_t1),
                                        filter=formatter('.*/(?P<SEG>cmm18)_hg38_(?P<COMP>HG_vs_Mo)_(?P<SCORE>penem)_hsp\-t1\.bed'),
                                        output=os.path.join(uc1_folder, 'uc1_{SEG[0]}_{SCORE[0]}_{COMP[0]}_p300_hsp.bed'),
                                        extras=[cmd, jobcall])
    btl_uc1_isect_p300 = btl_uc1_isect_p300.mkdir(uc1_folder)

    cmd = config.get('Pipeline', 'btl_uc1_p300_uniq')
    btl_uc1_p300_uniq = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='btl_uc1_p300_uniq',
                                       input=output_from(btl_uc1_isect_p300),
                                       filter=suffix('_hsp.bed'),
                                       output='_uniq.bed',
                                       extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'sciddo_uc1_enh_onoff')
    sciddo_uc1_enh_onoff = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='sciddo_uc1_enh_onoff',
                                          input=output_from(sciddo_deep_scan),
                                          filter=formatter('sciddo-run_hg38_(?P<SEG>cmm18)_(?P<SCORE>penem)_(?P<C1>HG)_vs_(?P<C2>Mo)\.h5'),
                                          output=os.path.join(uc1_folder, 'uc1_{SEG[0]}_{SCORE[0]}_{C1[0]}_vs_{C2[0]}_enh-on-off_splits.bed'),
                                          extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_uc1_enh_offon')
    sciddo_uc1_enh_offon = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='sciddo_uc1_enh_offon',
                                          input=output_from(sciddo_deep_scan),
                                          filter=formatter('sciddo-run_hg38_(?P<SEG>cmm18)_(?P<SCORE>penem)_(?P<C1>HG)_vs_(?P<C2>Mo)\.h5'),
                                          output=os.path.join(uc1_folder, 'uc1_{SEG[0]}_{SCORE[0]}_{C1[0]}_vs_{C2[0]}_enh-off-on_splits.bed'),
                                          extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'btl_uc1_isect_switch')
    btl_uc1_isect_switch = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='btl_uc1_isect_switch',
                                          input=output_from(sciddo_uc1_enh_onoff, sciddo_uc1_enh_offon),
                                          filter=formatter('(?P<BASENAME>\w+)_(?P<DIR>enh[\-onf]+)_splits\.bed'),
                                          output=os.path.join(uc1_folder, '{BASENAME[0]}_{DIR[0]}_p300.bed'),
                                          extras=[cmd, jobcall])
    btl_uc1_isect_switch = btl_uc1_isect_switch.follows(btl_uc1_p300_uniq)

    # run all task
    statediff_run_all = pipe.merge(task_func=touch_checkfile,
                                   name='statediff_run_all',
                                   input=output_from(sciddo_deep_conv_cmm18, sciddo_deep_conv_ecs,
                                                     sciddo_deep_stats, sciddo_deep_score,
                                                     sciddo_deep_scan, sciddo_deep_dump_seg,
                                                     sciddo_deep_dump_hsp_t1, sciddo_deep_dump_hsp_t100,
                                                     sciddo_deep_dump_raw_t1, smp_random_genomic_regions,
                                                     sciddo_deep_dump_trans,
                                                     btl_isect_rand_any,
                                                     btl_isect_any_pepr, btl_uc1_p300_uniq, btl_uc1_isect_switch,
                                                     btl_isect_enh_hsp, btl_isect_gene_hsp, btl_isect_hsp_rgb,
                                                     btl_isect_raw_other, btl_deep_isect_t100, btl_isect_raw_self,
                                                     btl_isect_hsp_any, btl_isect_hsp_gene),
                                   output=os.path.join(workdir, 'statediff_run_all.chk'))
    return pipe
