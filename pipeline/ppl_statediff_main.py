# coding=utf-8

import os as os
import re as re
import csv as csv
import operator as op
import collections as col
import itertools as itt

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
    header = ['chrom', 'start', 'end', 'name', 'score',
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
    with open(outputfile, 'w', newline='') as filtered:
        _ = filtered.write('#')
        writer = csv.DictWriter(filtered, delimiter='\t', fieldnames=header)
        writer.writeheader()
        writer.writerows(recs)
    statfile = outputfile.replace('sig.bed', 'stat.txt')
    with open(statfile, 'w') as stats:
        for k, v in zip(['total_num', 'shortest_peak', 'longest_peak',
                         'percentile_10', 'percentile_25'],
                        [total_peaks, shortest, longest, pct10_len, pct25_len]):
            _ = stats.write(k + '\t' + str(v) + '\n')
    return outputfile


def make_sciddo_scan_params(inputfolder, sub_folder, groups, cmd, jobcall):
    """
    :param inputfolder:
    :param sub_folder:
    :param groups:
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
                if fn.endswith('.h5') and fn.startswith('deep_sciddo_'):
                    in_path = os.path.join(root, fn)
                    out_dir = os.path.join(root, sub_folder)
                    os.makedirs(out_dir, exist_ok=True)
                    outfile = fn.replace('_sciddo_', '_hsp_')
                    for g1, g2 in itt.combinations(groups, 2):
                        tmp = cmd.format(**{'grp1': g1, 'grp2': g2})
                        out_tmp = outfile.replace('.h5', '_' + g1 + '_vs_' + g2 + '.h5')
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
    # deep_hg38_cmm18_HG_vs_He_hsp-emission.bed
    match_info = re.compile("deep_hg38_(?P<SEGMENT>[a-z0-9]+)_(?P<COMP>\w+)_(?P<DATA>(hsp|raw))-(?P<SCORING>[a-z]+)\.bed")
    inputfiles = []
    for root, dirs, datafiles in os.walk(inputroot):
        if root.endswith('data_dump') and datafiles:
            for bf in datafiles:
                if bf.endswith('.bed') and ('_hsp-' in bf or '_raw-' in bf):
                    inputfiles.append(os.path.join(root, bf))
    args = []
    for b1, b2 in itt.combinations_with_replacement(sorted(inputfiles), 2):
        bn1, bn2 = os.path.basename(b1), os.path.basename(b2)
        b1_info = match_info.match(bn1)
        assert b1_info is not None, 'Could not match {}'.format(bn1)
        if b1_info.group('DATA') != datatype:
            continue
        b2_info = match_info.match(bn2)
        if b2_info.group('DATA') != datatype:
            continue
        outname = '{}-{}-{}-ovl-{}-{}-{}.tsv'.format(b1_info.group('SEGMENT'),
                                                     b1_info.group('COMP'),
                                                     b1_info.group('SCORING'),
                                                     b2_info.group('SEGMENT'),
                                                     b2_info.group('COMP'),
                                                     b2_info.group('SCORING'))

        a = int(b1_info.group('SEGMENT') == b2_info.group('SEGMENT'))
        b = int(b1_info.group('COMP') == b2_info.group('COMP'))
        c = int(b1_info.group('SCORING') == b2_info.group('SCORING'))
        subfolder = subfolders[(a, b, c)]
        outpath = os.path.join(outroot, subfolder, outname)
        args.append([[b1, b2], outpath, cmd, jobcall])
    return args


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
    # Differential peak calling with replicates
    # Set folder for DEEP samples
    pepr_deep_folder = os.path.join(workdir, 'pepr', 'deep')
    pepr_deep_output = os.path.join(pepr_deep_folder, 'diff_out')

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('Py2EnvConfig')))
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

    # filter "passed" BED files for q-value
    # Ma_vs_Mo_H3K4me1__PePr_post-peak1.chk
    pepr_deep_filter = pipe.transform(task_func=filter_sort_pepr_peaks,
                                      name='pepr_deep_filter',
                                      input=output_from(pepr_deep_post),
                                      filter=formatter('(?P<SAMPLEID>\w+)__PePr_post\-peak(?P<REP>(1|2))\.chk'),
                                      output=os.path.join('{path[0]}', '{SAMPLEID[0]}__PePr_chip{REP[0]}_peak.sig.bed'))
    pepr_deep_filter = pepr_deep_filter.jobs_limit(4)

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
                                                                'deep_sciddo_hg38_cmm18.conv-chk'),
                                            extras=[cmd, jobcall])
    sciddo_deep_conv_cmm18 = sciddo_deep_conv_cmm18.mkdir(sciddo_deep_folder)

    # convert EpiCSeg
    cmd = config.get('Pipeline', 'sciddo_deep_conv_ecs')
    sciddo_deep_conv_ecs = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='sciddo_deep_conv_ecs',
                                          input=output_from(ecs_deep_seg),
                                          filter=formatter('ecs_deep_seg_(?P<STNUM>[0-9]+)\.chk'),
                                          output=os.path.join(sciddo_deep_folder, 'ecs{STNUM[0]}',
                                                              'deep_sciddo_hg38_ecs{STNUM[0]}.conv-chk'),
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

    # ==============================
    # Scan for differential regions
    # based on cell-type grouping
    cmd = config.get('Pipeline', 'sciddo_deep_scan')
    groups = config.get('Pipeline', 'sciddo_deep_groups').split()
    sciddo_deep_scan_pub_params = make_sciddo_scan_params(sciddo_deep_folder, 'hsp_run',
                                                          groups, cmd, jobcall)
    sciddo_deep_scan = pipe.files(sci_obj.get_jobf('in_out'),
                                  sciddo_deep_scan_pub_params,
                                  name='sciddo_deep_scan')
    sciddo_deep_scan = sciddo_deep_scan.follows(sciddo_deep_score)

    # ============================================
    # Generate HSP samples by randomized sampling
    # THIS TAKES AN INSANE AMOUNT OF TIME!
    cmd = config.get('Pipeline', 'sciddo_deep_smprand')
    sciddo_deep_smprand = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                         name='sciddo_deep_smprand',
                                         input=output_from(sciddo_deep_score),
                                         filter=formatter('scorechk$'),
                                         output=os.path.join('{path[0]}', 'smp_run', '{basename[0]}_smprand.h5'),
                                         extras=[cmd, jobcall])
    sciddo_deep_smprand = sciddo_deep_smprand.active_if(config.getboolean('Pipeline', 'weekend'))

    # ================================================================================
    # Dump data for visualization
    # and HSPs

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # dump fully annotated state segmentations
    deep_dset_re = 'deep_sciddo_hg38_(?P<SEGMENT>[cmes180]{5})\.conv-chk'
    cmd = config.get('Pipeline', 'sciddo_deep_dump_seg')
    sciddo_deep_dump_seg = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                                          name='sciddo_deep_dump_seg',
                                          input=output_from(sciddo_deep_conv_cmm18, sciddo_deep_conv_ecs),
                                          filter=formatter(deep_dset_re),
                                          output='{path[0]}/data_dump/{SEGMENT[0]}_*.bed',
                                          extras=['{path[0]}/data_dump', '{SEGMENT[0]}_*.bed', cmd, jobcall])

    deep_scan_re = 'deep_hsp_hg38_(?P<SEGMENT>[cmes180]{5})_CELLTYPE_(?P<C1>[HGeMoa]{2})_vs_CELLTYPE_(?P<C2>[HGeMoa]{2})\.h5'

    cmd = config.get('Pipeline', 'sciddo_deep_dump_emission')
    sciddo_deep_dump_emission = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                               name='sciddo_deep_dump_emission',
                                               input=output_from(sciddo_deep_scan),
                                               filter=formatter(deep_scan_re),
                                               output='{subpath[0][1]}/data_dump/deep_hg38_{SEGMENT[0]}_{C1[0]}_vs_{C2[0]}_emission-scores.bg',
                                               extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_replicate')
    sciddo_deep_dump_replicate = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                                name='sciddo_deep_dump_replicate',
                                                input=output_from(sciddo_deep_scan),
                                                filter=formatter(deep_scan_re),
                                                output='{subpath[0][1]}/data_dump/deep_hg38_{SEGMENT[0]}_{C1[0]}_vs_{C2[0]}_replicate-scores.bg',
                                                extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_hsp_ems')
    sciddo_deep_dump_hsp_ems = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='sciddo_deep_dump_hsp_ems',
                                              input=output_from(sciddo_deep_scan),
                                              filter=formatter(deep_scan_re),
                                              output='{subpath[0][1]}/data_dump/deep_hg38_{SEGMENT[0]}_{C1[0]}_vs_{C2[0]}_hsp-emission.bed',
                                              extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_hsp_rep')
    sciddo_deep_dump_hsp_rep = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='sciddo_deep_dump_hsp_rep',
                                              input=output_from(sciddo_deep_scan),
                                              filter=formatter(deep_scan_re),
                                              output='{subpath[0][1]}/data_dump/deep_hg38_{SEGMENT[0]}_{C1[0]}_vs_{C2[0]}_hsp-replicate.bed',
                                              extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_raw_ems')
    sciddo_deep_dump_raw_ems = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='sciddo_deep_dump_raw_ems',
                                              input=output_from(sciddo_deep_scan),
                                              filter=formatter(deep_scan_re),
                                              output='{subpath[0][1]}/data_dump/deep_hg38_{SEGMENT[0]}_{C1[0]}_vs_{C2[0]}_raw-emission.bed',
                                              extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_deep_dump_raw_rep')
    sciddo_deep_dump_raw_rep = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='sciddo_deep_dump_raw_rep',
                                              input=output_from(sciddo_deep_scan),
                                              filter=formatter(deep_scan_re),
                                              output='{subpath[0][1]}/data_dump/deep_hg38_{SEGMENT[0]}_{C1[0]}_vs_{C2[0]}_raw-replicate.bed',
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
                                     input=output_from(sciddo_deep_dump_emission, sciddo_deep_dump_replicate),
                                     filter=suffix('.bg'),
                                     output='.bw',
                                     extras=[cmd, jobcall])

    # intersect all vs all HSP files for various comparison plots
    # this operates on the level of merged HSPs
    btl_deep_folder = os.path.join(workdir, 'bedtools', 'deep')

    cmd = config.get('Pipeline', 'btl_deep_isect_hsp')
    btl_deep_isect_hsp_params = make_btools_all_intersect(sciddo_deep_folder,
                                                          os.path.join(btl_deep_folder, 'hsp_isect'),
                                                          'hsp',
                                                          cmd, jobcall)
    btl_deep_isect_hsp = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    btl_deep_isect_hsp_params,
                                    name='btl_deep_isect_hsp')
    btl_deep_isect_hsp = btl_deep_isect_hsp.mkdir(os.path.join(btl_deep_folder, 'hsp_isect'))
    btl_deep_isect_hsp = btl_deep_isect_hsp.follows(sciddo_deep_dump_hsp_ems)
    btl_deep_isect_hsp = btl_deep_isect_hsp.follows(sciddo_deep_dump_hsp_rep)

    cmd = config.get('Pipeline', 'btl_deep_isect_hsp')
    btl_deep_isect_raw_params = make_btools_all_intersect(sciddo_deep_folder,
                                                          os.path.join(btl_deep_folder, 'raw_isect'),
                                                          'raw',
                                                          cmd, jobcall)
    btl_deep_isect_raw = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    btl_deep_isect_raw_params,
                                    name='btl_deep_isect_raw')
    btl_deep_isect_raw = btl_deep_isect_raw.follows(btl_deep_isect_hsp)
    btl_deep_isect_raw = btl_deep_isect_raw.follows(sciddo_deep_dump_raw_ems)
    btl_deep_isect_raw = btl_deep_isect_raw.follows(sciddo_deep_dump_raw_rep)

    # now intersect merged segments with E < 1
    # with various annotations for later analysis

    # define regexp for dumped HSP BED files
    # deep_hg38_cmm18_HG_vs_He_hsp-emission.bed
    match_hsp_bed = 'deep_hg38_(?P<SEG>[a-z0-9]+)_(?P<COMP>\w+)_hsp\-(?P<SCORE>[a-z]+)\.bed'

    cmd = config.get('Pipeline', 'btl_deep_isect_genes')
    out_folder = os.path.join(btl_deep_folder, 'gene_isect')
    out_name = 'deep_hsp_ovl_degenes_{SEG[0]}_{COMP[0]}_{SCORE[0]}.tsv'
    btl_deep_isect_genes = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='btl_deep_isect_genes',
                                          input=output_from(sciddo_deep_dump_hsp_ems,
                                                            sciddo_deep_dump_hsp_rep),
                                          filter=formatter(match_hsp_bed),
                                          output=os.path.join(out_folder, out_name),
                                          extras=[cmd, jobcall])
    btl_deep_isect_genes = btl_deep_isect_genes.mkdir(out_folder)

    cmd = config.get('Pipeline', 'btl_deep_isect_prom')
    out_folder = os.path.join(btl_deep_folder, 'prom_isect')
    out_name = 'deep_hsp_ovl_degprom_{SEG[0]}_{COMP[0]}_{SCORE[0]}.tsv'
    btl_deep_isect_prom = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                         name='btl_deep_isect_prom',
                                         input=output_from(sciddo_deep_dump_hsp_ems,
                                                           sciddo_deep_dump_hsp_rep),
                                         filter=formatter(match_hsp_bed),
                                         output=os.path.join(out_folder, out_name),
                                         extras=[cmd, jobcall])
    btl_deep_isect_prom = btl_deep_isect_prom.mkdir(out_folder)

    cmd = config.get('Pipeline', 'btl_deep_isect_rgb')
    out_folder = os.path.join(btl_deep_folder, 'ensreg_isect')
    out_name = 'deep_hsp_ovl_ensreg_{SEG[0]}_{COMP[0]}_{SCORE[0]}.tsv'
    btl_deep_isect_ensreg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                           name='btl_deep_isect_ensreg',
                                           input=output_from(sciddo_deep_dump_hsp_ems,
                                                             sciddo_deep_dump_hsp_rep),
                                           filter=formatter(match_hsp_bed),
                                           output=os.path.join(out_folder, out_name),
                                           extras=[cmd, jobcall])
    btl_deep_isect_ensreg = btl_deep_isect_ensreg.mkdir(out_folder)

    cmd = config.get('Pipeline', 'btl_deep_isect_enh')
    out_folder = os.path.join(btl_deep_folder, 'enh_isect')
    out_name = 'deep_hsp_ovl_enh_{SEG[0]}_{COMP[0]}_{SCORE[0]}.tsv'
    btl_deep_isect_enh = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_deep_isect_enh',
                                        input=output_from(sciddo_deep_dump_hsp_ems,
                                                          sciddo_deep_dump_hsp_rep),
                                        filter=formatter(match_hsp_bed),
                                        output=os.path.join(out_folder, out_name),
                                        extras=[cmd, jobcall])
    btl_deep_isect_enh = btl_deep_isect_enh.mkdir(out_folder)

    # ======================================
    # SPECIAL TASKS FOR INDIVIDUAL USECASES
    # ======================================

    # use case 1: HepG2 / EP300 enhancer activity

    usecase_out = os.path.join(workdir, 'usecases')

    uc1_folder = os.path.join(usecase_out, 'uc1_p300')

    cmd = config.get('Pipeline', 'btl_uc1_isect_p300')
    btl_uc1_isect_p300 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                        name='btl_uc1_isect_p300',
                                        input=output_from(sciddo_deep_dump_hsp_ems),
                                        filter=formatter('.*/deep_hg38_(?P<SEG>cmm18)_(?P<COMP>HG_vs_Mo)_hsp\-emission\.bed'),
                                        output=os.path.join(uc1_folder, 'uc1_{SEG[0]}_{COMP[0]}_p300_hsp.bed'),
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
                                          filter=formatter('deep_hsp_hg38_(?P<SEG>cmm18)_CELLTYPE_(?P<C1>HG)_vs_CELLTYPE_(?P<C2>Mo)\.h5'),
                                          output=os.path.join(uc1_folder, 'uc1_{SEG[0]}_{C1[0]}_vs_{C2[0]}_enh-on-off_splits.bed'),
                                          extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'sciddo_uc1_enh_offon')
    sciddo_uc1_enh_offon = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                          name='sciddo_uc1_enh_offon',
                                          input=output_from(sciddo_deep_scan),
                                          filter=formatter('deep_hsp_hg38_(?P<SEG>cmm18)_CELLTYPE_(?P<C1>HG)_vs_CELLTYPE_(?P<C2>Mo)\.h5'),
                                          output=os.path.join(uc1_folder, 'uc1_{SEG[0]}_{C1[0]}_vs_{C2[0]}_enh-off-on_splits.bed'),
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





    # cmd = config.get('Pipeline', 'btl_deep_isect_pub')
    # btl_deep_isect_pub = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                                     name='btl_deep_isect_pub',
    #                                     input=output_from(sciddo_deep_dump_pub),
    #                                     filter=suffix('.bed'),
    #                                     output='.ovl.bed',
    #                                     extras=[cmd, jobcall])
    #
    # cmd = config.get('Pipeline', 'btl_deep_isect_deb_pub')
    # hsp_isect_out = os.path.join(sciddo_deep_folder, 'hsp_isect')
    # filter_re = 'hsp_(?P<C1>[A-Za-z]{2})_vs_(?P<C2>[A-Za-z]{2})_emission.bed'
    # btl_deep_isect_deb_pub = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                                         name='btl_deep_isect_deb_pub',
    #                                         input=output_from(sciddo_deep_dump_pub),
    #                                         filter=formatter(filter_re),
    #                                         output=os.path.join(hsp_isect_out,
    #                                                             'hsp_{C1[0]}_vs_{C2[0]}_emission_de-body.bed'),
    #                                         extras=[cmd, jobcall])
    # btl_deep_isect_deb_pub = btl_deep_isect_deb_pub.mkdir(hsp_isect_out)
    # btl_deep_isect_deb_pub = btl_deep_isect_deb_pub.follows(conv_deep_diff)
    #
    # cmd = config.get('Pipeline', 'btl_deep_isect_rgb_pub')
    # hsp_isect_out = os.path.join(sciddo_deep_folder, 'hsp_isect')
    # filter_re = 'hsp_(?P<C1>[A-Za-z]{2})_vs_(?P<C2>[A-Za-z]{2})_emission.bed'
    # btl_deep_isect_rgb_pub = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                                         name='btl_deep_isect_rgb_pub',
    #                                         input=output_from(sciddo_deep_dump_pub),
    #                                         filter=formatter(filter_re),
    #                                         output=os.path.join(hsp_isect_out,
    #                                                             'hsp_{C1[0]}_vs_{C2[0]}_emission_ensregb.bed'),
    #                                         extras=[cmd, jobcall])
    # btl_deep_isect_rgb_pub = btl_deep_isect_rgb_pub.mkdir(hsp_isect_out)


    ### below: full DEEP dataset, contains unpublished data, ignore

    # cmd = config.get('Pipeline', 'sciddo_deep_conv_prv')
    # sciddo_deep_conv_prv = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                                    name='sciddo_deep_conv_prv',
    #                                    input=output_from(cmm_deep_mkseg),
    #                                    filter=formatter('chk$'),
    #                                    output=os.path.join(sciddo_deep_folder, 'sciddo_deep_conv_prv.chk'),
    #                                    extras=[cmd, jobcall])
    # sciddo_deep_conv_prv = sciddo_deep_conv_prv.mkdir(sciddo_deep_folder)
    # sciddo_deep_conv_prv = sciddo_deep_conv_prv.active_if(False)

    # =====================================================================
    # Process state segmentation of ROADMAP samples
    # Set folders for ROADMAP samples

    # ======================================================================
    # ROADMAP dataset does not contain appropriate replicates, thus it is
    # less interesting to work with - ignore for now
    #
    # sciddo_remc_folder = os.path.join(workdir, 'sciddo', 'remc')
    #
    # sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('EnvConfig')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # cmd = config.get('Pipeline', 'sciddo_remc_conv')
    # sciddo_remc_conv = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                                name='sciddo_remc_conv',
    #                                input=output_from(sciddo_deep_conv_pub),
    #                                filter=formatter('chk$'),
    #                                output=os.path.join(sciddo_remc_folder, 'sciddo_remc_conv.chk'),
    #                                extras=[cmd, jobcall])
    # sciddo_remc_conv = sciddo_remc_conv.mkdir(sciddo_remc_folder)
    # sciddo_remc_conv = sciddo_remc_conv.active_if(os.path.isfile(config.get('References', 'remc_design')))

    return pipe
