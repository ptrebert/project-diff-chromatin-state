#!/usr/bin/env Rscript

suppressMessages(library(getopt))
suppressMessages(library(tximport))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))


make_sample_annotation = function(quant_path)
{
    sample_folders = list.dirs(quant_path, recursive=FALSE, full.names=FALSE)
    num_rows = length(sample_folders)

    sample_ids = c()
    karyotypes = c()
    organs = c()
    celltypes = c()
    for (sid in sample_folders)
    {
        parts = unlist(strsplit(sid, '_'))
        sample_ids = append(sample_ids, sid)
        organs = append(organs, substr(parts[3], 1, 2))
        celltypes = append(celltypes, substr(parts[3], 3, 4))
        karyotypes = append(karyotypes, substr(parts[2], 2, 2))
    }

    sample_annotation = data.frame(karyotype=karyotypes,
                                   organ=organs,
                                   celltype=celltypes,
                                   stringsAsFactors=TRUE)
    rownames(sample_annotation) = sample_ids
    return (sample_annotation)
}


collect_quant_files = function(quant_path)
{
    sample_folders = list.dirs(quant_path, recursive=FALSE, full.names=FALSE)
    paths = paste(quant_path, sample_folders, 'quant.sf', sep='/')
    names(paths) = sample_folders
    return (paths)

}

cmd_option_spec = matrix(c(
                        'loadpath', 'p', 1, 'character',
                        'genemap', 'g', 1, 'character',
                        'outpath', 'o', 1, 'character'
                        ), byrow=TRUE, ncol=4)

cmd_opts = getopt(cmd_option_spec)

if (is.null(cmd_opts$loadpath))
{
    stop('Need to specify full path to load Salmon quantification files')
}

genemap = read.table(cmd_opts$genemap, sep='\t', header=FALSE, row.names=NULL, col.names=c('TXNAME', 'GENEID'))

smp_ann = make_sample_annotation(cmd_opts$loadpath)
quant_files = collect_quant_files(cmd_opts$loadpath)

txi = tximport(quant_files, type="salmon", tx2gene=genemap)

dds = DESeqDataSetFromTximport(txi, smp_ann, ~ celltype)

dds = DESeq(dds)

# dump results
res_hg_vs_he = results(dds, contrast=c('celltype', 'HG', 'He'))
write.table(as.data.frame(res_hg_vs_he),
            paste(cmd_opts$outpath, 'deseq2_HG_vs_He.tsv', sep='/'),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

res_he_vs_mo = results(dds, contrast=c('celltype', 'He', 'Mo'))
write.table(as.data.frame(res_he_vs_mo),
            paste(cmd_opts$outpath, 'deseq2_He_vs_Mo.tsv', sep='/'),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

res_he_vs_ma = results(dds, contrast=c('celltype', 'He', 'Ma'))
write.table(as.data.frame(res_he_vs_ma),
            paste(cmd_opts$outpath, 'deseq2_He_vs_Ma.tsv', sep='/'),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

res_hg_vs_mo = results(dds, contrast=c('celltype', 'HG', 'Mo'))
write.table(as.data.frame(res_hg_vs_mo),
            paste(cmd_opts$outpath, 'deseq2_HG_vs_Mo.tsv', sep='/'),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

res_hg_vs_ma = results(dds, contrast=c('celltype', 'HG', 'Ma'))
write.table(as.data.frame(res_hg_vs_ma),
            paste(cmd_opts$outpath, 'deseq2_HG_vs_Ma.tsv', sep='/'),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

res_mo_vs_ma = results(dds, contrast=c('celltype', 'Ma', 'Mo'))
write.table(as.data.frame(res_mo_vs_ma),
            paste(cmd_opts$outpath, 'deseq2_Ma_vs_Mo.tsv', sep='/'),
            col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

quit(save="no", status=0)
