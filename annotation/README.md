# Reference data used for project

## Gene expression
GENCODE v21 main annotation file (gencode.v21.annotation.gtf.gz) and protein coding transcripts.

Salmon index built as follows:

```
salmon index -i idx_hg38_GENCv21_k31
 -t gencode.v21.pc_transcripts.fa.gz -k 31
 --gencode --threads 4 --type quasi --perfectHash
```

## Regulatory build

Ensembl Regulatory Build

Release 78 (corresponds to GENCODE v21 - MD5: 6f9c703caa5d23e6d8215d3368f4fa29)

ftp://ftp.ensembl.org/pub/release-78/regulation/homo_sapiens/RegulatoryFeatures_MultiCell.gff.gz

Release 65 (corresponds to GENCODE v10 - MD5: e0e9eb48d1ffef05d84b730e015aa94a)
used for ROADMAP Epigenomics data:

ftp://ftp.ensembl.org/pub/release-65/regulation/homo_sapiens/RegulatoryFeatures_MultiCell.gff.gz

## ENCODE data

#### Experiment ENCSR000EDV
- Snyder lab, EP300 in HepG2, hg38
- peak set: optimal idr thresholded peaks 
- URL: https://www.encodeproject.org/files/ENCFF806JJS/@@download/ENCFF806JJS.bed.gz

#### Experiment ENCSR000BLW
- Myers lab, EP300 in HepG2, hg38
- peak set: optimal idr thresholded peaks
- URL: https://www.encodeproject.org/files/ENCFF674QCU/@@download/ENCFF674QCU.bed.gz

#### Merging of peak files

Print to single file:

```
gunzip -c ENCFF* | sort -V -k1,3 > ENCODE_HepG2_EP300_peaks_ovl.bed
```

Merge overlapping peaks:

```
bedtools merge -i ENCODE_HepG2_EP300_peaks_ovl.bed -c 5,7,8,9,10 -prec 4 -o max,mean,max,mean,mean > ENCODE_HepG2_EP300_peaks.bed
```

Restrict to relevant chromosomes:

```
egrep "chr[0-9X]+\s" ENCODE_HepG2_EP300_peaks.bed > ENCODE_HepG2_EP300_peaks_1-22X.bed
```

## Misc annotation

Not important for project - used to reproduce some results from papers
URL for downloaded PAM120 scoring matrix:

http://www.quretec.com/u/vilo/edu/2002-03/Tekstialgoritmid_I/Loengud/Loeng3_Edit_Distance/bcorum_copy/seq_align5.htm