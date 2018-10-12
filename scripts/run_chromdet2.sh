#!/bin/bash

# Submit job 2018-03-02 20:00
# qsub -l h_vmem=80G,mem_free=40G,slots_free=2,tmp_free=20G -j y -N CDPUB2 -S /bin/bash run_chromdet2.sh

export PATH=/TL/epigenetics2/work/pebert/conda/envs/statediff/bin:/home/pebert/work/code/github/ChromDet/ChromDet-master/scripts:/bin:/usr/bin

SEGMENTDIR="/TL/deep/fhgfs/projects/pebert/thesis/projects/statediff/chromdet/deep_pub/pair_segment"

WORKDIR="/TL/deep/fhgfs/projects/pebert/thesis/projects/statediff/chromdet/deep_pub"

OUTDIR="/TL/deep/fhgfs/projects/pebert/thesis/projects/statediff/chromdet/deep_pub/chromdet_out"

LABEL_COLLAPSE="/home/pebert/work/code/mpggit/statediff/annotation/chromdet/chromhmm_18_state_collapse.tsv"
SAMPLE_LABEL="/home/pebert/work/code/mpggit/statediff/annotation/chromdet/deep_sample_labels2.tsv"

cd /home/pebert/work/code/github/ChromDet/ChromDet-master/scripts

date > ${WORKDIR}/deep_pub2.log
echo "" >> ${WORKDIR}/deep_pub2.log

run_S3det_analysis.pl -d ${SEGMENTDIR} -c ${LABEL_COLLAPSE} -h F \
                    -a ${SAMPLE_LABEL} -s /home/pebert/work/code/github/ChromDet/ChromDet-master/S3Det_modified \
                    -o cmm_states -v &>> ${WORKDIR}/deep_pub2.log

echo "" >> ${WORKDIR}/deep_pub2.log
date >> ${WORKDIR}/deep_pub2.log

mail -s "ChromDET finished" pebert@mpi-inf.mpg.de < ${WORKDIR}/deep_pub2.log

exit 0
