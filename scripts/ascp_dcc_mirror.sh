#!/bin/bash

# ascp command line documentation
# http://download.asperasoft.com/download/docs/ascp/3.0/html/index.html

SRCURL=dkfzaspera.dkfz-heidelberg.de

# Removed from list since HepG2 data were renamed
# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/01_HepG2_LiHG_Ct1/replicate1/paired/run140827_SN751_0199_AC52RUACXX/sequence
# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/01_HepG2_LiHG_Ct2/replicate1/paired/run140827_SN751_0199_AC52RUACXX/sequence

# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/43_Hm01_BlMo_Ct/replicate1/paired/run131016_SN471_0153_A_D263EACXX/sequence
# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/43_Hm03_BlMa_Ct/replicate1/paired/run140612_SN303_0188_A_C42GAACXX/sequence
# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/43_Hm03_BlMo_Ct/replicate1/paired/run140612_SN303_0188_A_C42GAACXX/sequence
# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/43_Hm05_BlMa_Ct/replicate1/paired/run140904_SN540_0225_B_C47EEACXX/sequence
# /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/43_Hm05_BlMo_Ct/replicate1/paired/run140904_SN540_0225_B_C47EEACXX/sequence

SRCFOLDERS=(/download/sequencing/strand_specific_mrna_sequencing/view-by-pid/41_Hf02_LiHe_Ct/replicate1/paired/run140212_SN758_0152_BC3ETGACXX/sequence
            /download/sequencing/strand_specific_mrna_sequencing/view-by-pid/41_Hf03_LiHe_Ct/replicate1/paired/run140827_SN751_0200_BC533PACXX/sequence)


TARGETFOLDER=/TL/deep/fhgfs/projects/pebert/thesis/projects/statediff/loaded_input/deep/rna_data
LOGFOLDER=/TL/deep/fhgfs/projects/pebert/thesis/projects/statediff/loaded_input/deep

# set appropriate username
USERNAME=

# set appropriate password
PASSWORD=

for SRCF in "${SRCFOLDERS[@]}"
do
    ASPERA_SCP_PASS=${PASSWORD} /TL/deep-share/archive00/software/bin/ascp -q -T -k1 -Q -P 33001 -l800M --overwrite=diff --symbolic-links=follow -L ${LOGFOLDER} ${USERNAME}@${SRCURL}:${SRCF} ${TARGETFOLDER}
done

exit $?
