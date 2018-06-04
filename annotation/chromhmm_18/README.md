# Annotation data for ChromHMM 18 state model

This model has been trained and used by the [ROADMAP initiative](http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state).
It is based on 6 histone modifications plus Input control:

  - H3K4me3
  - H3K4me1
  - H3K27ac
  - H3K27me3
  - H3K9me3
  - H3K36me3

## Downloaded files

The following set of files has been downloaded from the above stated reference.

File containing the model specification:
[model_18_core_K27ac.txt](http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/model_18_core_K27ac.txt)

File containing short labels to be displayed in a genome browser: 
[browserlabelmap_18_core_K27ac.tab](http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/browserlabelmap_18_core_K27ac.tab)

File containing state emission probabilities:
[emissions_18_core_K27ac.txt](http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/emissions_18_core_K27ac.txt)

File containing RGB colors for each state:
[colormap_18_core_K27ac.tab](http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/colormap_18_core_K27ac.tab)

## Manually created files

File containing scoring (state dissimilarity) scheme created for this project:

chromhmm_18_expert_scoring.tsv