# Project repository: Fast detection of differential chromatin domains with SCIDDO

## Unique identifier

This repository: [DOI: 10.17617/1.6K](https://doi.org/10.17617/1.6K)

bioRxiv preprint: [DOI: 10.1101/441766](https://doi.org/10.1101/441766)

SCIDDO manuscript: submitted

## Recreating figures

All figures in the publication can be recreated by running the relevant Jupyter notebook in `notebooks/plotting/active` with the following exceptions:

- figure 1: manually created

### Notebooks for figures

- figure 2: `plot_gumbel-fit_ka-params.ipynb`
- figure 3: `plot_compare_repscores.ipynb`
- figure 4: `plot_hsp_region_overlap.ipynb`
- figure 5: `plot_hsp-ovl-rgb_boxplot.ipynb`
- figure 6: `plot_hsp_body-enh_cdf.ipynb`
- figure 7: `plot_DEG_recovery_bars.ipynb`
- figure 8: `plot_hsp-gene-ovl_boxplot.ipynb`
- figure 9: `plot_hsp-t100-gene-ovl_bars.ipynb`
- figure 10: `plot_bgcov_DE-tpm-diff.ipynb`
- figure 11: `plot_uc1_p300_enh.ipynb`
- figure 12: `plot_hsp_pepr_perf.ipynb`

## The SCIDDO tool

Instruction how to install and use the SCIDDO tool are available here:

[SCIDDO repo: github.com/ptrebert/sciddo](https://github.com/ptrebert/sciddo)

## Steps to replicate all results

0) It is strongly recommended to run the analysis on a compute cluster with a DRMAA-compatible job scheduler. The expected runtime for the entire project is 2-3 days (assuming that a handful of servers is available to run the pipeline in parallel).
1) Create the software environment as specified in the Conda environment file (`/run_env/conda_env_statediff.yml`)
2) All pipelines of this project are implemented as [Ruffus pipelines](http://www.ruffus.org.uk/). It is recommended to install Ruffus to execute the pipelines. All (unformatted) command lines are also available in the `.ini` configuration files for each pipeline (see `/pipeline/ppl_statediff_main.ini`). Furthermore, you may want to use [PiedPiper](http://piedpiper.readthedocs.io) for convenience.
3) Download all DEEP raw data as listed in the Supplementary Material (see bioRxiv preprint). Note that this includes primary samples from human donors and access to these data files must by granted by the DEEP data access committee (see bioRxiv preprint for details).
4) Run complete pipeline `/pipeline/ppl_statediff_main`. For some of the annotation data, certain cleanup operations may be required. There are dedicated Jupyter notebooks to perform these operations in `notebooks/utils`.
5) After the pipeline has finished, you can execute the individual notebooks in `notebooks/plotting/active` to recreate the plots in the manuscript. Note that the notebooks commonly perform some caching of preprocessed data to speed up later layout changes of the plots, and thus may run for a couple of minutes upon the first execution.

