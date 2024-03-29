# pipelines2020
intro to population genomics

### 1. HPC intro, genotyping pipelines, PCA - *Acropora millepora* coral

[slides: TACC intro, genotyping](https://docs.google.com/presentation/d/1Po3J-SAM9Ju7l27Au2YeMmMK1d-0UBu6sLN_10VxWGI/edit?usp=sharing)

[slides: Pop structure](https://www.dropbox.com/s/l42knuvfsf3pif3/pop_structure.pptx?dl=0)

![*A. millepora* under fluorescent microscope](amillepora-redcyan.jpg)

Walkthrough: `pipelines_day1_QC_PCA_admixture.sh`

Scripts:
- `OKall_ibs.R`
- `OK_ibs.R`
- `admixturePlotting_pcangsd.R`
- `admixturePlotting_v5.R`

Files from `precomputed results` directory (in case you get stuck):
- `dd.pdf`: quality control plots
- `OKall_ibs.R` : initial IBS matrix
- `OKall_ibs.R` : final IBS matrix
- `bams.nc` : filtered bam list
- `pcangsd*` : several files containing `PCAngsd` output

### 2. Allele Frequency Spectra, demographic modeling

[slides](https://docs.google.com/presentation/d/1qvwG3MMP2xRPd4oGy6VyLXxzAKf99suKzCcMvUSFMrY/edit?usp=sharing)

Walkthrough: `pipelines_day2_AFS.sh`

### 3. Do-it-yourself! *Idotea baltica* project

![*Idotea baltica*](isopoda_idotea_balthica_01-10-15_1.jpg)

bams: `/work/01211/cmonstr/IB/recal_bams`

metadata: `https://www.dropbox.com/s/af5ggl5p7gtuccp/Idotea%20environment%20171005.csv?dl=0`
