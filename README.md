[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8354341.svg)](https://doi.org/10.5281/zenodo.8354341)

# RNAseq_workflow
This is a RNAseq workflow.
# Dependency
All the dependency and version are based on my current platform. The other version maybe compatible, but were untested.
## Basic
- perl v5.26.3
- Python 3.9.16
- R version 4.2.1
    - tidyverse_2.0.0
    - argparser_0.7.1
    - jsonlite_1.8.5
    - ggsci_3.0.0
    - cowplot_1.1.1
    - reshape2_1.4.4
    - rmarkdown_2.22
    - knitr_1.43
    - DT_0.28
## Prepare genomic data
- cufflinks-2.2.1
- Another GFF Analysis Toolkit (AGAT) - Version: v1.0.0
- HISAT2 version 2.2.1
- STAR 2.7.3a
## Filter raw sequencing data
- fastp 0.23.2
- FastQC v0.11.9
## Mapping
- HISAT2 version 2.2.1
- STAR 2.7.3a
- Python 3.9.16
    - RSeQC v5.0.1
## Quantitative
- R version 4.2.1
    - Rsubread_2.10.5
    - limma_3.52.4
    - edgeR_3.38.4
    - PCAtools_2.8.0
## DEG analysis
- R version 4.2.1
    - DESeq2_1.36.0
    - edgeR_3.38.4
    - ggrepel_0.9.3
    - ggtext_0.1.2 (optional)
## Enrichment
- R version 4.2.1
    - clusterProfiler_4.4.4
    - pathview_1.36.1
    - enrichplot_1.16.2
## Co-expression
- R version 4.2.1
    - WGCNA_1.72-1
# Preparation
## Sample information
Prepare `00.data/samples.txt`, this is a tab-separated file with four columns (groupName sampleName fq1 fq2), e.g.:
```
groupA	sampleA1	\<path to in1.fq of sampleA1>	\<path to in2.fq of sampleA1>
groupA	sampleA2	\<path to in1.fq of sampleA2>	\<path to in2.fq of sampleA2>
groupB	sampleB1	\<path to in1.fq of sampleB1>	\<path to in2.fq of sampleB1>
groupB	sampleB2	\<path to in1.fq of sampleB2>	\<path to in2.fq of sampleN2>
```
## Contrasts for DEG analysis
Prepare `04.DE_analysis/contrasts.txt`, this is a tab-separated file with two columns (treatment control), e.g.:
```
groupA	groupB
```
## Genome and annotation
- The genome file and a gff3 file should be exsited in `./db/` directory.
- For enrichment analysis, a R package `<orgdb>` for GO enrichment should be installed in `db/R_Library`, and `<organism>.kegg_info.RData` should be existed in `./db/`.
- `functionalAnnotation.txt` is a tab-separated file with first column `GeneID`.

# Usage
## Set variables
Some variables should be set, which is included in `./script/.conf`.
## Run RNAseq pipeline
If all the files and variables are prepared, execute `run_RNAseq.sh`.
```bash
cd ./script
nohup sh run_RNAseq.sh &
```
When the pipeline is finished without error, all result should be generated in corresponding directory.
## Generate report
There is a `report.Rmd` file, open it in Rstudio and click `knit`, you will get a analysis report in HTML format and a `result` directory containing all the result file. The gene functional information in `./db/functionalAnnotation.txt`was added to expression file and DEG file. You can packaging this files and delivery to your client.
```bash
tar zcvf RNAseq_result.tar.gz result/ report.html image/ libs/
```

# Citation
If you use this pipeline to processing transcriptome sequencing data, please cite:
> Wang Pengfei. (2023). laowang1992/RNAseq_workflow: a workflow for processing transcriptome sequencing data (1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8354341
