# RNAseq_workflow
This is a RNAseq workflow.
# 一 分析流程图
![image](./image/转录组分析流程图.png)
# 二 分析流程及结果
## 1	取样、mRNA提取、建库及测序
参考公司报告。
## 2	数据过滤
使用fastp<sup>[[1](#ref)]</sup>（version: 0.20.0）对raw data进行过滤得到clean data。统计过滤前后total bases、total reads、Q30、Q20、GC content以及有效数据比率（data_stat.csv/txt），同时使用FastQC（version: 0.11.9）对过滤前后的数据进行质量评估（QC/sample_fastqc.html）。
## 3	比对到参考基因组
使用HISAT2<sup>[[2](#ref)]</sup>（version: 2.1.0）将clean reads比对到甘蓝型油菜ZS11参考基因组<sup>[[3](#ref)]</sup>上，得到SAM（Sequence Alignment/Map）格式文件，然后使用SAMtools（version: 1.9）对比对结果（SAM文件）按照染色体和位置进行排序并转换为BAM（Binary Alignment/Map）格式文件<sup>[[4](#ref)]</sup>，可以将BAM文件导入IGV（Integrative Genomics Viewer）<sup>[[5](#ref)]</sup>对比对结果进行可视化。HISAT2可以使用更少资源的同时具有更快的速度，HISAT2比对时，对于非链特异性文库使用默认参数，链特异性文库需要指定文库类型（first使用--rna-strandness RF，second使用--rna-strandness FR）。比对完成后，我们对比对结果进行评估，统计比对率和唯一比对率。
## 4	表达量计算
根据比对结果（BAM文件），我们使用R（version: 4.0.2）软件的扩展包Rsubread<sup>[[6](#ref)]</sup>（version: 2.2.6）中的featureCounts函数计算每个基因的表达量（read count）并进行归一化处理（normalization），得到TPM（Transcripts Per Kilobase of exon model per Million mapped reads）和TMM（trimmed mean of M value）表达矩阵。
## 5	差异表达分析
根据基因表达矩阵（read count）文件，在有生物学重复的情况下，使用R（version: 4.0.2）软件的扩展包DESeq2<sup>[[7](#ref)]</sup>（version: 1.28.1）进行差异表达分析，在没有生物学重复的情况下，则使用R扩展包edgeR<sup>[[8](#ref)]</sup>（version: 3.30.3）进行差异表达分析，并推荐测生物学重复，也不算太贵。对于log2FoldChange绝对值大于1，并且padj小于0.05的基因则认为是差异表达基因（阈值需根据实际情况做出调整）。
## 6	功能富集分析
根据筛选出的差异表达基因，我们使用R（version: 4.0.2）软件的扩展包clusterProfiler<sup>[[9](#ref)]</sup>（version: 3.16.1）依据超几何分布检验来完成GO和KEGG富集分析（Over-representation analysis），设置参数pvalueCutoff和qvalueCutoff为0.05筛选显著富集的GO/KEGG term。在绘图时，如果GO或KEGG的term太多（一般是GO），建议取前10或15个term（GO中CC、BP和MF各选10或15各）进行绘图，如果相关term不在前10或15个内，也可以手动添加。<br/>
同时，我们使用clusterProfiler<sup>[[9](#ref)]</sup>进行GSEA（Gene Set Enrichment Analysis），其基本思想是将基因按照两组样本中差异表达程度排序，使用预先定义的基因集（GO或KEGG），检验预先定义的基因集是否排列在顶端或底端。GSEA可以充分利用基因差异表达程度的信息，排除人为筛选差异表达基因的主观性，是一种更先进富集方法，与over-representation富集方法互为补充。
## 7	WGCNA
加权基因共表达网络分析（WGCNA，Weighted correlation network analysis）可以用来鉴定样本间高度协同变化的基因集（模块），同时可以根据模块特征值（eigengene）将模块与外部性状信息相关联，以此鉴定与性状相关的模块并进一步挖掘关键基因<sup>[[10](#ref),[11](#ref)]</sup>。进行WGCNA至少需要15个样本，最好是20个及以上。在这里我们筛选差异表达基因使用R（version: 4.0.2）软件的扩展包WGCNA<sup>[[11](#ref)]</sup>（version: 1.69）进行基因模块的构建以及模块-样本、模块-形状的关联（其中各步骤参数均需根据实际情况决定）。
<br/>
<br/>

---
<div id="ref"></div>

#  参考文献
- [1] Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu. fastp: an ultra-fast all-in-one FASTQ preprocessor[J]. Bioinformatics, 2018, 34(17).
- [2] Daehwan Kim, Joseph M. Paggi, Chanhee Park, Christopher Bennett, Steven L. Salzberg. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype[J]. Nature Biotechnology: The Science and Business of Biotechnology, 2019, 37(8).
- [3] Jia-Ming Song, Zhilin Guan, Jianlin Hu, Chaocheng Guo, Zhiquan Yang, Shuo Wang, Dongxu Liu, Bo Wang, Shaoping Lu, Run Zhou, Wen-Zhao Xie, Yuanfang Cheng, Yuting Zhang, Kede Liu, Qing-Yong Yang, Ling-Ling Chen, Liang Guo. Eight high-quality genomes reveal pan-genome architecture and ecotype differentiation of Brassica napus[J]. Nature Plants, 2020, 6(1).
- [4] Li Heng, Handsaker Bob, Wysoker Alec, Fennell Tim, Ruan Jue, Homer Nils, Marth Gabor, Abecasis Goncalo, Durbin Richard. The Sequence Alignment/Map format and SAMtools.[J]. Bioinformatics (Oxford, England), 2009, 25(16).
- [5] Robinson James T, Thorvaldsdóttir Helga, Winckler Wendy, Guttman Mitchell, Lander Eric S, Getz Gad, Mesirov Jill P. Integrative genomics viewer.[J]. Nature biotechnology, 2011, 29(1).
- [6] Liao Yang, Smyth Gordon K, Shi Wei. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads[J]. Nucleic Acids Research, 2019, 47(8).
- [7] Michael I Love, Wolfgang Huber, Simon Anders. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2[J]. Genome Biology, 2014, 15(12).
- [8]	Mark D. Robinson, Davis J. McCarthy, Gordon K. Smyth. edgeR : a Bioconductor package for differential expression analysis of digital gene expression data[J]. Bioinformatics, 2010, 26(1).
- [9]	Yu Guangchuang, Wang Li-Gen, Han Yanyan, He Qing-Yu. clusterProfiler: an R package for comparing biological themes among gene clusters.[J]. Omics : a journal of integrative biology, 2012, 16(5).
- [10]	Bin Zhang, Steve Horvath. A General Framework for Weighted Gene Co-Expression Network Analysis[J]. Statistical Applications in Genetics and Molecular Biology, 2005,4(1).
- [11]	Peter Langfelder, Steve Horvath. WGCNA: an R package for weighted correlation network analysis[J]. BMC Bioinformatics, 2008, 9(2).
