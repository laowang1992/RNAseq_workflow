#  run_RNAseq.sh
#  
#  Copyright 2021 WangPF <wangpf0608@126.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
#!/bin/bash
# Program:
# 	RNAseq workflow
# History:
# 	20210118	First release
# Author:
# 	WangPF

## 加载配置文件
. ./.conf

# 默认设置会导致STAR align失败，因此修改设置
ulimit -n 1024000

#########################################################################
## 准备
cd ${work_dir}/db

gffread ${gff} -T -o ${gtf}
agat_convert_sp_gff2bed.pl --gff ${gff} --out ${bed}

if [ $aligner = STAR ]; then
	mkdir -p ${index}STARindex
	STAR --runThreadN $thread \
	     --runMode genomeGenerate \
	     --genomeDir ${index}STARindex \
	     --genomeFastaFiles ${genome} \
	     --sjdbGTFfile $gtf
elif [ $aligner = hisat2 ]; then
	hisat2-build -p ${thread} ${genome} ${index}
fi

#########################################################################

cat ${sampleInfo} | while read group sample fq1 fq2
do

# 过滤、质控
mkdir -p ${work_dir}/00.data/01.clean_data/QC
cd ${work_dir}/00.data/01.clean_data

fastp -i ${fq1} -o ./${sample}_1.clean.fastq.gz \
      -I ${fq2} -O ./${sample}_2.clean.fastq.gz \
      --json=./${sample}.json --html=${sample}.html --report_title="${sample} fastp report" \
      --thread=${thread}

# 比对
mkdir -p ${work_dir}/01.Mapping
cd ${work_dir}/01.Mapping

if [ $aligner = STAR ]; then
	mkdir -p ${sample}
	STAR --twopassMode Basic \
	     --runThreadN $thread --genomeDir ${index}STARindex \
	     --outSAMtype BAM SortedByCoordinate \
	     --outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA \
	     --outSAMmapqUnique 60 \
	     --outSAMmultNmax 1 \
	     --outFilterMismatchNoverReadLmax 0.04 --outSJfilterReads Unique \
	     --outFileNamePrefix ./$sample/$sample \
	     --readFilesCommand gunzip -c \
	     --readFilesIn $work_dir/00.data/01.clean_data/${sample}_1.clean.fastq.gz $work_dir/00.data/01.clean_data/${sample}_2.clean.fastq.gz
	ln -sf ./$sample/${sample}Aligned.sortedByCoord.out.bam ./$sample.sort.bam
	sambamba index --nthreads=${thread} ./$sample.sort.bam
	
elif [ $aligner = hisat2 ]; then
	if [ $strandSpecific = 0 ]; then
		hisat2 --new-summary -p ${thread} \
		       -x ${index} \
		       -1 ../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
		       -2 ../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		       -S ${sample}.sam \
		       1> ${sample}.log 2>&1
	elif [ $strandSpecific = 1 ]; then
		hisat2 --new-summary -p ${thread} \
		       -x ${index} --rna-strandness RF \
		       -1 ../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
		       -2 ../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		       -S ${sample}.sam \
		       1> ${sample}.log 2>&1
	elif [ $strandSpecific = 2 ]; then
		hisat2 --new-summary -p ${thread} \
		       -x ${index} --rna-strandness FR \
		       -1 ../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
		       -2 ../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		       -S ${sample}.sam \
		       1> ${sample}.log 2>&1
	fi
	
	#samtools sort -@ ${thread} -O BAM -o ${sample}.sort.bam ${sample}.sam
	# 使用sambamba代替samtools排序
	sambamba view --format=bam --with-header --sam-input --nthreads=${thread} --output-filename ${sample}.bam ${sample}.sam
	sambamba sort -t ${thread} -m 20GB --tmpdir=./ -o ${sample}.sort.bam ${sample}.bam
	rm ${sample}.sam ${sample}.bam
fi

# 定量
mkdir -p  ${work_dir}/02.Quantification
cd ${work_dir}/02.Quantification
Rscript run-featurecounts.R \
	-b ../01.Mapping/${sample}.sort.bam \
	-g ${gtf} -o ${sample} \
	--nthread ${thread} \
	--featureType ${featureType} \
	--attrType ${attrType} \
	--strandSpecific ${strandSpecific}
done


cd ${work_dir}/03.Merge_result
chmod u+x script/support_scripts/run_TMM_scale_matrix.pl
cut -f2 ${sampleInfo} | sed 's/^/..\/02.Quantification\//' | sed 's/$/.count/' > genes.quant_files.txt
perl script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes

# 差异表达
cd ${work_dir}/04.DE_analysis
#perl ${trinity}/Analysis/DifferentialExpression/run_DE_analysis.pl \
#	--matrix ../03.Merge_result/genes.counts.matrix \
#	--method ${de_method} \
#	--samples_file ../00.data/samples.txt \
#	--min_reps_min_cpm 2,1 \
#	--output ./ \
#	--contrasts contrasts.txt

# contrasts.txt: one contrast one line, "treatment	control"
if [ $de_method = DESeq2 ]; then
	Rscript run_DESeq2.R \
		--matrix ../03.Merge_result/genes.counts.matrix \
		--samples_file ../00.data/samples.txt \
		--min_reps 2 \
		--min_cpm 1 \
		--contrasts ./contrasts.txt
elif [ $de_method = edgeR ]; then
	Rscript run_edgeR.R \
		--matrix ../03.Merge_result/genes.counts.matrix \
		--samples_file ../00.data/samples.txt \
		--min_reps 1 \
		--min_cpm 1 \
		--dispersion $dispersion \
		--contrasts ./contrasts.txt
fi
# upset plot
#Rscript upset.R --de_log2FoldChange ${de_log2FoldChange} --de_padj ${de_padj}

# volcano plot
for i in ./*DE_results.txt
do
Rscript Volcano_plot.R \
	--de_result ${i} \
	--padj_cutoff ${de_padj} \
	--log2FC_cutoff ${de_log2FoldChange} \
	--heatmap TRUE
done
	
# 富集分析
cd ${work_dir}/05.ORA
for de_result in ${work_dir}/04.DE_analysis/*DE_results.txt
do
Rscript enrich.R --de_result ${de_result} \
	--de_log2FoldChange ${de_log2FoldChange} \
	--de_padj ${de_padj} \
	--enrich_pvalue ${enrich_pvalue} \
	--enrich_qvalue ${enrich_qvalue} \
	--orgdb ${orgdb} \
	--species ${species}
done

cd ${work_dir}/06.GSEA
for de_result in ${work_dir}/04.DE_analysis/*DE_results.txt
do
Rscript GSEA.R --de_result ${de_result} \
	--gsePadj ${gsePadj} \
	--orgdb ${orgdb} \
	--drawPdf ${pdf} \
	--species ${species}
done


## 统计信息
#
cd ${work_dir}/00.data/00.raw_data
fastqc -o ./QC --nogroup --threads ${thread} `cut -f3,4 ../samples.txt`
cd ${work_dir}/00.data/01.clean_data
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
# 
cd ${work_dir}/00.data/
Rscript dataStat.R

cd ${work_dir}/01.Mapping
if [ $aligner = STAR ]; then
	perl STARalignStat.pl ${sampleInfo} ./ align_stat.tsv
elif [ $aligner = hisat2 ]; then
	echo -e "Sample	Total read	Mapping read	Mapping rate	Unique mapping read	Unique mapping rate" > align_stat.tsv
	for i in $(cut -f2 ${sampleInfo}); do perl Hisat2alignStat.pl $i; done >> align_stat.tsv
fi

# 将所有bam一起geneBody_coverage.py实在是太慢了
#awk '{print $2".sort.bam"}' ${sampleInfo} > bamFiles.txt
#geneBody_coverage.py -r ${bed} -i bamFiles.txt -o geneBody_coverage
# 由所有bam一起geneBody_coverage.py改为分开并行以提升速度
#awk '{print $2}' ${sampleInfo} | \
#	parallel -j ${thread} -I% --max-args 1 \
#	geneBody_coverage.py -r ${bed} -i %.sort.bam -o %

#awk '{print $2}' ${sampleInfo} | \
#	parallel -j ${thread} -I% --max-args 1 \
#	RPKM_saturation.py -r ${bed} -i %.sort.bam -o %
#awk '{print $2}' ${sampleInfo} | \
#	parallel -j ${thread} -I% --max-args 1 \
#	pdftoppm -singlefile -r 500 -png %.saturation.pdf %.saturation

# 并行统计覆盖度和饱和度太占内存了，改为单个计算。20230807
cat ${sampleInfo} | while read group sample fq1 fq2
do
	geneBody_coverage.py -r ${bed} -i $sample.sort.bam -o $sample &
	RPKM_saturation.py -r ${bed} -i $sample.sort.bam -o $sample &
	wait
done

for file in `ls *pdf`
do
	pdftoppm -singlefile -r 500 -png $file ${file%.pdf*}
done

cd ${work_dir}/02.Quantification
Rscript QuantificationStat.R

cd ${work_dir}/03.Merge_result
Rscript stat_plot.R
