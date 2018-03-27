#!/usr/bin/bash

STAR_index=/lustre/rdi/user/songx/tools/genome/hg19/STARindex_150/
#STAR=/lustre/rdi/user/songx/tools/software/GTEx-Pipeline/STAR/source/STAR
STAR=/lustre/rdi/user/songx/tools/software/GTEx-Pipeline/gtex-pipeline-master/rnaseq/src/run_STAR.py
GTF=/lustre/rdi/user/songx/tools/genome/hg19/gencode/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
rRNA=/lustre/rdi/user/songx/tools/genome/hg19/rRNA/RSeQC.hg19_rRNA.bed
ref_bed=/lustre/rdi/user/songx/tools/genome/hg19/hg19_RefSeq.bed
rsem_reference=/lustre/rdi/user/songx/tools/genome/hg19/rsem_reference_gencode/rsem_reference_gencode
QCmerge=/lustre/rdi/user/songx/RNA-seq/bin/mapstat.py
RSeQC=/lustre/rdi/user/songx/tools/anaconda2/bin
RSEM=/lustre/rdi/user/songx/tools/software/GTEx-Pipeline/RSEM-1.3.0
Trim=/lustre/rdi/user/songx/tools/software/GTEx-Pipeline/Trimmomatic-0.36/trimmomatic-0.36.jar
adapter=/lustre/rdi/user/songx/tools/software/GTEx-Pipeline/Trimmomatic-0.36/adapters
bowtie2=/lustre/rde/user/guoxc/software/bin/bowtie2
rRNA=/lustre/rdi/user/songx/tools/Database/rRNA/rRNA
DESeq2=/lustre/rdi/user/songx/RNA-seq/R/DESeq2.R

cat $PWD/example.txt |while read PID SID dir; do
sample=$PID"_"$SID

mkdir -p $PWD/Mapping/$sample
mkdir -p $PWD/QC/$sample
mkdir -p $PWD/Quantification/
mkdir -p $PWD/Trimmomatic/$sample/rRNA

##Trimmomatic (raw2clean)
##phred33 :设置碱基的质量格式,如果不设置，默认的是-phred64;
#trimlog file :就是产生日志，包括如下部分内容：read的名字，留下来的序列的长度，第一个碱基的起始位置，从开始trimmed的长度，最后的一个碱基位于初始read的位置，trimmed掉的碱基数量;
#LEADING :3 切除首端碱基质量小于3的碱基;
#TRAILING :3 切除末端碱基质量小于3的碱基;
#ILLUMINACLIP: 1.adapter.lis:2:30:10 1.adapter.list为adapter文件，允许的最大mismatch数，palindrome模式下匹配碱基数阈值：simple模式下的匹配碱基数阈值;
#SLIDINGWINDOW:4:15  Windows的size是4个碱基，其平均碱基质量小于15，则切除;
#MINLEN:36  最低reads长度为36;
#CROP:  保留的reads长度;
#HEADCROP: 在reads的首端切除指定的长度。
java -jar $Trim PE -threads 4 -phred33 -trimlog $PWD/Trimmomatic/$sample/trim.log $dir/$SID"_"1.fastq.gz $dir/$SID"_"2.fastq.gz $PWD/Trimmomatic/$sample/$sample"_"1.paired.fq.gz $PWD/Trimmomatic/$sample/$sample"_"1.unpaired.fq.gz $PWD/Trimmomatic/$sample/$sample"_"2.paired.fq.gz $PWD/Trimmomatic/$sample/$sample"_"2.unpaired.fq.gz ILLUMINACLIP:$adapter/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

mv $PWD/Trimmomatic/$sample/$sample"_"1.paired.fq.gz $PWD/Trimmomatic/$sample/$sample"_"1.clean.fq.gz
mv $PWD/Trimmomatic/$sample/$sample"_"2.paired.fq.gz $PWD/Trimmomatic/$sample/$sample"_"2.clean.fq.gz

##remove rRNA
less $PWD/Trimmomatic/$sample/$sample"_"1.clean.fq.gz |head -n 4000000 |gzip > $PWD/Trimmomatic/$sample/rRNA/$sample"_"1.test.fq.gz
less $PWD/Trimmomatic/$sample/$sample"_"2.clean.fq.gz |head -n 4000000 |gzip > $PWD/Trimmomatic/$sample/rRNA/$sample"_"2.test.fq.gz
$bowtie2 -x $rRNA -1 $PWD/Trimmomatic/$sample/rRNA/$sample"_"1.test.fq.gz -2 $PWD/Trimmomatic/$sample/rRNA/$sample"_"2.test.fq.gz --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 -S $PWD/Trimmomatic/$sample/rRNA/$sample"_"test.sam 2> $PWD/Trimmomatic/$sample/rRNA/$sample.rRNA.stat.txt
rm $PWD/Trimmomatic/$sample/rRNA/$sample"_"1.test.fq.gz $PWD/Trimmomatic/$sample/rRNA/$sample"_"2.test.fq.gz $PWD/Trimmomatic/$sample/rRNA/$sample"_"test.sam

##判断rRNA含量,如果rRNA含量大于等于1%，则剔除rRNA；
rRNA_rate=`tail -n1 $PWD/Trimmomatic/$sample/rRNA/$sample.rRNA.stat.txt|awk -F '%' '{print $1}'`
rRNA_rate=`echo $rRNA_rate/1|bc`
if [ $rRNA_rate -ge 1 ]; then
    $bowtie2 -x $rRNA -1 $PWD/Trimmomatic/$sample/$sample"_"1.clean.fq.gz -2 $PWD/Trimmomatic/$sample/$sample"_"2.clean.fq.gz --end-to-end --sensitive -p 8 --phred33 --un-conc-gz $PWD/Trimmomatic/$sample/rRNA/$sample.unmapped.gz -X 600 -S $PWD/Trimmomatic/$sample/rRNA/$sample.sam 2> $PWD/Trimmomatic/$sample/rRNA/$sample.rRNA.stat.txt
    mv $PWD/Trimmomatic/$sample/rRNA/$sample.unmapped.1.gz $PWD/Trimmomatic/$sample/$sample"_"1.clean.fq.gz
    mv $PWD/Trimmomatic/$sample/rRNA/$sample.unmapped.2.gz $PWD/Trimmomatic/$sample/$sample"_"2.clean.fq.gz
    rm $PWD/Trimmomatic/$sample/rRNA/$sample.sam 
fi
  
##STAR mapping
python $STAR $STAR_index $PWD/Trimmomatic/$sample/$sample"_"1.clean.fq.gz $PWD/Trimmomatic/$sample/$sample"_"2.clean.fq.gz $sample -o Mapping/$sample --annotation_gtf $GTF
rm -rf $PWD/Mapping/$sample/$sample._STARpass1 

##用STAR进行基因定量（HTseq-count）
less $PWD/Mapping/$sample/$sample.ReadsPerGene.out.tab |sed '1,4d'|cut -f 1,2|awk -v var="$sample" '{print "GeneID\t"var"\n"$0}'|sort -k1 -n -r|uniq > $PWD/Quantification/$sample/$sample.star.readcount.txt
paste $PWD/Quantification/*/*.star.readcount.txt |awk -F '\t' '{for(i=1;i<=NF;i++){if (i<=1){printf $1"\t"}else if(i%2==0){printf $i"\t"}}printf "\n"}' $PWD/Quantification/readcount.txt >$PWD/Quantification/readcount.txt

##用RSEM对transcriptome进行定量（可选）
$RSEM/rsem-calculate-expression --alignments --num-threads 4 --estimate-rspd  --forward-prob 0.5 --fragment-length-max 1000 --paired-end $PWD/Mapping/$sample/$sample.Aligned.toTranscriptome.out.bam $rsem_reference $sample -o $PWD/Quantification/$sample/$sample.rsem
cut -f 1,2,5 $PWD/Quantification/$sample/$sample.rsem.isoforms.results > $PWD/Quantification/$sample/$sample.rsem.isforms.txt
paste $PWD/Quantification/*/*.rsem.isforms.txt |awk -F '\t' '{for(i=1;i<=NF;i++){if (i<=2){printf $1"\t"$2"\t"}else if(i%3==0){printf $i"\t"}}printf "\n"}' >$PWD/Quantification/isform.readcount.txt

##QC_Merge
python $QCmerge --rootdir $PWD

###RSeQC 
#Calculate the distributions of clipped nucleotides across reads
#NOTE:
#out.clipping_profile.r is the R script file used to generate pdf file(s).
#out.clipping_profile.xls: contains 3 columns: 
#the first column is position (starting from 0) of read in 5’->3’ direction; 
#the second columnis the number of reads clipped at this position; 
#the third column is the number of reads non-clipped at this position.
$RSeQC/clipping_profile.py -i $PWD/Mapping/$sample/$sample"_"mapped_ref.bam -o $PWD/QC/$sample/$sample -s PE

#Calculate inner distance between read pairs
#This module is used to calculate the inner distance (or insert size) between two paired RNA reads. The distance is the mRNA length between two paired fragments. We first determine the genomic (DNA) size between two paired reads: D_size = read2_start - read1_end
#NOTE:
#if two paired reads map to the same exon: inner distance = D_size
#if two paired reads map to different exons:inner distance = D_size - intron_size
#if two paired reads map non-exonic region (such as intron and intergenic region): inner distance = D_size
#The iner_distance might be a negative value if two fragments were overlapped.
$RSeQC/inner_distance.py -i $PWD/Mapping/$sample/$sample"_"mapped_ref.bam -o $PWD/QC/$sample/$sample -r $ref_bed

##GC 
$RSeQC/read_GC.py -i $PWD/Mapping/$sample/$sample"_"mapped_ref.bam -o $PWD/QC/$sample/$sample
$RSeQC/read_NVC.py -i $PWD/Mapping/$sample/$sample"_"mapped_ref.bam -o $PWD/QC/$sample/$sample
$RSeQC/read_quality.py -i $PWD/Mapping/$sample/$sample"_"mapped_ref.bam -o $PWD/QC/$sample/$sample

#read_distribution
#Total Reads: This does NOT include those QC fail,duplicate and non-primary hit reads
#Total Tags: reads spliced once will be counted as 2 tags, reads spliced twice will be counted as 3 tags, etc. And because of this, “Total Tags” >= “Total Reads”
#Total Assigned Tags: number of tags that can be unambiguously assigned the 10 groups (see below table).
#Tags assigned to “TSS_up_1kb” were also assigned to “TSS_up_5kb” and “TSS_up_10kb”, tags assigned to “TSS_up_5kb” were also assigned to “TSS_up_10kb”. Therefore, “Total Assigned Tags” = CDS_Exons + 5’UTR_Exons + 3’UTR_Exons + Introns + TSS_up_10kb + TES_down_10kb.
#When assign tags to genome features, each tag is represented by its middle point.
$RSeQC/read_distribution.py -i $PWD/Mapping/$sample/$sample"_"mapped_ref.bam -r $ref_bed > $PWD/QC/$sample/$sample"_"region_stat.txt

##对transcriptome进行定量
$RSEM/rsem-calculate-expression --alignments --num-threads 4 --estimate-rspd  --forward-prob 0.5 --fragment-length-max 1000 --paired-end $PWD/Mapping/$sample/$sample.Aligned.toTranscriptome.out.bam $rsem_reference $sample -o $PWD/Quantification/$sample/$sample.rsem
cut -f 1,2,5 $PWD/Quantification/$sample/$sample.rsem.isoforms.results > $PWD/Quantification/$sample/$sample.rsem.isforms.txt
paste $PWD/Quantification/*/*.rsem.isforms.txt |awk -F '\t' '{for(i=1;i<=NF;i++){if (i<=2){printf $1"\t"$2"\t"}else if(i%3==0){printf $i"\t"}}printf "\n"}' >$PWD/Quantification/isform.readcount.txt

done

###DESeq2差异分析
genename=/lustre/rdi/user/songx/tools/genome/hg19/gencode/gencode.V19.genename.txt
foldchange=$1
pvalue=$2
padj=$3

sed '1d' $PWD/compare.txt |while read compare design; do

output=$PWD/Differential/$compare

if [ ! -d "$output" ];then
  mkdir -p $output
fi

Rscript $DESeq2 --rawcount $PWD/Quantification/readcount.txt --condition $PWD/condition.txt --compare $compare --design $design --foldchange $foldchange --pvalue $pvalue --padj $padj --genename $genename --outdir $PWD/Differential

done
