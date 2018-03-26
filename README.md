# RNA-Seq 分析流程

该流程为全基因组RNA-Seq分析流程，主要包括Trim、去rRNA、Mapping、QC、定量（gene level 和 transcriptom level）、差异分析
；

路径：/lustre/rdi/user/songx/RNA-seq/bin/pipeline.sh
   
step1：prepare文件

运行此流程需提供bin/中example.txt格式文件，文件第一列为Patient ID，第二列为sample ID，第三列为样本的fastq文件绝对路径；

step2：running

bash bin/pipeline.sh foldchange值 pvalue值 qvalue值（若以pvalue为差异分析的阈值，qvalue则为NA）
