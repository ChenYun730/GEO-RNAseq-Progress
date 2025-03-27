# GEO-RNAseq-Progress
Daily Updates：Recording of project progress

## 任务背景：
基于指定的GEO编号GSE255647，通过RNA-seq数据的下载、预处理、比对、定量、差异表达分析和功能富集分析，完成转录组数据的全流程分析，并以Markdown形式记录过程，最终以GitHub仓库展示。

**任务要求**： 分析SARS-CoV-2以1 MOI感染Calu-3/2B4细胞系12h后，转录组水平的变化。

# ✔已完成
### 模块1：任务准备

**实用操作**

ls（英文全拼：list files）: 列出目录及文件名

cd（英文全拼：change directory）：切换目录

pwd（英文全拼：print work directory）：显示目前的目录

mkdir（英文全拼：make directory）：创建一个新的目录

rmdir（英文全拼：remove directory）：删除一个空的目录

cp（英文全拼：copy file）: 复制文件或目录

rm（英文全拼：remove）: 删除文件或目录

mv（英文全拼：move file）: 移动文件与目录，或修改文件与目录的名称



**完成**：

基于`micromamba`配置环境，安装并熟悉以下工具：

1. fastqc、multiqc、fastp：质控与清洗

2. hisat2、stringtie：比对与定量

3. samtools：BAM文件操作

4. R与DESeq2、clusterProfiler、ggplot2

5. 注册GitHub账号，掌握基本操作

**思考题**：

1. 如何在没有管理员权限下安装fastqc？

   答：有多种方法，包括micromamba、手动下载二进制文件以及在Java环境下通过源码编译等方式。这里使用micromamba、bioconda（或镜像）：
```
   micromamba create -p ~/micromamba_envs/fastqc_env -c bioconda -c conda-forge fastqc
```
2. GitHub的README文件有什么作用？
   
   项目的“门面”和“使用说明书”，像这样👀。

# 🤷‍♀️未完成

### 模块2：数据下载与预处理

**目标**：从GEO下载原始数据，进行质控和清洗。

**任务**：

1. 使用wget、curl或aspera下载指定GEO编号的FASTQ数据(以SRR27961778为例）
```
#下载SRA文件
    wget -O SRR27961778.sra \"https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR27961778/SRR27961778"
#转换为FASTQ文件
   fastq-dump --split-files SRR27961778.sra
#可以存为压缩包
   gzip SRR27961779_1.fastq
```

2. 质控：fastqc和multiq
```
#生成原始数据质控报告
   fastqc SRR27961778_1.fastq.gz SRR27961778_2.fastq.gz -o ./raw_fastqc_results 
   cd /mnt/alamo01/users/chenyun730/program/test/raw_fastqc_results
   multiqc ./SRR*
```

3. 清洗：fastp（去低质、去接头）
```
#先在test里新建一个文件夹clean——data用于存储清洗后的数据
$ cd /mnt/alamo01/users/chenyun730/program/test
$ mkdir clean_data
$ for F in /mnt/alamo01/users/chenyun730/program/test/raw_data/*_1.fastq.gz; do    
 R=${F%_*}_2.fastq.gz;
 BASE=${F##*/}; SAMPLE=${BASE%_*};
 time fastp -i $F -I $R
   -o /mnt/alamo01/users/chenyun730/program/test/clean_data/${SAMPLE}_cleaned_1.fp.gz
   -O /mnt/alamo01/users/chenyun730/program/test/clean_data/${SAMPLE}_cleaned_2.fp.gz
   -e 25 -q 30 -u 10 -r -M 30 -W 5 -w 64 -f 10 -F 10
   --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
   --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
   -h /mnt/alamo01/users/chenyun730/program/test/clean_data/${SAMPLE}_report.html
   -j /mnt/alamo01/users/chenyun730/program/test/clean_data/${SAMPLE}_report.json &
done
```
6. 输出：质控报告、清洗后数据
```
multiqc ./ -n "Cleaned_Data_QC"
```

**思考题**：

如果质控报告中碱基质量偏低，如何调整fastp的参数？
