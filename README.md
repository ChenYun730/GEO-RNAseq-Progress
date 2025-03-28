# GEO-RNAseq-Progress
Daily Updates：Recording of project progress

## 任务背景：
基于指定的GEO编号GSE255647，通过RNA-seq数据的下载、预处理、比对、定量、差异表达分析和功能富集分析，完成转录组数据的全流程分析，并以Markdown形式记录过程，最终以GitHub仓库展示。

**任务要求**： 分析SARS-CoV-2以1 MOI感染Calu-3/2B4细胞系12h后，转录组水平的变化。

# ✔已完成
## 模块1：任务准备

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


## 模块2：数据下载与预处理

**目标**：从GEO下载原始数据，进行质控和清洗。

### 任务：

**1. 使用wget、curl或aspera下载指定GEO编号的FASTQ数据(以SRR27961778为例）**
```
#下载SRA文件
    wget -O SRR27961778.sra \"https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR27961778/SRR27961778"
#转换为FASTQ文件
   fastq-dump --split-files SRR27961778.sra
#可以存为压缩包
   gzip SRR27961779_1.fastq
```

**2. 质控：fastqc和multiq**
```
#生成原始数据质控报告
   fastqc SRR27961778_1.fastq.gz SRR27961778_2.fastq.gz -o ./raw_fastqc_results 
   cd /mnt/alamo01/users/chenyun730/program/test/raw_fastqc_results
   multiqc ./SRR*
```

**3. 清洗：fastp（去低质、去接头）**
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
#e ——平均质量阈值； q ——质量阈值（默认为15，这里提高到30为更严格）；
u ——允许达标碱基的比例为10%； r ——启用 polyG 修剪（适用于 NovaSeq）；
M ——polyG 的最小长度阈值；
W 5 ——滑动窗口大小（5bp）； w 64 ——使用64个线程；
f 10 F 10 ——切除 Read 1/2 前 10 bp
```

**4. 输出：质控报告、清洗后数据**
```
multiqc ./ -n "Cleaned_Data_QC"
```

**思考题**：

如果质控报告中碱基质量偏低，如何调整fastp的参数？
 
1.全局调整： 提高质量阈值q、降低允许未达标碱基的比例u、提高平均质量阈值e；
   
2.若5' 或 3' 端质量差，序列开头或末尾质量明显下降：提高首位切除的bp数（-f、-F、-t）；
   
3.若接头序列残留（MultiQC 的 Adapter Content 部分显示超标）：自动检测接头（即使已指定序列）或者额外切除5bp（-trim_front1 5);
   
4.若高比例的重复序列或低复杂度序列：过滤掉长度< 50bp 的 reads（默认 15）-l 50；移除低复杂度序列-low_complexity_filter；
   
5.若局部区域质量波动大：减小滑动窗口大小-W；滑动窗口平均质量阈值从 30 降至 25（更宽松）-M
   
# 🤷‍♀️未完成   

## 模块3：比对与定量

**目标**：基于参考基因组构建索引并进行比对与定量。

### 任务：

下载物种参考基因组和基因注释（GTF文件）
```
 #下载参考基因组的FASTQ和GTF文件
 cd /mnt/alamo01/users/chenyun730/program/test/compare
 mkdir homo_sapiens
 wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
 wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
 # 解压文件
 gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
 gunzip Homo_sapiens.GRCh38.104.gtf.gz
```

用hisat2-build构建索引，hisat2比对
```
#先构建索引再对比
hisat2-build -p 64 \
  /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/GRCh38.fa \     /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/genome_index/genome_index \                                  2> /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/index_build.log

# 检查索引文件
ls -lh /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/genome_index/

hisat2 -k 1 -p 64 -x /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/genome_index \
  -S SRR27961779.sam \
  --novel-splicesite-outfile SRR27961779_junction.bed \
  --no-unal --dta \
  --un-conc-gz SRR27961779_ummapped.fq.gz \
  -1 /mnt/alamo01/users/chenyun730/program/test/clean_data/SRR27961779_cleaned_1.fp.gz \
  -2 /mnt/alamo01/users/chenyun730/program/test/clean_data/SRR27961779_cleaned_2.fp.gz

#deepseek版本
hisat2 -k 1 -p 64 \
  -x /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/genome_index/genome_index \
  -S SRR27961779.sam \
  --novel-splicesite-outfile SRR27961779_junction.bed \
  --no-unal --dta \
  --un-conc-gz SRR27961779_unmapped.fq.gz \
  -1 /mnt/alamo01/users/chenyun730/program/test/clean_data/SRR27961779_cleaned_1.fp.gz \
  -2 /mnt/alamo01/users/chenyun730/program/test/clean_data/SRR27961779_cleaned_2.fp.gz \
  2> SRR27961779.align.stats

#转换SAM为BAM并排序
samtools sort -@ 16 \
  -o /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.sorted.bam \
  /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.sam

samtools index /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.sorted.bam
rm /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.sam  # 清理中间文件

#使用GTF文件进行转录本定量（StringTie）
stringtie /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.sorted.bam \
  -G /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/Homo_sapiens.GRCh38.104.gtf \
  -o /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.gtf \
  -p 64 \
  -B  # 生成Ballgown兼容文件

```

使用stringtie进行转录本组装和基因定量
```
stringtie /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.sorted.bam \
  -G /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/Homo_sapiens.GRCh38.104.gtf \
  -o /mnt/alamo01/users/chenyun730/program/test/alignment/SRR27961779.gtf \
  -p 64 \
  -B  # 生成Ballgown兼容文件

```

输出：比对结果（BAM）、定量矩阵


**思考题**：

为什么需要构建参考基因组索引？

1. 加速对比，减少对比时间；
   
2. 占用更少的内存；
   
3. 支持复杂比对模式：
     （a）剪接比对：RNA-seq数据需检测外显子连接，索引会预存剪切位点信息（HISAT2的snp_tran索引）；
     （b）突变容忍：包含SNP/突变的索引（如genome_snp_tran）可提高多态性样本的比对率。
