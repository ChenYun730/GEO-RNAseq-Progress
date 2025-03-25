# GEO-RNAseq-Progress
Daily Updates：Recording of project progress

## 任务背景：
基于指定的GEO编号GSE255647，通过RNA-seq数据的下载、预处理、比对、定量、差异表达分析和功能富集分析，完成转录组数据的全流程分析，并以Markdown形式记录过程，最终以GitHub仓库展示。

**任务要求**： 分析SARS-CoV-2以1 MOI感染Calu-3/2B4细胞系12h后，转录组水平的变化。

# ✔已完成
### 模块1：任务准备

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

   micromamba create -p ~/micromamba_envs/fastqc_env **-c bioconda -c conda-forge** fastqc

2. GitHub的README文件有什么作用？
   
   项目的“门面”和“使用说明书”，像这样👀。

# 🤷‍♀️未完成

### 模块2：数据下载与预处理

**目标**：从GEO下载原始数据，进行质控和清洗。

**任务**：

1. 使用wget、curl或aspera下载指定GEO编号的FASTQ数据

2. 质控：fastqc和multiqc

3. 清洗：fastp（去低质、去接头）

4. 输出：质控报告、清洗后数据

**思考题**：

如果质控报告中碱基质量偏低，如何调整fastp的参数？
