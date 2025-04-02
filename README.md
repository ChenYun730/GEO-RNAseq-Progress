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
# 参考基因组 GRCh38.fa 构建索引（提交脚本并投递任务）
mkdir -p /mnt/alamo01/users/chenyun730/program/test/scripts
mkdir -p /mnt/alamo01/users/chenyun730/program/test/logs

#新建shell脚本
vim hisat2_index.sh
#! /usr/bin/sh
#micromamba activate R441
cd /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data
hisat2-build -p 64 GRCh38.fa GRCh38
2> /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/home_data

#赋予脚本执行权限
chmod 777 hisat2_index.sh

#提交任务(无需指定队列）
qsub -cwd -V -l cpu=64:mem=64G -q fast -N hisat2_index /mnt/alamo01/users/chenyun730/program/test/scripts/hisat2_index.sh
# qstat查看任务队列；qstat -f ID 查看任务； qdel 删除任务

# 检查索引文件
ls -lh /mnt/alamo01/users/chenyun730/program/test/compare/homo_sapiens/genome_index/

# 进行比对
#! /bin/bash
#source /mnt/alamo01/users/chenyun730/bin/micromamba
#micromamba activate R441
cd /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/
hisat2 -k 1 -p 64 \
  -x /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data/GRCh38 \
  -S SRR27961779.sam \
  --novel-splicesite-outfile SRR27961779_junction.bed \
  --no-unal --dta \
  --un-conc-gz SRR27961779_unmapped.fq.gz \
  -1 /mnt/alamo01/users/chenyun730/program/test/clean_data/SRR27961779_cleaned_1.fp.gz \
  -2 /mnt/alamo01/users/chenyun730/program/test/clean_data/SRR27961779_cleaned_2.fp.gz \
  2> SRR27961779.align.stats

#多个文件进行比对
#! /bin/bash
#source /mnt/alamo01/users/chenyun730/bin/micromamba
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
cd /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/
for SAMPLE in "${SAMPLES[@]}"; do
hisat2 -k 1 -p 64 \
  -x /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data/GRCh38 \
  -S ${SAMPLE}.sam \
  --novel-splicesite-outfile ${SAMPLE}_junction.bed \
  --no-unal --dta \
  --un-conc-gz ${SAMPLE}_unmapped.fq.gz \
  -1 /mnt/alamo01/users/chenyun730/program/test/clean_data/${SAMPLE}_cleaned_1.fp.gz \
  -2 /mnt/alamo01/users/chenyun730/program/test/clean_data/${SAMPLE}_cleaned_2.fp.gz \
  2> ${SAMPLE}.align.stats
echo "Submitted alignment job for ${SAMPLE}"
done
```
使用samtools和stringtie进行转录本组装和基因定量
```
#转换SAM为BAM并排序
$ vim samtools.sh
#! /bin/bash
#source /mnt/alamo01/users/chenyun730/bin/micromamba
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
cd /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/
for SAMPLE in "${SAMPLES[@]}"; do
echo "Processing ${SAMPLE}..."
  samtools sort -@ 8 -o ${SAMPLE}_sorted.bam ${SAMPLE}.sam
  samtools index ${SAMPLE}_sorted.bam
  echo "Finished processing ${SAMPLE}"
done

#使用stringtie进行转录本组装和基因定量
$ vim stringtie.sh
#! /bin/bash
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
        echo "Processing ${SAMPLE}..."
        stringtie /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/${SAMPLE}_sorted.bam \
  -G /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf \
  -o /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/${SAMPLE}.gtf \
  -p 64 \
  -B
        echo "Finished processing ${SAMPLE}"
done

 ls /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/*.gtf > gtf_list.txt
 stringtie --merge -G /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf   -o merged.gtf   gtf_list.txt
 
 #用sh.提交，计算样本的基因表达量（以合并后的基因组做参考）
#! /bin/bash
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing ${SAMPLE}..."
    stringtie /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/${SAMPLE}_sorted.bam \
        -G /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/merged.gtf \
         -o /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/quantify/${SAMPLE}.gtf \
         -e \
         -B \
         -A /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/quantify/${SAMPLE}_gene_abundance.tab
     echo "Finished processing ${SAMPLE}"
done
#下载prepDE.py3运行得到 .csv文件
$ python3 /mnt/alamo01/users/chenyun730/program/test/scripts/prepDE.py3 -i /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/quantify/quantify_gtf_list.txt -v

```

输出：比对结果（BAM）、定量矩阵


**思考题**：

为什么需要构建参考基因组索引？

1. 加速对比，减少对比时间；
   
2. 占用更少的内存；
   
3. 支持复杂比对模式：
     （a）剪接比对：RNA-seq数据需检测外显子连接，索引会预存剪切位点信息（HISAT2的snp_tran索引）；
     （b）突变容忍：包含SNP/突变的索引（如genome_snp_tran）可提高多态性样本的比对率。


# 🤷‍♀️未完成  
### 模块4：数据分析与可视化

**目标**：对定量数据进行质控并与GEO原作者的表达矩阵比较。

**任务**：

1. 表达矩阵的log2转换和标准化
```
#加载R包和.csv数据
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(FactoMineR)
library(factoextra)
>gene_counts <-read.csv("/mnt/alamo01/users/chenyun730/program/test/homo_sapiens/quantify/gene_count_matrix.csv", row.names = 1, header = TRUE)
>transcript_counts <- read.csv("/mnt/alamo01/users/chenyun730/program/test/homo_sapiens/quantify/transcript_count_matrix.csv", row.names = 1, header = TRUE)
>head(gene_counts)
>head(transcript_counts)

# log2 转换（加 1 避免 log(0)）
log_gene_counts <- log2(gene_counts + 1)
log_transcript_counts <- log2(transcript_counts + 1)

# Z-score 归一化
zscore_gene_counts <- t(scale(t(log_gene_counts)))
zscore_transcript_counts <- t(scale(t(log_transcript_counts)))

# 检查标准化结果
summary(zscore_gene_counts)
summary(zscore_transcript_counts)
```

2. 热图
```
# 基因水平热图
 # 仅保留至少 1 个非 NA 值的基因
filtered_counts <- zscore_gene_counts[rowSums(is.na(zscore_gene_counts)) == 0, ]
filtered_counts <- filtered_counts[rowSums(is.nan(filtered_counts)) == 0, ]
 pheatmap(filtered_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Gene Expression Heatmap ")

pheatmap(zscore_gene_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Gene Expression Heatmap")

#转录本水平热图
  # 仅保留至少 1 个非 NA 值的基因
filtered_counts_transcript <- zscore_transcript_counts[rowSums(is.na(zscore_transcript_counts)) == 0, ]
filtered_counts_transcript <- filtered_counts_transcript[rowSums(is.nan(filtered_counts_transcript)) == 0, ]

pheatmap(filtered_counts_transcript, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Transcript Expression Heatmap")
```
3.PCA分析
```
# 计算 PCA
gene_pca_res <- PCA(t(zscore_gene_counts), graph = FALSE)

# PCA 可视化
fviz_pca_ind(gene_pca_res,
             col.ind = "cos2", # 根据 cos2 贡献度着色
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE, # 避免标签重叠
             title = "PCA of Gene Expression")
# 计算 PCA
transcript_pca_res <- PCA(t(zscore_transcript_counts), graph = FALSE)

# PCA 可视化
fviz_pca_ind(transcript_pca_res,
             col.ind = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE,
             title = "PCA of Transcript Expression")


```
4. 输出：可视化结果（热图、PCA图）
```
# 保存基因水平热图
pdf("gene_heatmap.pdf", width = 8, height = 6)
pheatmap(zscore_gene_counts, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Gene Expression Heatmap")
dev.off()

# 保存转录本水平热图
pdf("transcript_heatmap.pdf", width = 8, height = 6)
pheatmap(zscore_transcript_counts, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Transcript Expression Heatmap")
dev.off()

# 保存基因水平 PCA 图
pdf("gene_PCA.pdf", width = 8, height = 6)
fviz_pca_ind(gene_pca_res, col.ind = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE, title = "PCA of Gene Expression")
dev.off()

# 保存转录本水平 PCA 图
pdf("transcript_PCA.pdf", width = 8, height = 6)
fviz_pca_ind(transcript_pca_res, col.ind = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE, title = "PCA of Transcript Expression")
dev.off()

```

**思考题**：

如果PCA中样本未能按组分离，原因可能是什么？
