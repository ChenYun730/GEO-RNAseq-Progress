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
for SAMPLE in "${SAMPLES[@]}"; do
echo "Processing ${SAMPLE}..."
  mv /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/alignment/${SAMPLE}.sam /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}
  cd /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}
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
 cd /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}
        stringtie /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}/${SAMPLE}_sorted.bam \
  -G /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf \
  -o /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}.gtf \
  -p 64 \
  -B
        echo "Finished processing ${SAMPLE}"
done

find /mnt/alamo01/users/chenyun730/program/test/alignment/ -type f -name "*.gtf" > gtf_list.txt
 stringtie --merge -G /mnt/alamo01/users/chenyun730/program/test/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf
-o merged.gtf   gtf_list.txt
 
 #用sh.提交，计算样本的基因表达量（以合并后的基因组做参考）
#! /bin/bash
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing ${SAMPLE}..."
    stringtie /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}/${SAMPLE}_sorted.bam \
        -G /mnt/alamo01/users/chenyun730/program/test/alignment/merged.gtf \
         -o /mnt/alamo01/users/chenyun730/program/test/quantify/${SAMPLE}.gtf \
         -e \
         -B \
         -A /mnt/alamo01/users/chenyun730/program/test/quantify/${SAMPLE}/${SAMPLE}_gene_abundance.tab
     echo "Finished processing ${SAMPLE}"
done
```

提取FPKM到一个矩阵
```
vim fpkm.sh
#! /bin/bash
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
    awk 'BEGIN {OFS="\t"} {if (NR>1) print $9, $6, $12}' /mnt/alamo01/users/chenyun730/program/test/quantify/${SAMPLE}/${SAMPLE}/t_data.ctab > /mnt/alamo01/users/chenyun730/program/test/quantify/${SAMPLE}/${SAMPLE}_FPKM.txt
done
OUTPUT_FILE="/mnt/alamo01/users/chenyun730/program/test/quantify/all_samples_FPKM.txt"
echo -n "gene_id" > $OUTPUT_FILE
for SAMPLE in "${SAMPLES[@]}"; do
    echo -n -e "\tFPKM_${SAMPLE}" >> $OUTPUT_FILE
done
echo "" >> $OUTPUT_FILE
FIRST_SAMPLE="${SAMPLES[0]}"
FIRST_FILE="/mnt/alamo01/users/chenyun730/program/test/quantify/${FIRST_SAMPLE}/${FIRST_SAMPLE}_FPKM.txt"
cut -f1,2 "$FIRST_FILE" | tail -n +2 > /tmp/combined.tmp
for ((i=1; i<${#SAMPLES[@]}; i++)); do
    SAMPLE="${SAMPLES[$i]}"
    SAMPLE_FILE="/mnt/alamo01/users/chenyun730/program/test/quantify/${SAMPLE}/${SAMPLE}_FPKM.txt"
    cut -f2 "$SAMPLE_FILE" | tail -n +2 > "/tmp/${SAMPLE}_col.tmp"
    paste /tmp/combined.tmp "/tmp/${SAMPLE}_col.tmp" > /tmp/combined_new.tmp
    mv /tmp/combined_new.tmp /tmp/combined.tmp
done
cat /tmp/combined.tmp >> "$OUTPUT_FILE"
rm /tmp/*_col.tmp /tmp/combined.tmp
echo "FPKM data for all samples has been merged into $OUTPUT_FILE"

#下载prepDE.py3运行得到 .csv文件（使用count可进行）这里不用
python3 /mnt/alamo01/users/chenyun730/program/test/scripts/prepDE.py3 -i /mnt/alamo01/users/chenyun730/program/test/quantify/quantify_gtf_list.txt -v
python3 /mnt/alamo01/users/chenyun730/program/test/scripts/prepDE.py3 -i .  -v  #最佳版本

```

输出：比对结果（BAM）、定量矩阵


**思考题**：

为什么需要构建参考基因组索引？

1. 加速对比，减少对比时间；
   
2. 占用更少的内存；
   
3. 支持复杂比对模式：
     （a）剪接比对：RNA-seq数据需检测外显子连接，索引会预存剪切位点信息（HISAT2的snp_tran索引）；
     （b）突变容忍：包含SNP/突变的索引（如genome_snp_tran）可提高多态性样本的比对率。

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
fpkm <- read.table("/mnt/alamo01/users/chenyun730/program/test/quantify/all_samples_FPKM.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE) #运行后由于有重复的gene_id而报错
fpkm <- read.table("/mnt/alamo01/users/chenyun730/program/test/quantify/all_samples_FPKM.txt",header=TRUE, sep="\t", check.names=FALSE)
library(dplyr)#取均值合并
sum(duplicated(fpkm[, 1]))
fpkm_unique <- fpkm %>%
group_by(gene_id = fpkm[,1]) %>%
summarise(across(.cols = everything(), .fns = mean))
fpkm_matrix <- as.data.frame(fpkm_unique)
rownames(fpkm_matrix) <- fpkm_matrix$gene_id
fpkm_matrix$gene_id <- NULL

#再检查一下
 sum(duplicated(rownames(fpkm_matrix))) #此时输出为0，说明已经没有重复的gene_id

#取log2对数
fpkm_log2 <- log2(fpkm_matrix + 1)
```

2. 热图
```
/mnt/alamo01/users/chenyun730/program/test/results/
pheatmap(fpkm_log2[top500, ],
         scale="row",
         show_rownames=FALSE,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         main="Top 500 Most Variable Genes")
dev.off()

# 修改文件名
expr <- read.csv("your_expression_matrix.csv", row.names = 1)  # or .tsv with read.delim   # 读取表达矩阵
srr_ids <- c("SRR27961778", "SRR27961779", "SRR27961780", "SRR27961787", "SRR27961788", "SRR27961789") #指定你想提取的样本
new_names <- c("Mock-1", "Mock-2", "Mock-3", "SARS-2-1", "SARS-2-2", "SARS-2-3") #对应的目标命名
expr_subset <- expr[, srr_ids] # 提取感兴趣的6个样本
colnames(expr_subset) <- new_names  #重命名列
write.csv(expr_subset, file = "subset_6samples_renamed.csv") #保存结果

```

3. PCA图
```
fpkm_log2_filtered <- fpkm_log2[apply(fpkm_log2, 1, function(x) sd(x) != 0), ] #过滤掉NA
cat("Retained genes for PCA:", nrow(fpkm_log2_filtered), "\n")
fpkm_log2[is.na(fpkm_log2)] <- 0
pca_input <- t(fpkm_log2)
pca <- prcomp(pca_input, scale.=TRUE)
pca_df <- data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
library(ggplot2)
pdf("PCA_plot.pdf", width=8, height=6)
ggplot(pca_df, aes(PC1, PC2, label=Sample)) +
  geom_point(size=4, color="steelblue") +
  geom_text(vjust=-1) +
  theme_minimal() +
  labs(title="PCA of Samples", x="PC1", y="PC2")
dev.off()
```
下载原作者的矩阵进行比对

1.自己的counts矩阵先绘图检查异质性
```
micromamba install -c conda-forge -c bioconda bioconductor-geoquery
# 检查重复的 gene_id
sum(duplicated(rownames(counts)))  # 返回重复基因数量，返回0
sum(is.na(counts)))
counts[is.na(counts)] <- 0# 替换 NA 为 0（或根据需求过滤低表达基因）
keep <- rowSums(counts > 1) >= 3
counts <- counts[keep, ]# 保留至少在 3 个样本中 count > 1 的基因
colData <- data.frame(
  group = factor(rep(c("Mock", "SARS-2"), each = 3)),
  row.names = colnames(counts)
) # 定义样本分组（Mock vs SARS-2）

# 构建 DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ group
)
# 转换数据（blind=TRUE 适用于无监督分析）
rld <- rlog(dds, blind = TRUE)  #rlog 更适合小数据集（样本数 < 30）
rld_matrix <- assay(rld)  # 提取转换后的矩阵
#修改命名里的“-”
colData$group <- gsub("-", "_", colData$group)
ann_colors <- list(Group = c(Mock = "blue", SARS_2 = "red")) #热图的步骤
#绘制 PCA 图（分组着色）
pca_data <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

morandi_blue <- "#89CFF0" 
morandi_orange <- "#FFB347"
 ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Mock" = "#89CFF0", "SARS-2" = "#FFB347")) +
  xlab(paste0("PC1 (", percent_var[1], "%)")) +
  ylab(paste0("PC2 (", percent_var[2], "%)")) +
  ggtitle("PCA Plot of my_matrix") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  geom_text(aes(label = name), vjust = 1.5, size = 4)
ggsave("/mnt/alamo01/users/chenyun730/program/test/results/PCA_plot.png", width = 8, height = 6, dpi = 300)
#绘制热图
library(pheatmap)
top_genes <- head(order(rowVars(rld_matrix), decreasing = TRUE), 500)
annotation_col <- data.frame(
  Group = colData$group,
  row.names = colnames(counts))
ann_colors <- list(
  Group = c(Mock = "#89CFF0", SARS_2 = "#FFB347"))
pdf("/mnt/alamo01/users/chenyun730/program/test/results/heatmap_plot.pdf",width = 10, height = 12, pointsize = 15 )
pheatmap(
  rld_matrix[top_genes, ],
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Top 500 Genes by Variance (rlog-transformed)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
dev.off()
```
2.下载GEO原作者的矩阵进行分析（GSE255647_all.counts.tsv.gz）
```
# 读取 count 矩阵（gz 格式的 TSV 文件）
counts <- read.table("GSE255647_all.counts.tsv.gz", 
                     header = TRUE, 
                     sep = "\t", 
                     row.names = 1, 
                     check.names = FALSE)
head(counts[, 1:6])
str(counts)
# 提取目标样本
selected_samples <- c("SRR27961778", "SRR27961779", "SRR27961780", 
                      "SRR27961787", "SRR27961788", "SRR27961789")

counts_subset <- counts[, selected_samples]

# 重命名为更易懂的格式
colnames(counts_subset) <- c("Mock-1", "Mock-2", "Mock-3", 
                             "SARS-2-1", "SARS-2-2", "SARS-2-3")
# 选择目标样本
target_samples <- c("MOCK_12hpi_rep1", "MOCK_12hpi_rep2", "MOCK_12hpi_rep3",
                    "SARS-2_12hpi_rep1", "SARS-2_12hpi_rep2", "SARS-2_12hpi_rep3")

counts_subset <- counts[, target_samples]

# 重命名列
colnames(counts_subset) <- c("Mock-1", "Mock-2", "Mock-3",
                             "SARS-2-1", "SARS-2-2", "SARS-2-3")
library(DESeq2)

# 构建样本信息
group <- factor(c("Mock", "Mock", "Mock", "SARS-2", "SARS-2", "SARS-2"))
col_data <- data.frame(row.names = colnames(counts_subset),
                       condition = group)

# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = col_data,
                              design = ~ condition)

# 过滤低表达基因（推荐）
dds <- dds[rowSums(counts(dds)) > 10, ]

# 运行标准 DESeq2 流程
dds <- DESeq(dds)

# 绘制热图+PCA图
rld <- rlog(dds, blind = TRUE) #另一个 vsd <- vst(dds, blind = TRUE)
png("/mnt/alamo01/users/chenyun730/program/test/results/PCA_plot.png", 
    width = 800, height = 800)
plotPCA(rld, intgroup = "condition")
dev.off()

library(pheatmap)
library(matrixStats)
top_var_genes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 500) #rld改为vst
mat <- assay(rld)[top_var_genes, ]
mat <- t(scale(t(mat)))
pdf("/mnt/alamo01/users/chenyun730/program/test/results/Heatmap_geo.pdf", 
    width = 8, height = 10)
pheatmap(mat, 
         annotation_col = col_data,
         show_colnames = TRUE,
         show_rownames = FALSE,
         fontsize_row = 6，
          main = "Top 500 Genes by Variance of geo_matrix", )
dev.off()

```

3.将12个样本进行合并分析
**发现命名问题并矫正——symbol**
```
count_matrix <- read.csv("/mnt/alamo01/users/chenyun730/program/test/my_matrix/count_matrix.csv", row.names = NULL, check.names = FALSE)
count_matrix$symbol <- ifelse(grepl("\\|", count_matrix$gene_id),
                              sub(".*\\|", "", count_matrix$gene_id),
                              count_matrix$gene_id)
 library(dplyr)
 symbol_matrix <- count_matrix %>%
  select(-gene_id) %>%
  group_by(symbol) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()
symbol_matrix <- symbol_matrix %>%
  relocate(symbol, .before = everything())
 write.csv(symbol_matrix, "/mnt/alamo01/users/chenyun730/program/test/my_matrix/symbol_count_matrix.csv", row.names = FALSE)
```
 获得原作者的矩阵count_matrix_geo.csv
 ```
 counts_all <- read.table("/mnt/alamo01/users/chenyun730/program/test/geo_matrix/GSE255647_all.counts.tsv.gz",header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
selected_samples <- c("MOCK_24hpi_rep1", 
                      "MOCK_24hpi_rep2", 
                      "MOCK_24hpi_rep3", 
                      "SARS-2_24hpi_rep1", 
                      "SARS-2_24hpi_rep2", 
                      "SARS-2_24hpi_rep3")
count_matrix_geo <- counts_all[, selected_samples]
count_matrix_geo <- cbind(gene_id = rownames(count_matrix_geo), count_matrix_geo)
colnames(count_matrix_geo)[2:7] <- c("Mock-1", "Mock-2", "Mock-3", 
                                     "SARS-2-1", "SARS-2-2", "SARS-2-3")
write.csv(count_matrix_geo,
          file = "/mnt/alamo01/users/chenyun730/program/test/geo_matrix/count_matrix_geo.csv",
          row.names = FALSE,
          quote = FALSE)
```
 构建合并矩阵（merged_counts)
```
df1 <- read.csv("/mnt/alamo01/users/chenyun730/program/test/my_matrix/symbol_count_matrix.csv", check.names = FALSE)
df2 <- read.csv("/mnt/alamo01/users/chenyun730/program/test/geo_matrix/symbol_count_matrix_geo.csv", check.names = FALSE)
colnames(df2)[2:ncol(df2)] <- paste0("geo.", colnames(df2)[2:ncol(df2)])
 merged_df <- full_join(df1, df2, by = "symbol")
 merged_df[is.na(merged_df)] <- 0
 write.csv(merged_df, "/mnt/alamo01/users/chenyun730/program/test/results/merged_count_matrix.csv", row.names = FALSE)
```
PCA图
```
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(FactoMineR)
library(factoextra)
data <- read.csv("/mnt/alamo01/users/chenyun730/program/test/results/merged_count_matrix.csv", check.names = FALSE)
rownames(data) <- data$symbol
data <- data[, -1]
log_data <- log2(data + 1)
log_data_filtered <- log_data[apply(log_data, 1, function(x) sd(x) > 0), ]
pca_data <- t(log_data_filtered)
pca_result <- prcomp(pca_data, scale. = TRUE)
 group <- rownames(pca_data)
group <- ifelse(grepl("^geo\\.Mock", group), "geo.Mock",
         ifelse(grepl("^geo\\.SARS2", group), "geo.SARS2",
         ifelse(grepl("^Mock", group), "Mock", "SARS2")))
group <- factor(group, levels = c("Mock", "SARS2", "geo.Mock", "geo.SARS2"))
 pca_plot <- fviz_pca_ind(pca_result,
                         geom.ind = "point",
                         pointshape = 21,
                         col.ind = group,
                         palette = c("#1f77b4", "#d62728", "#2ca02c", "#9467bd"),
                         addEllipses = TRUE,
                         legend.title = "Group") +
  theme_minimal() +
  ggtitle("PCA of RNA-seq Samples")
ggsave("PCA_plot.png", plot = pca_plot, width = 6, height = 5, dpi = 300)
```
热图
```
 library(RColorBrewer)
> data <- read.csv("/mnt/alamo01/users/chenyun730/program/test/results/merged_count_matrix.csv", header = TRUE, row.names = 1)
> count_data <- data[, 1:12]
sum(is.na(scaled_data)) #去除na和全是0的
 scaled_data[is.na(scaled_data)] <- 0
scaled_data[is.infinite(scaled_data)] <- 0
nonzero_rows <- rowSums(count_data) > 0
filtered_data <- scaled_data[nonzero_rows, ]
sum(is.infinite(scaled_data))
 log_data <- log2(count_data + 1)
> scaled_data <- t(scale(t(log_data)))
> sample_groups <- c(rep("Mock", 3), rep("SARS2", 3), rep("geo.Mock", 3), rep("geo.SARS2", 3))
> annotation_col <- data.frame(Group = factor(sample_groups))
rownames(annotation_col) <- colnames(scaled_data)
> group_colors <- list(
  Group = c(
    Mock = "#1F77B4",
    SARS2 = "#FF7F0E",
    geo.Mock = "#2CA02C",
    geo.SARS2 = "#D62728"
  )
)
pheatmap(
  scaled_data,
  main = "Heatmap of merged_matrix",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = group_colors,
  scale = "row",
  border_color = NA,
  fontsize_col = 12,
  angle_col = 45,
  filename = "/mnt/alamo01/users/chenyun730/program/test/results/heatmap_merged_matrix.pdf",
  width = 10,
  height = 14
)
```
**将outcome和geo进行分组进行差异分析**
```
 library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
count_data <- read.csv("/mnt/alamo01/users/chenyun730/program/test/results/merged_count_matrix.csv", row.names = 1)
 condition <- factor(c(rep("outcome", 6), rep("geo", 6)))
 coldata <- data.frame(condition)
rownames(coldata) <- colnames(count_data)
 stopifnot(identical(rownames(coldata), colnames(count_data)))
 dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)
 dds <- dds[rowSums(counts(dds)) > 10, ]
 dds <- DESeq(dds)
 res <- results(dds, contrast = c("condition", "geo", "outcome"))
 res <- res[order(res$padj), ]
 write.csv(as.data.frame(res), file = "DESeq2_diff_expr_results.csv")
```

**思考题**：

如果PCA中样本未能按组分离，原因可能是什么？

数据为标准化（这里做了）：log2 转换以及scale.=TRUE；太多零表达或低变异的基因，可能掩盖了真实的组间差异；组间差异太小导致的生物差异不明显等原因。

### 模块5：差异表达分析与功能富集
**目标**：差异表达分析与功能富集，理解生物学意义。

**任务**：

DESeq2筛选差异表达基因(outcome与geo）
```
# outcome_matrix
count_data <- read.csv("/mnt/alamo01/users/chenyun730/program/test/my_matrix/symbol_count_matrix.csv",row.names = 1)
 condition <- factor(c(rep("mock", 3), rep("sars2", 3)))
count_data[is.na(count_data)] <- 0
 coldata <- data.frame(condition)
rownames(coldata) <- colnames(count_data)
 stopifnot(identical(rownames(coldata), colnames(count_data)))
 dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)
 dds <- dds[rowSums(counts(dds)) > 10, ]
 dds <- DESeq(dds)
 res <- results(dds, contrast = c("condition","mock","sars2"))
 res <- res[order(res$padj), ]
 write.csv(as.data.frame(res), file = "/mnt/alamo01/users/chenyun730/program/test/my_matrix/outcome_DESeq2_diff_expr_results.csv")

# geo_matrix
count_data <- read.csv("/mnt/alamo01/users/chenyun730/program/test/geo_matrix/count_matrix_geo.csv",row.names = 1)
 condition <- factor(c(rep("mock", 3), rep("sars2", 3)))
count_data[is.na(count_data)] <- 0
 coldata <- data.frame(condition)
rownames(coldata) <- colnames(count_data)
 stopifnot(identical(rownames(coldata), colnames(count_data)))
 dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)
 dds <- dds[rowSums(counts(dds)) > 10, ]
 dds <- DESeq(dds)
 res <- results(dds, contrast = c("condition","mock","sars2"))
 res <- res[order(res$padj), ]
 write.csv(as.data.frame(res), file = "/mnt/alamo01/users/chenyun730/program/test/geo_matrix/geo_DESeq2_diff_expr_results.csv")

```
**clusterProfiler进行GO/KEGG富集分析**

```
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(DOSE)

# 筛选显著差异基因
deg <- as.data.frame(res)
deg <- na.omit(deg)
deg_filtered <- deg[deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, ]

# 提取 gene symbol 或 ENSEMBL id
gene_symbols <- rownames(deg_filtered)
# symbol或 "ENSEMBL"，根据gene ID 类型 转 Entrez ID
gene_entrez <- bitr(gene_symbols,
                    fromType = "SYMBOL", 
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# 提取 ENTREZID 向量
entrez_ids <- gene_entrez$ENTREZID
ego <- enrichGO(gene         = entrez_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "ALL",       # BP/CC/MF/ALL
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

# 保存结果
write.csv(as.data.frame(ego), file = "GO_enrichment_results.csv")

# 可视化（柱状图 / 点图）
png("GO_barplot_geo.png", width = 800, height = 600)
barplot(ego, showCategory = 20, title = "GO Enrichment_geo")
dev.off()
png("GO_dotplot_geo.png", width = 800, height = 600)
dotplot(ego, showCategory = 20, title = "GO Enrichment_geo")
dev.off()
```
**GSEA**
```
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)  # 如果分析 KEGG 可以用到
library(tibble)
library(enrichplot)
# 已运行DEseq2得到res
# 去掉 NA
res <- na.omit(res)

# 提取 SYMBOL 和 log2FoldChange，转为 Entrez ID
res$symbol <- rownames(res)
gene_info <- bitr(res$symbol,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# 合并
res_merged <- merge(res, gene_info, by.x = "symbol", by.y = "SYMBOL")

# 构造排序向量：ENTREZID 为 names，log2FC 为值
gene_list <- res_merged$log2FoldChange
names(gene_list) <- res_merged$ENTREZID

# 按 log2FC 降序排序（非常重要！）
gene_list <- sort(gene_list, decreasing = TRUE)

# 运行 GSEA（以 GO 为例）
ego_gsea <- gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENTREZID",
                  nPerm = 1000,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  verbose = TRUE)

# 保存富集图
png("GSEA_dotplot.png", width = 1000, height = 800)
dotplot(ego_gsea, showCategory = 20, title = "GSEA-GO BP")
dev.off()

# 也可以画 GSEA 单条路径图（比如 top1）
png("GSEA_pathway1.png", width = 1000, height = 600)
gseaplot2(ego_gsea, geneSetID = ego_gsea$ID[1], title = ego_gsea$Description[1])
dev.off()

```
输出：差异表达结果、富集结果

思考题：

为什么要进行多重检验校正？

### 模块六 更换比对方法验证跟geo的差异
**1. 下载STAR文件并构建参考基因组，比对定量得到count矩阵star_count_matrix.txt**
```
#新建一个test_star文件储存该模块的文件
#download.sh用于下载STAR的GRCh38（2.7.0）及注释文件
vim download.sh
#! /bin/bash
#micromamba activate R441
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz
gunzip gencode.v42.primary_assembly.annotation.gtf.gz

qsub -cwd -V -l cpu=64:mem=64G -q fast -N homodata_download /mnt/alamo01/users/chenyun730/program/test_star/scripts/download.sh

#star_index.sh用于构建索引
 micromamba create -n star_env -c bioconda -c conda-forge star=2.7.10a
 micromamba activate star_env
 STAR --version
vim star_index.sh
#! /bin/bash
#micromamba activate R441
module load STAR/2.7.10a
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v42.primary_assembly.annotation.gtf \
     --sjdbOverhang 100

 qsub -cwd -V -l cpu=64:mem=64G -q fast -N star_index /mnt/alamo01/users/chenyun730/program/test_star/scripts/star_index.sh

#进行比对
vim star_align.sh
#! /bin/bash
#micromamba activate R441
module load STAR/2.7.10a
INPUT_DIR=/mnt/alamo01/users/chenyun730/program/test_hisat2/clean_data
OUTPUT_DIR=/mnt/alamo01/users/chenyun730/program/test_star/alignment
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"
STAR --runThreadN 8 \
     --genomeDir /mnt/alamo01/users/chenyun730/program/test_star/homodata/star_index \
     --readFilesIn ${INPUT_DIR}/${SAMPLE}_cleaned_1.fp.gz ${INPUT_DIR}/${SAMPLE}_cleaned_2.fp.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --twopassMode Basic
done

# 使用FeatureCount定量
vim  featurecounts.sh
#! /bin/bash
#micromamba activate R441
BAM_DIR=/mnt/alamo01/users/chenyun730/program/test_star/alignment
GTF_FILE=/mnt/alamo01/users/chenyun730/program/test_star/homodata/gencode.v42.primary_assembly.annotation.gtf
OUTPUT_DIR=/mnt/alamo01/users/chenyun730/program/test_star/quantification
COUNT_FILE=${OUTPUT_DIR}/gene_counts.txt
SUMMARY_FILE=${OUTPUT_DIR}/summary.txt
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
BAM_FILES=()
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILES+=("${BAM_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam")
done
featureCounts \
    -T 8 \
    -a ${GTF_FILE} \
    -o ${COUNT_FILE} \
    -g gene_name \
    -p \
     --countReadPairs \
     -s 1 \
     --verbose \
    ${BAM_FILES[@]} > ${SUMMARY_FILE} 2>&1
cut -f1,7- ${COUNT_FILE} | grep -v '^#' | \
awk 'BEGIN{OFS="\t"}
     NR==1 {print "gene_id", "mock_1", "mock_2", "mock_3", "sars2_1", "sars2_2", "sars2_3"; next}
     {print $0}' > ${OUTPUT_DIR}/star_count_matrix.csv
```

**2.构建合并矩阵star_merged.csv**
```
#先读入生成的star_count_matrix.txt转成symbol再以.csv输出为symbol_star_count.csv
 df<- read.delim("/mnt/alamo01/users/chenyun730/program/test_star/quantification/star_count_matrix.txt", row.names = 1, check.names = FALSE)
 df$gene_id <- rownames(df)
rownames(df) <- NULL
 is_ensembl <- grepl("^ENSG", df$gene_id)
 ensembl_ids <- as.character(df$gene_id[is_ensembl])
df$symbol[is_ensembl] <- id_map$SYMBOL[match(df$gene_id[is_ensembl], id_map$ENSEMBL)]
df$symbol <- df$gene_id
df$symbol[is_ensembl] <- id_map$SYMBOL[match(df$gene_id[is_ensembl], id_map$ENSEMBL)]
library(dplyr)
df_clean <- df %>%
filter(!is.na(symbol) & symbol != "") %>%
select(-gene_id) %>%
group_by(symbol) %>%
summarise(across(everything(), sum)) %>%
ungroup() %>%
relocate(symbol)
write.csv(df_clean, "/mnt/alamo01/users/chenyun730/program/test_star/my_matrix/symbol_star_count.csv", row.names = FALSE)

# 合并两个symbol矩阵为star_merged.csv
df1<- read.csv("/mnt/alamo01/users/chenyun730/program/test_star/my_matrix/symbol_star_count.csv",row.names = 1, check.names = FALSE)
df2<- read.csv("/mnt/alamo01/users/chenyun730/program/test_hisat2/geo_matrix/symbol_count_matrix_geo.csv",row.names = 1, check.names = FALSE)
merged_df <- full_join(df1, df2, by = "symbol")
merged_df[is.na(merged_df)] <- 0
write.csv(merged_df, "/mnt/alamo01/users/chenyun730/program/test_star/results/star_merged.csv", row.names = FALSE)
```
PCA
```
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(FactoMineR)
library(factoextra)
data <- read.csv("/mnt/alamo01/users/chenyun730/program/test_star/my_matrix/symbol_star_count.csv", check.names = FALSE)
rownames(data) <- data$symbol
data <- data[, -1]
log_data <- log2(data + 1)
log_data_filtered <- log_data[apply(log_data, 1, function(x) sd(x) > 0), ]
pca_data <- t(log_data_filtered)
pca_result <- prcomp(pca_data, scale. = TRUE)
 group <- rownames(pca_data)
group = factor(rep(c("mock", "sars2"), each = 3))
 pca_plot <- fviz_pca_ind(pca_result,
                         geom.ind = "point",
                         pointshape = 21,
                         col.ind = group,
                         palette = c("#1f77b4", "#d62728", "#2ca02c", "#9467bd"),
                         addEllipses = TRUE,
                         legend.title = "Group") +
  theme_minimal() +
  ggtitle("PCA of RNA-seq Samples")
ggsave("PCA_plot.png", plot = pca_plot, width = 6, height = 5, dpi = 300)
```

**用Stingtie也一样处理STAR比对后的bam文件，然后得到stingtie_merged.csv数据**
```
#生成.gtf文件
vim stringtie_star.sh
#! /bin/bash
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample ${SAMPLE}..."
     stringtie /mnt/alamo01/users/chenyun730/program/test_star/alignment/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
        -G /mnt/alamo01/users/chenyun730/program/test_star/homodata/gencode.v42.primary_assembly.annotation.gtf \
        -o ./${SAMPLE}/${SAMPLE}.gtf  \
         -e \
         -B \
         -p 8
done
python3 /mnt/alamo01/users/chenyun730/program/test_hisat2/scripts/prepDE.py3 .

#生成合并矩阵
 df1<- read.csv("symbol_gene_count.csv", check.names = FALSE)  
 df2<- read.csv("/mnt/alamo01/users/chenyun730/program/test_hisat2/geo_matrix/symbol_count_matrix_geo.csv", check.names = FALSE)
 merged_df <- full_join(df1, df2, by = "symbol")                                 merged_df[is.na(merged_df)] <- 0
write.csv(merged_df, "/mnt/alamo01/users/chenyun730/program/test_star/results/stringtie_star_merged.csv", row.names = FALSE)

#PCA
data<- read.csv("/mnt/alamo01/users/chenyun730/program/test_star/results/stringtie_gene_star.csv",row.names = 1, check.names = FALSE)
rownames(data) <- data$symbol
data <- data[, -1]
log_data <- log2(data + 1)
log_data_filtered <- log_data[apply(log_data, 1, function(x) sd(x) > 0), ]
pca_data <- t(log_data_filtered)
pca_result <- prcomp(pca_data, scale. = TRUE)
 group <- rownames(pca_data)
group = factor(rep(c("mock", "sars2", each = 3))
 pca_plot <- fviz_pca_ind(pca_result,
                         geom.ind = "point",
                         pointshape = 21,
                         col.ind = group,
                         palette = c("#1f77b4", "#d62728", "#2ca02c", "#9467bd"),
                         addEllipses = TRUE,
                         legend.title = "Group") +
  theme_minimal() +
  ggtitle("PCA of RNA-seq Samples")
ggsave("/mnt/alamo01/users/chenyun730/program/test_star/results/PCA_stringtie_star.png", plot = pca_plot, width = 6, height = 5, dpi = 300)
Error: unexpected symbol in:
"group = factor(rep(c("mock", "sars2","geo_mock","geo_sars2", each = 3))
 pca_plot"

```
**检查多映射reads比例（SAM/BAM文件中的XS标签），过高（>10%）可能提示同源/重复序列问题**
```
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961778/SRR27961778_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 1415394
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961779/SRR27961779_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18106826
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961780/SRR27961780_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 19704911
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961787/SRR27961787_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18680847
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961788/SRR27961788_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18875024
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961789/SRR27961789_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18252563
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$
```
