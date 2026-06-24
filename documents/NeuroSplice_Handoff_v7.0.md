## 🇬🇧 English Version: NeuroSplice Team Handoff Document

# 🧬 NeuroSplice – Team Technical Handoff Document

**Date:** May 2026 **Domain:** Computational Neuroscience & Alternative Splicing **Primary Compute Environment:** Yale HPC (Bouchet Cluster / Gibbs Farnam, Open OnDemand)

## 1\. Project Summary

**Goal:** Investigate the transcriptomic footprint (specifically Alternative Splicing \[AS\] and Differential Gene Expression \[DEG\]) of an ADHD mouse model. **Core Finding:** The knockout of a major ADHD risk gene causes minimal transcriptional changes (DGE) but **widespread splicing dysregulation (306-341 events per region)** across multiple brain regions. **Current Phase:** Translational Evidence Consolidation & Figure Generation.

## 2\. Dataset Metadata

* **Dataset Accession:** GSE117357 (NCBI GEO)  
* **Model:** *Adgrl3* Knockout (KO) vs. Wild Type (WT)  
* **Background:** C57BL/6J  
* **Tissue/Regions:** Hippocampus, Prefrontal Cortex (PFC), Striatum  
* **Sample Size:** $n=60$ total (10 WT vs 10 KO per region)  
* **Age/Sex:** Adult (8-12 weeks), All Male  
* **Library Prep:** Single-end, strand-specific RNA-seq  
* **Read Length:** 76bp

## 3\. Pipeline Specifications & Parameters

**Alignment & Gene Counting**

* **Aligner:** STAR v2.7.10b  
* **Gene Quantification:** `featureCounts`  
  * **Command:** `featureCounts -a /path/to/annotations.gtf -o counts.txt -s 2 alignments/*.bam`  
  * *(Note: `-s 2` confirmed reverse-stranded via metadata audit)*  
* **Reference Genome:** mm10  
* **Annotation:** Ensembl 102 (GTF)

**Alternative Splicing (AS) Analysis**

* **Tool:** rMATS v4.3.0  
* **Key Command Parameters:** `--readLength 76 --t single --novelSS 1 --allow-clipping --libType fr-firststrand`  
* **Filtering Thresholds:**  
  * Significance: FDR $\< 0.05$  
  * Effect Size (Category A/Systemic): $|\\Delta\\Psi| \\ge 0.1$  
  * Effect Size (Category B/Tissue-Specific): $|\\Delta\\Psi| \\ge 0.3$

**Differential Gene Expression (Complete)**

* **Tool:** DESeq2 v1.38.3  
* **Filtering Thresholds:** padj $\< 0.05$, $|log\_2FC| \> 0.5$  
* **Design:** `~genotype` (Run independently per brain region; cohort confirmed to be all-male)

## 4\. Key Results & Candidate Categorization

*Total unique core genes under investigation: 9*

### 4.1 AS Analysis Categories

**Category A: Pan-Tissue Splicing Instability (Gene-Level)** *Definition: Splicing disrupted in all 3 regions, validating a systemic effect across independent biological replicates.*

* **Genes (3):** *Pts*, *Lrp8*, *Myo9b*

**Category B1: High-Confidence Region-Specific Drivers** *Definition: $|\\Delta\\Psi| \\ge 0.3$ (single tissue) OR $|\\Delta\\Psi| \\ge 0.2$ in $\\ge 2$ tissues AND strong functional prior (synaptic/neuronal survival genes).*

* **Genes (3):** *Neil2* (PFC), *Unc13b* (Stri/Hippo), *Tox3* (PFC)  
* *(Note: The Pts-Striatum event is counted under Category A for gene-level accounting; its B1 classification is functional, not additive to the gene count).*

**Category C: Contextual/Regulatory Modulators (Ambiguous)** *Definition: Significant but features direction-flipping across tissues.*

* **Genes (2):** *Arid5a*, *Bcl2l11*

**Category D: Technical Controls** *Definition: Non-coding RNA, used purely as pipeline sanity checks.*

* **Genes (1):** *Gm10419* (lincRNA)

### 4.2 DEG Analysis Status

* \[COMPLETE\] Hippocampus: 19 significant DEGs
* \[COMPLETE\] PFC: 69 significant DEGs
* \[COMPLETE\] Striatum: 22 significant DEGs

### 4.3 Cross-Species Validation Protocols

**1\. Virtual PCR & Conservation Check**

* **Protocol:**  
  1. Extract coordinates from rMATS (e.g., *Pts* RI target coordinates: `chr9:50524722-50526949` mm10). *(Note: extracted from RI.MATS.JCEC.txt, confirm before use).*  
  2. **Human Ortholog Lookup:** Map mouse gene symbols to human orthologs via Ensembl REST API (cross-species coordinate mapping via REST API is unsupported).  
  3. **Transcript Verification:** Verify human orthologs have known Ensembl exon structures to exclude pseudogenes.  
  4. Conservation threshold: PhyloP100 score $\> 2.0$ or PhastCons $\> 0.8$ (Requires HPC UCSC liftOver binary).  
  5. GWAS Check: Open Targets Genetics query for psychiatric associations.
  * **Results:** 6/6 candidate genes have confirmed human protein-coding orthologs. 6/6 genes have psychiatric GWAS associations, including a direct ADHD hit for *NEIL2*.

**2\. Human Digital Phenotype Intersection**

* **Reference:** Liu & Borsari et al., *Cell* 2025 (DOI: 10.1016/j.cell.2025.01.014). Identifies 37 GWAS-associated psychiatric/digital-phenotype genes.  
* **Intersection Strategy:**  
  * *Primary:* Intersect regional **AS Genes** ($FDR \< 0.05$) vs. 37-gene list.  
  * *Secondary:* Intersect regional **DEG lists** ($padj \< 0.05$) vs. 37-gene list.  
* **Statistical Test:** Fisher's Exact Test. **Background \= all expressed genes, strictly defined as `baseMean ≥ 10` in DESeq2.**  
* **Significance:** Apply Bonferroni correction across the 6 independent tests (3 regions $\\times$ 2 layers). Target threshold: $p \< 0.0083$.
* **Results:** 0/6 tests significant. Sole overlap was *Dlg4* in Hippocampus AS (not enriched, OR < 1). Expected result as GWAS and KO models probe different biological tiers.

## 5\. Candidate Gene Summary Table

| Gene | Category | Highest $|\\Delta\\Psi|$ | Region | Event Type | Biological Notes / Mechanism | | :--- | :--- | :--- | :--- | :--- | :--- | | **Neil2** | B1 | \+0.460 | PFC | SE | Oxidative DNA repair $\\rightarrow$ neuronal genome integrity $\\rightarrow$ expressed in PFC interneurons. ADHD link is indirect, but highest effect size warrants investigation. | | **Pts** | A & B1 | \+0.301 | Striatum | RI | Predicted to disrupt BH4 synthesis, a rate-limiting step in dopamine biosynthesis. Functional consequence requires validation. | | **Tox3** | B1 | \+0.336 | PFC | A3SS | Neuronal survival $\\rightarrow$ calcium-dependent transcription $\\rightarrow$ expressed in PFC pyramidal neurons. ADHD link indirect; large PFC effect warrants investigation. | | **Unc13b** | B1 | \-0.242 | Hippo/Stri | SE | Synaptic vesicle priming. *(Note: Effect size just above category threshold; inclusion justified by strong prior as a core synaptic priming factor).* | | **Myo9b** | A | \-0.217 | PFC | SE | RhoGAP $\\rightarrow$ dendritic spine dynamics $\\rightarrow$ spine morphology defects reported in ADHD models. | | **Lrp8** | A | \-0.163 | Striatum | SE | Reelin receptor $\\rightarrow$ synaptic plasticity $\\rightarrow$ ADHD-associated pathway (*Reln* also dysregulated in Striatum). |

## 6\. Failed Approaches & Negative Results (Do Not Repeat)

1. **30v30 Pooled Analysis:** NOT RECOMMENDED. Pooling hippocampus, PFC, and striatum as replicates violates biological independence and artificially inflates statistical power. Use per-region 10v10 comparisons to demonstrate systemic effects.  
2. **Initial HISAT2 Alignment:** Initial alignment used HISAT2 v2.2.2 (without `--dta` flag). This was superseded by STAR for full rMATS compatibility. HISAT2 BAMs remain on the HPC (same paths as STAR BAMs), are marked as deprecated, and should NOT be used for any downstream analysis.  
3. **DaPars Analysis:** Failed due to insufficient 3' coverage for APA (Alternative Polyadenylation) calling.  
4. **CAGE Intersection:** Failed because no public mouse brain CAGE data is available to validate Transcription Start Sites (TSS).  
5. **LincRNA Excluded:** *Gm10419* showed consistent splicing but was excluded from mechanistic claims as it lacks translatability to human protein-coding paradigms.  
6. **Directionality Flips:** *Arid5a* and *Bcl2l11* demonstrated tissue-specific direction flipping. Demoted from "Core" to Category C; interpreted as regulatory complexity rather than simple loss-of-function.

## 7\. File System & Data Locations

**Yale HPC (Bouchet / Gibbs Cluster):**

/gpfs/gibbs/project/neurosplice\_gse117357/

├── raw\_data/                  (FASTQ files)

├── alignments/                (STAR BAMs)

└── rmats\_out/                 (rMATS v4.3.0 raw output)

/palmer/scratch/\[PI\_NAME\]/neurosplice\_backup/  \<-- Target backup location for external drive

**Local Drive (Mac Desktop):**

~/Desktop/WXY/NeuroSplice/results/

├── rmats\_gse117357\_analysis/  (AS results, filtered CSVs)

├── deseq2/                    (DEG output \- COMPLETE)

├── qc/multiqc\_align.html      (QC report)

└── figures/                   (Current Sashimi plots, Volcano plots)

**Local/Shared Drive \[🚨 CRITICAL WARNING \- EXTERNAL DRIVE\]:**

/Volumes/Untitled/NeuroSplice/  \<-- Contains 60GB raw BAMs; MUST BE BACKED UP TO HPC

## 8\. Action Items & Next Steps

| Task | Description | Owner | Deadline | Status |
| :---- | :---- | :---- | :---- | :---- |
| **Data Security** | Back up external drive BAMs to HPC `palmer_scratch` immediately. | Drive Owner (Coord: Jason) | ASAP | 🔴 To Do |
| **Metadata Audit** | Verify library strandedness (`fr-firststrand`) via GEO/SRA raw metadata. | Jason | May 11 | ✅ Complete |
| **DEG Run** | Execute DESeq2 pipeline (`~genotype`) and `featureCounts` for 3 regions. | Jason | May 12 | ✅ Complete |
| **Fisher's Test** | Run overlap test using `baseMean ≥ 10` background and Bonferroni correction. | Jason | May 14 | ✅ Complete |
| **Virtual PCR** | Human ortholog + GWAS complete. PhyloP scores pending HPC liftOver. | Sarah | May 15 | 🟡 Partial |
| **Figure 1 Draft** | Generate Sashimi plots for Pts and Neil2. *(Inputs: `alignments/[sample].bam` \+ `rmats_out/[Event].MATS.JCEC.txt`)* | Sarah | May 18 | 🔴 To Do |
| **Figure 2 Draft** | Generate Sashimi plots/metrics for Unc13b and Tox3. *(Inputs: same as Fig 1\)* | Sarah | May 20 | 🔴 To Do |

---

---

### 🇨🇳 中文版本：NeuroSplice 团队技术交接文档 (v7.0)

# 🧬 NeuroSplice – 团队技术交接文档

**日期:** 2026年5月 **领域:** 计算神经科学与可变剪接 (Alternative Splicing, AS) **主要计算环境:** 耶鲁大学 HPC (Bouchet Cluster / Gibbs Farnam, Open OnDemand)

## 1\. 项目摘要

**目标:** 研究 ADHD 小鼠模型的转录组学特征（重点关注可变剪接 \[AS\] 和差异基因表达 \[DEG\]）。 **核心发现:** 该主要 ADHD 风险基因的敲除几乎没有引起转录水平的变化（DGE极少），但在多个脑区引发了**广泛的剪接失调（每个脑区 306-341 个事件）**。 **当前阶段:** 转化证据整理与作图阶段 (Translational Evidence Consolidation & Figure Generation)。

## 2\. 数据集元数据

* **数据集编号:** GSE117357 (NCBI GEO)  
* **小鼠模型:** *Adgrl3* 敲除 (KO) vs. 野生型 (WT)  
* **遗传背景:** C57BL/6J  
* **组织/脑区:** 海马体 (Hippocampus)、前额叶皮层 (PFC)、纹状体 (Striatum)  
* **样本量:** 总计 $n=60$（每个脑区 10 WT vs 10 KO）  
* **年龄/性别:** 成年（8-12周），全雄性 (All Male)  
* **建库策略:** 单端测序 (Single-end)，链特异性 RNA-seq  
* **读长:** 76bp

## 3\. 流程规范与参数

**比对与基因定量**

* **比对软件:** STAR v2.7.10b  
* **基因定量:** `featureCounts`  
  * **完整命令示例:** `featureCounts -a /path/to/annotations.gtf -o counts.txt -s 2 alignments/*.bam`  
  * *（注：通过元数据审计已确认测序方向为 reverse-stranded，因此使用 `-s 2`）*  
* **参考基因组:** mm10  
* **注释文件:** Ensembl 102 (GTF)

**可变剪接 (AS) 分析**

* **分析工具:** rMATS v4.3.0  
* **核心命令参数:** `--readLength 76 --t single --novelSS 1 --allow-clipping --libType fr-firststrand`  
* **过滤阈值:** \* 显著性: FDR $\< 0.05$  
  * 效应量 (类别 A/系统性): $|\\Delta\\Psi| \\ge 0.1$  
  * 效应量 (类别 B/组织特异性): $|\\Delta\\Psi| \\ge 0.3$

**差异基因表达 (DEG) 分析 (已完成)**

* **分析工具:** DESeq2 v1.38.3  
* **过滤阈值:** padj $\< 0.05$, $|log\_2FC| \> 0.5$  
* **模型设计:** `~genotype`（针对每个脑区独立运行；已确认为全雄性队列，因此移除 sex 协变量）

## 4\. 核心结果与候选基因分类

*总计研究的独特核心基因数量: 9个*

### 4.1 AS 分析类别

**类别 A: 泛组织剪接不稳定性 (基因水平)** *定义: 在所有3个脑区中均发现剪接异常，在独立的生物学重复中验证了系统性效应。*

* **基因 (3个):** *Pts*, *Lrp8*, *Myo9b*

**类别 B1: 高置信度区域特异性驱动因子** *定义: $|\\Delta\\Psi| \\ge 0.3$ (单一组织) 或 $|\\Delta\\Psi| \\ge 0.2$ (在 $\\ge 2$ 个组织中)，且具有强烈的生物学功能先验（如突触/神经元存活基因）。*

* **基因 (3个):** *Neil2* (PFC), *Unc13b* (海马/纹状体), *Tox3* (PFC)  
* *（注：Pts-纹状体事件在基因水平的统计中归入类别A；它在B1的分类是基于功能的，不增加总基因计数）。*

**类别 C: 语境/调控调节因子 (存在歧义)** *定义: 具有统计学显著性，但在不同组织间存在“方向翻转”现象。*

* **基因 (2个):** *Arid5a*, *Bcl2l11*

**类别 D: 技术质控对照** *定义: 非编码 RNA，纯粹用作分析流程的合理性检验。*

* **基因 (1个):** *Gm10419* (lincRNA)

### 4.2 DEG 分析状态

* \[已完成\] 海马体: 19 个显著差异基因 (DEGs)  
* \[已完成\] PFC: 69 个显著差异基因 (DEGs)  
* \[已完成\] 纹状体: 22 个显著差异基因 (DEGs)

### 4.3 跨物种验证方案

**1\. 虚拟 PCR 与保守性检查**

* **操作步骤:**  
  1. 从 rMATS 输出文件中提取坐标（例如，Pts RI 目标坐标: `chr9:50524722-50526949` mm10）。*（注：从 `RI.MATS.JCEC.txt` 中提取，使用前请确认）。*  
  2. **人类直系同源映射:** 通过 Ensembl REST API 将小鼠基因符号映射到人类直系同源基因（不支持跨物种坐标映射）。  
  3. **转录本核实:** 验证人类直系同源基因具有已知的 Ensembl 外显子结构，以排除假基因。  
  4. 保守性阈值: PhyloP100 得分 $\> 2.0$ 或 PhastCons $\> 0.8$（需要在 HPC 上运行 UCSC liftOver）。  
  5. GWAS 检查: 在 Open Targets Genetics 查询精神疾病相关性。
  * **结果:** 6/6 个候选基因具有确认的人类蛋白质编码直系同源物。6/6 个基因在 GWAS 中有精神疾病相关性，其中 *NEIL2* 直接与 ADHD 相关。

**2\. 人类数字表型 (Digital Phenotype) 交集分析**

* **参考文献:** Liu & Borsari et al., *Cell* 2025 (DOI: 10.1016/j.cell.2025.01.014)。鉴定了 37 个与 GWAS 相关的精神病学/数字表型基因。  
* **交集策略:** \* *主要测试:* 将各脑区的 **AS 基因** ($FDR \< 0.05$) 与 37 基因列表求交集。  
  * *次要测试:* 将各脑区的 **DEG 列表** ($padj \< 0.05$) 与 37 基因列表求交集。  
* **统计测试:** Fisher 精确检验。**背景集 \= 所有表达基因（严格定义为 DESeq2 中 `baseMean ≥ 10`）。**  
* **显著性标准:** 对 6 个独立测试（3 个脑区 $\\times$ 2 个分子层）应用 Bonferroni 校正。目标阈值: $p \< 0.0083$。
* **结果:** 0/6 项测试显著。唯一的重叠是海马体 AS 中的 *Dlg4*（未富集，OR < 1）。这是符合预期的，因为 GWAS 与敲除模型探究的是不同层级的生物学现象。

## 5\. 候选基因总结表

| 基因 | 类别 | 最高 $|\\Delta\\Psi|$ | 脑区 | 事件类型 | 生物学注释 / 潜在机制 | | :--- | :--- | :--- | :--- | :--- | :--- | | **Neil2** | B1 | \+0.460 | PFC | SE | 氧化 DNA 修复 $\\rightarrow$ 神经元基因组完整性 $\\rightarrow$ 在 PFC 中间神经元中表达。与 ADHD 的联系是间接的，但因其拥有数据集中最大的效应量，值得深入研究。 | | **Pts** | A & B1 | \+0.301 | 纹状体 | RI | 预测会破坏 BH4 的合成，这是多巴胺生物合成的限速步骤。其功能后果需要后续验证。 | | **Tox3** | B1 | \+0.336 | PFC | A3SS | 神经元存活 $\\rightarrow$ 钙依赖性转录 $\\rightarrow$ 在 PFC 锥体神经元中表达。与 ADHD 的联系是间接的，但因其在 PFC 有巨大的效应量，值得深入研究。 | | **Unc13b** | B1 | \-0.242 | 海马/纹状体 | SE | 突触小泡启动 (Synaptic vesicle priming)。*（注：其效应量刚好超过类别阈值；但由于其作为核心突触启动因子的强烈先验支持，被优先纳入）。* | | **Myo9b** | A | \-0.217 | PFC | SE | RhoGAP $\\rightarrow$ 树突棘动态 $\\rightarrow$ 在 ADHD 模型中已有树突棘形态缺陷的相关报道。 | | **Lrp8** | A | \-0.163 | 纹状体 | SE | Reelin 受体 $\\rightarrow$ 突触可塑性 $\\rightarrow$ ADHD 相关通路（*Reln* 基因在纹状体中也同样发生失调）。 |

## 6\. 失败的尝试与避坑指南（请勿重复）

1. **30v30 混合池化分析:** 强烈不推荐。将海马、PFC 和纹状体作为重复样本混合，违反了生物学独立性假设，并会人为夸大统计功效。必须使用各脑区独立的 10v10 比较来证明系统性效应。  
2. **早期的 HISAT2 比对:** 最初比对使用了 HISAT2 v2.2.2（未加 `--dta` 参数）。为了与 rMATS 完全兼容，该方案已被 STAR 取代。HPC 上仍残留有 HISAT2 的 BAM 文件（路径同 STAR BAMs），已标记为过时版本，请勿用于任何下游分析。  
3. **DaPars 分析:** 失败。由于 3' 端覆盖度不足，无法进行 APA（可变多聚腺苷酸化）调用。  
4. **CAGE 交集验证:** 失败。因为目前缺乏公开的小鼠大脑 CAGE 数据来验证转录起始位点 (TSS)。  
5. **LincRNA 排除:** *Gm10419* 虽然显示出一致的剪接变化，但已被排除在机制推演之外，因为它无法转化为人类蛋白质编码模型。  
6. **方向性翻转:** *Arid5a* 和 *Bcl2l11* 在不同组织间表现出剪接方向的翻转。已从“核心”降级为类别 C；这被解释为复杂的调控环境差异，而非单纯的功能丧失 (loss-of-function)。

## 7\. 文件系统与数据位置

**Yale HPC (Bouchet / Gibbs Cluster):**

/gpfs/gibbs/project/neurosplice\_gse117357/

├── raw\_data/                  (FASTQ 原始文件)

├── alignments/                (STAR BAM 文件)

└── rmats\_out/                 (rMATS v4.3.0 原始输出)

/palmer/scratch/\[PI\_NAME\]/neurosplice\_backup/  \<-- 外部硬盘数据的最终备份目标位置

**本地硬盘 (Mac Desktop):**

~/Desktop/WXY/NeuroSplice/results/

├── rmats\_gse117357\_analysis/  (AS 结果, 过滤后的 CSV)

├── deseq2/                    (DEG 输出 \- 已完成)

├── qc/multiqc\_align.html      (QC 质量报告)

└── figures/                   (现有的 Sashimi 图, 火山图)

**本地/共享硬盘 \[🚨 紧急警告 \- 这是外部硬盘\]:**

/Volumes/Untitled/NeuroSplice/  \<-- 包含 60GB 原始 BAM 文件；必须尽快备份到 HPC

## 8\. 待办事项与下一步计划

| 任务 | 描述 | 负责人 | 截止日期 | 状态 |
| :---- | :---- | :---- | :---- | :---- |
| **数据安全** | 立即将外部硬盘的 BAM 文件备份到 HPC 的 `palmer_scratch` 目录。 | 硬盘所有者 (协调: Jason) | 尽快 | 🔴 To Do |
| **元数据审计** | 通过 GEO/SRA 原始元数据验证文库链特异性 (`fr-firststrand`)。 | Jason | 5月11日 | ✅ Complete |
| **DEG 运行** | 针对 3 个脑区执行 `featureCounts` 和 DESeq2 分析 (`~genotype`)。 | Jason | 5月12日 | ✅ Complete |
| **Fisher 测试** | 使用 `baseMean ≥ 10` 作为背景，应用 Bonferroni 校正运行交集测试。 | Jason | 5月14日 | ✅ Complete |
| **虚拟 PCR** | 人类直系同源基因 + GWAS 已完成。PhyloP 得分需等待 HPC liftOver。 | Sarah | 5月15日 | 🟡 Partial |
| **图 1 草图** | 为 Pts 和 Neil2 生成 Sashimi 图。*(输入文件: `alignments/[sample].bam` \+ `rmats_out/[Event].MATS.JCEC.txt`)* | Sarah | 5月18日 | 🔴 To Do |
| **图 2 草图** | 为 Unc13b 和 Tox3 生成 Sashimi 图及相关验证指标。*(输入文件: 同图 1\)* | Sarah | 5月20日 | 🔴 To Do |

