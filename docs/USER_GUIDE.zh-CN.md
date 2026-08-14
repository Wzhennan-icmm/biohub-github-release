# BioHub v0.4 功能说明书

> 文档版本：0.4.0<br>
> 最后更新：2026-08-14<br>
> 适用对象：植物基因组组装、注释、比较基因组、群体分析和表达分析流程的使用者与维护者

BioHub 把历史上分散的 Python、R、AWK 和 shell 工具整合为单一命令行程序。核心数据处理由 Rust 实现；绘图和科研流程通过明确后端契约调用 Rscript、Snakemake 及领域软件。本说明书描述当前 `0.4.0` 代码的真实接口、输入输出、依赖、限制和验证状态。

重要：目录中的 `implemented` 表示脚本 ID 已有可调用实现，不等于算法已在所有物种、文件变体和历史数据上完成科研等价性验证。发表或生产使用前，应使用本项目夹具、自己的去标识化代表数据和经领域专家批准的 golden output 复核结果。

阅读导航：先看[安装与启动](#2-安装与启动)和[命令模型](#3-命令模型)；按场景上手可看[快速教程](#4-快速教程)；查参数与边界看[57 个目录命令](#6-57-个目录命令)和[直接命令组](#7-直接命令组)；逐项运行科研流程请看[13 个 Recipe 使用与发表复现手册](RECIPES.zh-CN.md)；发表前检查[验证状态](#8-验证状态)与[可重复性、数据安全与引用](#10-可重复性数据安全与引用)。GitHub 页面右侧大纲可直接定位具体命令。

## 1. 软件定位与功能边界

BioHub 当前提供：

- 57 个可通过目录发现的迁移命令 ID；
- 13 个带 schema、输入清单、日志、摘要和归档的 Snakemake recipe；
- FASTA、GFF/GTF、BLAST6、VCF、BAM、PAF、MUMmer coords、表达矩阵等常见文本或二进制生物信息格式处理；
- 最长转录本、注释格式转换、互惠最佳命中、四倍简并位点、深度图和点阵图等工作流；
- 统一机器可读 JSON 目录、按 recipe 外部依赖预检、容器构建和跨平台 CI；
- 旧命令入口兼容层，便于历史流程逐步迁移。

BioHub 当前不提供：

- 自动推断集群账户、分区、内存、运行时间或研究设计；
- FASTA/GFF/VCF 的完整标准校验器；
- `annotation-vcf` 中的密码子变化、氨基酸变化或变异致病性预测；
- 对所有历史脚本的逐字节输出顺序保证；
- 未经项目数据验证即可直接用于发表结论的科研保证。

## 2. 安装与启动

### 2.1 从源码构建

要求 Rust stable 工具链和 Cargo。仓库根目录执行：

```bash
git clone https://github.com/Wzhennan-icmm/biohub-github-release.git
cd biohub-github-release/biohub-rs
cargo build --release --locked
./target/release/biohub --version
./target/release/biohub --help
```

后文使用 `biohub` 作为命令名。可将 `biohub-rs/target/release` 加入 `PATH`，或把示例中的 `biohub` 替换为完整路径。

`--locked` 强制使用仓库提交的 `Cargo.lock`，避免依赖解析随时间漂移。发布资产可用时，应同时核对其 SHA256 文件。

### 2.2 Docker

仓库根目录执行：

```bash
docker build -t biohub:0.4.0 .
docker run --rm biohub:0.4.0 --version
docker run --rm biohub:0.4.0 catalog --format json
docker run --rm -v "$PWD:/work" -w /work biohub:0.4.0 doctor
```

核心容器包含 Rscript、samtools、MAFFT、`pal2nal.pl`、recipe 定义和示例，但不宣称包含全部领域软件。运行 recipe 应使用目录记录的领域镜像或经过锁定的本地环境。处理本地文件时必须挂载目录；容器退出后，仅挂载目录内的输出会保留。

### 2.3 外部依赖

| 依赖 | 使用功能 | 是否必需 |
| --- | --- | --- |
| Rust/Cargo | 源码构建 | 仅源码安装必需 |
| Rscript | `dotplot` | 按命令必需 |
| samtools | 3 个双端 BAM 合并命令 | 按命令必需 |
| mafft | `orthofinder-to-pal2nal` | 按命令必需 |
| pal2nal.pl | `orthofinder-to-pal2nal` | 按命令必需 |
| hamstr | 预留适配器 | 可选；当前目录命令不直接调用 |
| Snakemake | 13 个 recipe | recipe 必需 |
| Python 3 | recipe 校验、汇总和 provenance 辅助脚本 | recipe 必需 |
| cafe5、codeml、mcmctree、minimap2、syri、plink2、vcftools、Python 2 | 对应比较/群体 recipe | 按 recipe 必需 |
| DESeq2、vegan | RNA-seq、RDA recipe | 按 recipe 必需 |

运行预检：

```bash
biohub doctor
biohub doctor --json
biohub doctor --strict
biohub doctor --recipe selection-branch-site
```

普通模式报告缺失项但可返回成功；`--strict` 在任何检查项缺失时返回非零，适合 CI。`--json` 便于流程解析。

## 3. 命令模型

### 3.1 推荐入口

```text
biohub catalog [--format table|json] [--kind all|command|recipe]
biohub run <script-id> [options] [--force]
biohub doctor [--strict] [--json] [--recipe <recipe-id>]
biohub recipe list|describe|init|validate|run|report
biohub r list
biohub r run dotplot ...
```

发现所有迁移命令：

```bash
biohub catalog
biohub catalog --format json
biohub run annotation-vcf --help
```

JSON 目录每项包含 `id`、`kind`、领域、描述、状态、后端、依赖、版本和许可证；recipe 还含工作流、schema 和容器。用 `--kind command` 或 `--kind recipe` 过滤。`biohub run <id> --help` 当前显示元数据，不展开全部参数；参数以本说明书和实际调用错误信息为准。

### 3.2 兼容入口

```text
biohub scripts catalog
biohub scripts run <script-id> ...
```

兼容入口计划保留至 v1.x。新流程应使用 `catalog` 和 `run`，便于获得输出保护。历史直接命令组仍可用，见第 7 节。

### 3.3 输入、输出和覆盖规则

- 文本输入通常接受空白分隔字段；GFF、VCF 和 PAF 的关键路径要求制表符分隔。
- FASTA 支持多行序列，输出通常归并为单行序列。
- 许多单输出命令省略 `-o/--output` 时写标准输出；可使用 `-o -` 明确写标准输出。
- `biohub run` 对经过统一输出写入器的主输出启用“已存在即失败”；传 `--force` 允许替换。
- 兼容入口 `biohub scripts run` 保留历史覆盖行为。
- 多文件或历史固定文件名命令可能直接创建/覆盖文件，不全部受统一保护。运行 `compare-busco-results`、`convert-lastz2-jcvi`、NextGenomics 坐标转换、BAM 合并及目录输出命令时，优先使用新的空目录。
- 基于哈希表聚合的部分命令不保证行顺序。若下游需要确定性比较，按业务键排序后再比对。
- 路径支持 `~` 展开；输入文本应为 UTF-8 或 ASCII。

### 3.4 返回码

| 返回码 | 含义 |
| --- | --- |
| `0` | 成功，包括帮助或非严格依赖报告 |
| `1` | 输入、解析、I/O、外部程序、未知命令或科研前置验证失败 |
| `2` | 部分命令的用法或参数枚举错误 |

批处理必须检查返回码和标准错误，不能仅检查输出文件是否存在。`orthofinder-to-pal2nal` 可能产生部分结果后因某一 orthogroup 失败而整体返回非零。

## 4. 快速教程

### 4.1 提取每个基因最长转录本

仓库含合成夹具 `examples/longest-transcript/input.fa`：

```bash
biohub fasta longest-transcript \
  -f examples/longest-transcript/input.fa \
  -o longest.fa
```

该稳定入口优先读取 FASTA 标题中的 `gene=` 或 `gene_name=`；缺失时再按标题结构回退分组。每组选择序列长度最大记录，等长时保留先出现记录。示例结果包含 `txA.2` 和 `txB.1`。

核对：

```bash
rg '^>' longest.fa
```

目录中的 `get-longest-transcript` 使用不同的“点号后缀”分组规则，不能与此入口无条件互换。

### 4.2 从 PAF 或 MUMmer coords 绘制点阵图

先检查 Rscript：

```bash
biohub doctor
```

PAF：

```bash
biohub run dotplot \
  --input examples/dotplot/input.paf \
  --output dotplot-paf.pdf \
  --format paf
```

MUMmer `show-coords` 风格输入：

```bash
biohub run dotplot \
  --input examples/dotplot/input.coords \
  --output dotplot-coords.png \
  --format coords
```

省略 `--format` 时，`.paf` 后缀推断为 PAF，其他后缀推断为 coords。输出只允许 `.pdf` 或 `.png`。PAF 第 5 列 `-` 显示为反向颜色；当前 coords 解析仅取每行前四个数值坐标，统一按正向着色，因此方向敏感的 MUMmer 解释应先转为 PAF 或人工复核。

### 4.3 PASA/GeMoMa 注释整理与坐标转换

提取 PASA 结果中的蛋白记录和每个 gene 块内记录数最多的 mRNA 块：

```bash
biohub run extract-pasa-results \
  --input pasa.gff3 \
  --out-seq pasa.protein.fa \
  --out-gff pasa.longest.gff3
```

转换 GeMoMa 记录语义：

```bash
biohub run convert-gemoma-gff3 \
  --input gemoma.gff \
  --output gemoma.converted.gff3
```

将 PASA contig 坐标映射到染色体：

```bash
biohub run convert-gene-annotation-contigs2chr-PASA \
  -gff pasa.gff3 \
  --background placement.tsv \
  --output pasa.chromosome.gff3
```

这些转换依赖历史背景表列序和 GFF feature 顺序，不是通用格式转换器。首次处理新来源时，应抽取正链、负链、边界跨越、多转录本和未映射 contig 各一例，人工核对坐标、链方向和 Parent/ID 关系。

### 4.4 从双向 BLAST6 获取互惠最佳命中

推荐按 BLAST6 第 12 列 bit score 分别选择双向最佳命中，再保留互惠对：

```bash
biohub run get-best-hit-by-score \
  --query A_vs_B.blast6 \
  --refs B_vs_A.blast6 \
  --output A_B.rbh.tsv
```

两个输入每行至少 12 列；第 1、2、12 列分别视为 query、target、score。输出两列 `A_gene<TAB>B_gene`。分数并列时保留先读到的记录。若输入已经预先压缩为每个 query 一行，也可用：

```bash
biohub run compare-two-blast \
  --input A_to_B.best.tsv \
  --resBlastRes B_to_A.best.tsv \
  --output A_B.reciprocal.tsv
```

`compare-two-blast` 不按得分选择；同一 query 重复时最后一行覆盖先前值。

### 4.5 按 GFF 区间注释 VCF

```bash
biohub run annotation-vcf \
  --reference reference.fa \
  --gff genes.gff3 \
  --vcf calls.vcf \
  --format tsv \
  --output calls.annotated.tsv
```

输出列为：

```text
chrom pos ref alt type status gene_ids transcript_ids feature
```

`type` 为 SNV、INS、DEL 或 MNV；`status` 为 CDS、Exon、Intron、Gene-body 或 Intergenic。多等位 ALT 当前只使用第一个等位。参考 FASTA 用于限定合法染色体名，不核对 REF 碱基，也不计算同义/非同义或氨基酸变化。

JSON 输出：

```bash
biohub run annotation-vcf \
  -r reference.fa -g genes.gff3 -v calls.vcf \
  -f json -o calls.annotation.json
```

JSON 模式会创建 JSON 文件；部分扩展名组合会同时创建 TSV/文本主文件和同前缀 `.json`。建议明确使用 `.json`，并在下游记录实际产物清单。

### 4.6 MAFFT + PAL2NAL 密码子比对

输入：一个目录内每个文件代表一个蛋白 orthogroup；一个包含对应 CDS 的 FASTA。

```bash
biohub doctor --strict
biohub run orthofinder-to-pal2nal \
  --pathOfprot orthogroup_proteins/ \
  --outPutPath codon_alignments/ \
  --nuclOfcds all.cds.fa
```

BioHub 在调用外部程序前检查：

- CDS FASTA ID 唯一；
- 每个蛋白 ID 能找到对应 CDS；
- CDS 非空且长度可被 3 整除；
- 每组蛋白数和选中 CDS 数一致。

每组成功时输出蛋白比对、`*.cds.fasta` 和 `*.codon.paml`。目录还包含 `validation_summary.tsv` 与 `README.txt`。MAFFT 参数固定为 `--maxiterate 1000 --localpair`；PAL2NAL 参数固定为 `-output paml -nogap`。任一组跳过或失败时命令返回非零；保留诊断表和已完成组。推荐每次使用空输出目录。

### 4.7 Recipe 科研工作流

Recipe 把多步分析固定为可检查的 Snakemake DAG。模板中的 `REQUIRED` 和
`null` 是待填写占位符，不是默认科研阈值。先复制模板，再填写项目路径、
样本/物种清单、阈值、重复次数、随机种子和资源参数：

本节提供统一入口和流程概览。13 个 recipe 的 manifest 表头、全部配置字段、产物表、
恢复方式和论文 Methods 核对项见[专用复现手册](RECIPES.zh-CN.md)。

```bash
biohub recipe list
biohub recipe describe selection-branch-site
biohub recipe init selection-branch-site --workdir analysis-config
# 编辑 analysis-config/config.yaml
biohub doctor --recipe selection-branch-site --strict
biohub recipe validate selection-branch-site \
  --config analysis-config/config.yaml
biohub recipe run selection-branch-site \
  --config analysis-config/config.yaml \
  --workdir runs/selection-20260814 --cores 8
biohub recipe report --workdir runs/selection-20260814
```

默认 `local` profile。Slurm 使用：

```bash
biohub recipe run <recipe-id> --config config.yaml --workdir run-dir \
  --profile slurm --cores 8
```

内置 Slurm profile 依赖 Snakemake 8.6+ 与
`snakemake-executor-plugin-slurm`，为安全起见默认只提交一个作业。生产运行前
复制 profile，填写集群批准的并发数、账户、分区、内存和运行时间。也可把
自定义 profile 目录传给 `--profile`。

每个 run 目录至少包含：

- `config.resolved.yaml`、`config.sha256`、`recipe.id` 和 `command.sh`；
- `run.json`、`provenance.json`、`versions.tsv`、`recipe.sources.sha256` 和 `inputs.manifest.tsv`；
- `logs/`、`results/`、`report/`、压缩归档及成功运行后的 `checksums.sha256`。

现有 run 目录默认拒绝覆盖。`--resume` 仅允许 recipe ID 和配置 SHA256 与原运行
一致；修改配置必须新建 run 目录。失败状态写入 `run.json`，修复原因后可用原配置
恢复。`validate` 执行 Snakemake dry-run，只检查 DAG/输入和 schema，不代替领域
结果验证。`recipe report` 已存在 `report/report.html` 时同样拒绝覆盖；确认替换后传
`--force`。

`versions.tsv` 对目录声明的每项依赖记录版本、不可用状态，或“已安装但无安全版本
探针”；不会通过启动 codeml/MCMCTree 等交互程序猜测版本。`provenance.json` 记录
workflow SHA256 与推荐容器；容器运行时应设置 `BIOHUB_CONTAINER_DIGEST`。成功后的
`recipe.sources.sha256` 覆盖 recipe 内工作流、schema、脚本及共享 provenance helper；
`checksums.sha256` 覆盖稳定 run 文件，故意排除后续状态可能变化的 `run.json`。

| Recipe ID | 主要输入与检查 | 核心产物 | 领域依赖 |
| --- | --- | --- | --- |
| `comparative-orthology-codon` | orthogroup 蛋白目录、CDS FASTA、预期 taxa；检查 ID、数量和 CDS 三联体 | 蛋白/CDS/无 gap PAML 密码子比对、验证汇总 | MAFFT、PAL2NAL |
| `gene-family-cafe` | gene count、Newick 树；检查 taxa 一致、二叉树、分支长度、超度量容差和 family 过滤 | 重复 CAFE 运行、似然/参数稳定性选择、扩张收缩表 | CAFE5 |
| `selection-branch-site` | test manifest、PAML 比对、仅含一个 `#1` 的前景树 | 配对 Model A/null、LRT、50:50 混合卡方 p、全局 BH、BEB 位点 | codeml |
| `dating-mcmctree` | replicate manifest、stage1/stage2 ctl；检查 `ndata` 和 `usedata=3→2` | 隔离重复、节点后验摘要、分位数和 ESS | mcmctree |
| `assembly-t2t-evaluate` | assembly manifest、参考 FASTA；显式预期染色体和端粒 motif 参数 | N50/N/长度审计、端粒末端证据、PAF 覆盖/identity | minimap2 |
| `synteny-sv` | 成对 assembly manifest；可要求序列 ID 集一致 | 原始/过滤 PAF、SyRI 结果、SV 类型汇总 | minimap2、SyRI |
| `population-gwas` | VCF/PGEN/BED、trait manifest、可选 covariate；检查样本覆盖 | 每 trait PLINK2 结果、状态与机械最小 p 汇总 | PLINK2 |
| `population-selection` | VCF、两群体比较 manifest；检查样本存在且集合不相交 | windowed Weir FST、各群体 nucleotide diversity 汇总、候选窗口 | VCFtools |
| `kmer-gwas` | trait manifest、预制 k-mer table prefix、显式 legacy 脚本路径 | 每 trait 关联结果、显著行和候选 k-mer FASTA | Python 2、外部 k-mer GWAS 脚本 |
| `family-denovo-rate` | family/pair manifest、candidate TSV、callable BED；检查坐标和 callable 覆盖 | candidate 审计、pair/combined rate、Garwood exact Poisson CI | Python 3 标准库 |
| `rnaseq-deseq2` | 原始整数 count、sample design、contrast manifest；检查样本集合和设计满秩 | normalized counts、全量差异表、contrast 汇总、PCA | R、DESeq2 |
| `functional-enrichment` | foreground set、background、显式 term association/source | 单侧超几何检验、按配置范围 BH、覆盖审计和图 | base R |
| `microbiome-rda` | feature table、metadata、constraint/Condition 列；检查样本和缺失值 | 过滤审计、RDA scores、overall/term/axis permutation tests | R、vegan |

这些 recipe 当前状态为 `experimental`。统计参数可复现，不表示数据设计、校准点、
前景分支、群体定义、候选阈值或因果解释已获领域批准。发表前必须用去标识化代表
数据、独立工具输出和专家批准的 golden result 核对。

## 5. 全量命令参考阅读方法

本节每项采用相同字段：

- **用途**：实现做什么；
- **调用**：推荐的最短完整命令；
- **输入/输出**：必需列、默认输出和格式；
- **依赖/注意**：外部程序、排序、兼容行为和科研限制。

所有条目均可由 `biohub catalog --format json` 发现。未显式写外部依赖的条目使用 Rust 标准实现。

## 6. 57 个目录命令

### 6.1 可视化

### `dotplot`

- **用途**：把 PAF 或 MUMmer coords 比对坐标绘制为静态点阵图。
- **调用**：`biohub run dotplot --input <align.paf|coords> --output <plot.pdf|png> [--format paf|coords] [--force]`
- **输入/输出**：PAF 至少 9 列；coords 每个有效行至少能解析 4 个数值。输出扩展名决定 PDF 或 PNG。
- **依赖/注意**：依赖 Rscript，仅使用 base R。空输入、无可解析行、非法格式或非 PDF/PNG 输出失败。PAF 支持正反向着色；coords 当前不解析方向。

### `psmc-plot`

- **用途**：把 `biohub psmc merge` 生成的多个样本轨迹绘制为静态 PSMC 人口史图。
- **调用**：`biohub run psmc-plot --input <merged.tsv> --output <plot.pdf|png> [--x-scale linear|log10] [--y-scale linear|log10] [--stages stages.tsv] [--force]`
- **输入/输出**：主表必需列为 `Sample Time Ne`。可选 stages 表必需列为 `label start end color`；区间只做背景注释。输出扩展名决定 PDF/PNG。
- **依赖/注意**：依赖 Rscript，仅使用 base R；两个坐标轴默认 `log10`。log10 模式要求正值。旧脚本中的地质年代、颜色和坐标范围已移除；不会自动解释人口史事件。

### `plot-depth-pandepth`

- **用途**：解析 Pandepth 窗口表，按染色体长度筛选并绘制 GC—深度散点图。
- **调用**：`biohub run plot-depth-pandepth -i <pandepth.tsv> -o <out-dir> [-l <min-Mb>]`
- **输入/输出**：默认列为 `Chr Start End MeanDepth GC(%)`；也可用 `#` 开头的同名表头定位列。输出 `chromosome_stats.tsv`、`filtered_depth.tsv`、`depth_gc_scatter.svg`。
- **依赖/注意**：`-l` 默认 `10.0` Mb。没有有效行或没有达到长度阈值的染色体时失败；绘图最多抽样约 25,000 点。

### `plot-depth-pandepth2`

- **用途**：使用与 `plot-depth-pandepth` 相同数据，输出更密集采样和不同点样式。
- **调用**：`biohub run plot-depth-pandepth2 -i <pandepth.tsv> -o <out-dir> [-l <min-Mb>]`
- **输入/输出**：统计表与过滤表同上，图文件为 `depth_gc_styled.svg`。
- **依赖/注意**：默认最短染色体 10 Mb；绘图最多抽样约 40,000 点。它不是新统计模型，仅是样式变体。

### `plot-mosdepth-point`

- **用途**：将 mosdepth 风格 `chrom start end coverage` 窗口转换为累计染色体坐标散点图。
- **调用**：`biohub run plot-mosdepth-point -i <regions.tsv> -o <out-dir> [-l <min-bp>]`
- **输入/输出**：输入至少 4 列。输出 `mosdepth_points.tsv` 和 `mosdepth_scatter.svg`。
- **依赖/注意**：`-l` 默认 `0` bp；按染色体名称排序，染色体间加入显示间隔；图最多抽样约 60,000 点。

### 6.2 序列与标识符处理

### `change-scaffolds-name`

- **用途**：替换制表文本第一列的 scaffold/contig ID。
- **调用**：`biohub run change-scaffolds-name -i <table.tsv> -l <map.tsv> [-o <output>]`
- **输入/输出**：映射表为 `旧ID 新ID`；输入第一列匹配时替换，未匹配行原样保留。省略输出时写 stdout。
- **依赖/注意**：只改第一列；空行跳过。该目录命令的映射方向与直接命令 `biohub rename scaffolds` 相反，迁移旧流程时必须核对。

### `change-scaffolds-name-fasta`

- **用途**：按映射表筛选并重命名 FASTA 记录。
- **调用**：`biohub run change-scaffolds-name-fasta -i <input.fa> -l <map.tsv> [-o <output.fa>]`
- **输入/输出**：映射表为 `新ID 旧ID`。仅输出映射表中存在且 FASTA 可找到的旧 ID，顺序跟随映射表；序列写为单行。
- **依赖/注意**：会丢弃未映射 FASTA 记录和原标题描述。若要求保留未映射记录，使用直接入口 `biohub rename fasta-scaffolds` 并采用其 `旧ID 新ID` 映射方向。

### `change-seqname-for-fasta`

- **用途**：批量重命名目录内所有 `*.aln` FASTA 标题。
- **调用**：`biohub run change-seqname-for-fasta -i <input-dir> -l <map.tsv> -o <output-dir>`
- **输入/输出**：映射表为 `新ID 旧ID`；每个 `<name>.aln` 输出为 `<name>.newName.fa`。未映射标题保留原 ID。
- **依赖/注意**：只扫描输入目录第一层，区分文件名大小写；现有输出可能被覆盖，建议空目录运行。

### `filter-seq-by-length`

- **用途**：保留长度大于或等于阈值的 FASTA 记录。
- **调用**：`biohub run filter-seq-by-length -i <input.fa> -l <minimum-length> [-o <output.fa>]`
- **输入/输出**：长度按拼接后的序列字符数计算；输出保留标题并把序列写为单行。
- **依赖/注意**：阈值必须是非负整数；不会识别 gap 或模糊碱基的特殊含义，它们同样计入长度。

### `genome-gc`

- **用途**：计算整个 FASTA 的 GC 比例。
- **调用**：`biohub run genome-gc -f <genome.fa> [-o <ratio.txt>]`
- **输入/输出**：输出单个 `GC/(A+T+C+G)` 小数；忽略标题、N 和其他字符。
- **依赖/注意**：没有 ATCG 时输出 `0`。这是全基因组汇总，不提供逐序列或窗口统计。

### `get-the-longest-seq`

- **用途**：按历史点号标题规则选择每个基因最长记录，并输出名称对应关系。
- **调用**：`biohub run get-the-longest-seq -i <protein.fa> [-o <mapping.tsv>]`
- **输入/输出**：标题预期形如 `transcript.gene...`；使用第 2 个点号字段作为 gene、第 1 个字段作为 transcript，输出 `transcript<TAB>gene`，不输出序列。
- **依赖/注意**：点号不足的标题被忽略；聚合输出顺序不保证。名称规则不符合时应使用稳定入口或先标准化标题。

### `get-longest-transcript`

- **用途**：删除标题最后一个点号后缀作为 gene 分组，输出每组最长 FASTA 记录。
- **调用**：`biohub run get-longest-transcript -i <input.fa> [-o <longest.fa>]`
- **输入/输出**：例如 `GeneA.iso1` 和 `GeneA.iso2` 分为同组；输出完整原标题和单行序列。
- **依赖/注意**：无点号标题被忽略；等长保留先出现记录；组和输出顺序不保证。不要与能识别 `gene=` 的 `biohub fasta longest-transcript` 混淆。

### `extract-longest-pep`

- **用途**：按历史 Ensembl 下载标题规则选取最长肽序列。
- **调用**：`biohub run extract-longest-pep -f <proteins.fa> [-o <longest.fa>]`
- **输入/输出**：把标题第一个 `-` 前内容视为分组名，输出该组最长序列，标题改为分组名。
- **依赖/注意**：该规则适合原脚本数据，不是通用 Ensembl header 解析。输出顺序不保证，重复标题处理保留历史兼容逻辑。

### `extract-gene-family-info-alt`

- **用途**：`extract-longest-pep` 的兼容别名。
- **调用**：`biohub run extract-gene-family-info-alt -f <proteins.fa> [-o <longest.fa>]`
- **输入/输出**：完全复用 `extract-longest-pep` 的输入、输出和 `-` 分组规则。
- **依赖/注意**：新流程应使用语义更明确的 `extract-longest-pep`；此 ID 仅用于旧脚本兼容。

### `trim-ttaggg-fastq`

- **用途**：丢弃序列头端或尾端匹配端粒 motif 及其反向互补循环移位的 FASTQ reads。
- **调用**：`biohub run trim-ttaggg-fastq -i <reads.fastq> -o <filtered.fastq> [-s <motif>]`
- **输入/输出**：默认 motif 为 `TTAGGG`；按四行 FASTQ 记录读取，保留未命中的完整记录。
- **依赖/注意**：当前只比较 read 首尾各 6 bp，因此自定义非 6 bp motif 不会得到完整泛化行为。`.gz` 输入或输出明确不支持；截断的最后一条记录会被忽略。

### 6.3 表格、集合、BUSCO 和 GO

### `convert-3line2one`

- **用途**：把连续三行结构变异深度记录合并为一行。
- **调用**：`biohub run convert-3line2one -i <depth.tsv> [-o <windows.tsv>]`
- **输入/输出**：每个有效输入行至少 4 列；每组三行输出首行前 3 列，再依次输出三行第 4 列，共 6 列。
- **依赖/注意**：不足三行的尾组丢弃；格式不合格行跳过并会改变分组计数语义，运行前应先校验每行列数。

### `merge-two-txt`

- **用途**：按第一列去重，保留第三列数值最大的一整行。
- **调用**：`biohub run merge-two-txt -i <combined.tsv> [-o <best.tsv>]`
- **输入/输出**：输入至少 3 列；第三列无法解析时按 `0` 处理。输出字段以制表符连接。
- **依赖/注意**：并列保留先出现行；输出键顺序不保证。名称源自历史脚本，实际只接收一个已合并输入文件。

### `compare-as-and-noAS`

- **用途**：按 noAS 文件第一列筛选 AS 表中具有对应条目的行。
- **调用**：`biohub run compare-as-and-noAS -nA <noAS.tsv> -AS <AS.tsv> [-o <matched.tsv>]`
- **输入/输出**：两个输入至少 3 列；输出满足第一列存在于 noAS 且双方第三列非空的 AS 原始行。
- **依赖/注意**：空白切分后第三列天然非空，因此核心条件近似为第一列交集。需用历史样例确认该兼容语义符合预期。

### `compare-busco-results`

- **用途**：比较两组 BUSCO 条目集合。
- **调用**：`biohub run compare-busco-results -a <ancestor.txt> -o <offspring.txt>`
- **输入/输出**：每个非空整行视为一个集合元素；固定在当前目录写 `justHave.txt`（后代独有）和 `justLost.txt`（祖先独有）。此处 `-o` 是第二个输入，不是输出路径。
- **依赖/注意**：固定文件名可能覆盖已有结果，且 `--force` 不保护它们。必须在独立工作目录运行。

### `save-go`

- **用途**：从注释表每行提取包含字符串 `GO` 的字段。
- **调用**：`biohub run save-go -i <annotations.tsv> [-o <go.tsv>]`
- **输入/输出**：第一列为 gene，其余字段中含 `GO` 的项原样输出；没有 GO 时仍输出 gene 和空第二部分。
- **依赖/注意**：不校验 `GO:NNNNNNN` 格式，不拆分字段内部的多个 term，也不去重。

### `merge-gos`

- **用途**：合并 Swiss-Prot、NR 和 TrEMBL 来源的 gene—GO 二列表。
- **调用**：`biohub run merge-gos -s <swiss.tsv> -n <nr.tsv> -T <trembl.tsv> [-o <merged.tsv>]`
- **输入/输出**：每行使用前两列；以 `EC` 开头的第二列被排除；输出每个 gene—term 一行。
- **依赖/注意**：同一来源内重复项去重，但跨来源重复项可能保留；输出顺序不保证。输入打开失败的来源会被静默跳过，发表流程应先检查文件存在性。

### `orthogenes`

- **用途**：清理 `group:member1,member2,...` 列表中的重复成员。
- **调用**：`biohub run orthogenes -i <groups.txt> -o <deduplicated.txt>`
- **输入/输出**：保持成员首次出现顺序；仅输出去重后至少 2 个成员且 group 名非空的行。
- **依赖/注意**：单成员组被主动删除；输出参数必需。

### `zhouxiaoxuan-mergexls`

- **用途**：按历史基因名规范把两个空白分隔的 XLS-like 表连接。
- **调用**：`biohub run zhouxiaoxuan-mergexls -a <metadata.tsv> -b <names.tsv> [-o <merged.tsv>]`
- **输入/输出**：第一表首列建立键与附加列；第二表每个整行按第一个点号前缀查找。匹配时输出 `第二表原行、第一表完整名、附加列`；未匹配行原样输出。
- **依赖/注意**：键解析依赖首表名称中 `:` 与 `.` 的历史命名结构，非通用 join。

### `check-duplication-gene-pairs`

- **用途**：对称去重基因对。
- **调用**：`biohub run check-duplication-gene-pairs -i <pairs.tsv> [-o <unique.tsv>]`
- **输入/输出**：只使用前两列；首次出现 `A B` 后，后续 `A B` 和 `B A` 都丢弃。输出两列。
- **依赖/注意**：丢弃输入第 3 列及以后信息；结果保留首次方向。

### 6.4 BLAST、同源关系与 JCVI

### `compare-two-blast`

- **用途**：从两个方向的二列最佳命中表中保留互惠关系。
- **调用**：`biohub run compare-two-blast -i <A-to-B.tsv> -r <B-to-A.tsv> [-o <rbh.tsv>]`
- **输入/输出**：每行至少两列 query、target；输出互惠的两列关系。
- **依赖/注意**：同一 query 重复时最后一条覆盖之前值；命令本身不比较 score。先为每个 query 选择最佳命中。

### `get-best-idy`

- **用途**：对 BLAST-like 表按 query 汇总最大 identity。
- **调用**：`biohub run get-best-idy -i <blast.tsv> [-o <best-identity.tsv>]`
- **输入/输出**：使用第 1 列 query 和第 3 列 identity；输出两列 query、最大 identity。
- **依赖/注意**：不输出 target，因此不能单独构造同源对；非法 identity 按 `0`。输出顺序不保证。

### `get-best-hit-based-on-idy`

- **用途**：按 BLAST6 第 12 列 bit score 为每个 query 保留整条最佳记录。
- **调用**：`biohub run get-best-hit-based-on-idy -i <blast6.tsv> [-o <best.tsv>]`
- **输入/输出**：输入至少 12 列；输出原始最佳行。
- **依赖/注意**：名称保留历史叫法，但实际选择依据是第 12 列，不是第 3 列 identity；并列保留先出现行，输出顺序不保证。

### `get-best-hit-genes`

- **用途**：基于 BLAST6 第 11 列 e-value 统计双向最佳和次佳命中，并输出互惠最佳对的诊断字段。
- **调用**：`biohub run get-best-hit-genes -i <combined-blast6.tsv> [-o <hits.tsv>]`
- **输入/输出**：同一输入同时构建两个方向；输出 8 列，包括 gene、最佳 target、双方最佳/次佳 e-value 和次佳 target。
- **依赖/注意**：要求制表符且至少 11 列。一个输入被视为无向候选集合；输出顺序不保证，使用前应确认与 BLAST 生成方式匹配。

### `merge-blastp-best-jcvi`

- **用途**：把目录内多个二列表合并成 JCVI 风格 best-hit 矩阵。
- **调用**：`biohub run merge-blastp-best-jcvi -p <tables-dir> [-o <matrix.tsv>]`
- **输入/输出**：每文件前两列为 orthogene、gene；文件按名称排序，各文件形成一列，缺失填 `.`。
- **依赖/注意**：表头不自动生成；目录内所有普通文件都会读取。gene 首次出现顺序决定行顺序。

### `convert-lastz2-jcvi`

- **用途**：把 LASTZ 块切分为 JCVI 可消费的伪基因 BED-like 表和 simple 配对表。
- **调用**：`biohub run convert-lastz2-jcvi -i <blocks.bed> -r <ref.len> -q <query.len> [--refName Ref] [--queryName Query]`
- **输入/输出**：块表至少 7 列：query chr/start/end、reference chr/start/end、strand。固定在当前目录写 `query.bed`、`ref.bed`、`ref_query.simple`。
- **依赖/注意**：当前实现为兼容历史逻辑，`-r/-q` 文件参数被要求但内容未使用；以 100 bp 片段和 500 bp gap 构造伪基因。固定输出可能覆盖，必须空目录运行。

### `get-best-hit-by-score`

- **用途**：按两个方向 BLAST6 的第 12 列 score 选择互惠最佳命中。
- **调用**：`biohub run get-best-hit-by-score -i <A-vs-B.blast6> -r <B-vs-A.blast6> [-o <rbh.tsv>]`
- **输入/输出**：输入各至少 12 列；输出 `A_query<TAB>B_target`。
- **依赖/注意**：非法 score 视为负无穷；并列保留先出现记录；输出顺序不保证。

### `get-best-hit-by-score-one-file`

- **用途**：从单个 BLAST6 表同时构建 query→target 和 target→query 最佳 score，输出互惠对及 identity。
- **调用**：`biohub run get-best-hit-by-score-one-file -i <blast6.tsv> -p <prefix> [-o <output>]`
- **输入/输出**：默认输出 `<prefix>.idy.txt`，三列 query、target、第 3 列 identity。
- **依赖/注意**：默认输出由函数直接创建，覆盖保护不完整；建议总是传 `-o` 且使用新路径。

### `get-best-hit-from-blast`

- **用途**：从多物种二列关系文件目录筛选每个基因恰有 `物种数-1` 个命中且关系互相可追溯的集合。
- **调用**：`biohub run get-best-hit-from-blast -i <relations-dir> -n <species-count> [-o <groups.tsv>]`
- **输入/输出**：读取目录内所有普通文件的前两列；输出中心 gene 后接各 target。
- **依赖/注意**：不识别物种列，只按命中数量判断；重复关系会计数；代表组去重和输出顺序受哈希迭代影响。需用已知 orthogroup 验证。

### `get-best-hit-for-scan`

- **用途**：`get-best-hit-from-blast` 的兼容别名。
- **调用**：`biohub run get-best-hit-for-scan -i <relations-dir> -n <species-count> [-o <groups.tsv>]`
- **输入/输出**：与主命令完全一致。
- **依赖/注意**：新流程应使用 `get-best-hit-from-blast`，避免把别名误解为独立算法。

### 6.5 基因注释与坐标转换

### `filter-gff-by-id`

- **用途**：按 ID 列表保留 GFF 行。
- **调用**：`biohub run filter-gff-by-id -gff <input.gff3> -id <ids.txt> [-o <filtered.gff3>]`
- **输入/输出**：ID 文件每个非空整行作为匹配值；GFF 注释列按分号切分后只取第一个属性 token 做完全匹配；注释行原样保留。
- **依赖/注意**：若属性为 `ID=gene1;...`，列表必须含 `ID=gene1`，而不是 `gene1`。不递归保留 Parent/child；需提前准备完整目标属性集合。

### `filter-gtf-ctg`

- **用途**：从 GTF/GFF-like 表中排除指定 contig。
- **调用**：`biohub run filter-gtf-ctg -i <input.gtf> -id <exclude.txt> -o <filtered.gtf>`
- **输入/输出**：排除列表每行一个完整 contig 名；匹配输入第一列即删除，其余行原样输出。
- **依赖/注意**：输出必需；注释行若第一列不在排除列表会保留。比较区分大小写。

### `filter-gff-by-fasta`

- **用途**：使用 FASTA `.fai` 或二列长度表删除越过序列边界的 GFF3 feature，并传播失效模型。
- **调用**：`biohub run filter-gff-by-fasta --gff <input.gff3> --fai <reference.fa.fai> [-o <filtered.gff3>]`
- **输入/输出**：长度文件前两列为 sequence ID、正整数长度。合法 GFF3 必须 9 列；注释、空行和非 9 列原样保留。输出统计写 stderr。
- **依赖/注意**：未知 sequence、`start=0`、`end<start` 或 `end>length` 视为失效。越界 gene/transcript 的子记录删除；越界 exon/CDS 会使其 Parent transcript 及后代失效。多 Parent 任一失效时保守删除整条记录。数值格式错误直接失败。

### `extract-pasa-results`

- **用途**：从 PASA 风格结果中提取 `#PROT` 蛋白和每个 gene 内记录数最多的 mRNA 块。
- **调用**：`biohub run extract-pasa-results -i <pasa.gff3> -s <proteins.fa> -g <longest.gff3>`
- **输入/输出**：GFF 输出保留 gene 及选中 mRNA 块；`#PROT name isoform sequence` 行生成 FASTA `>name-isoform`。
- **依赖/注意**：最长依据是 mRNA 块的记录行数，不是 CDS 碱基长度；要求记录按 gene/mRNA 分块排序。常规注释行前的非 `#PROT` 注释不写入 GFF。

### `convert-gemoma-gff3`

- **用途**：把历史 GeMoMa GFF 重新编号并展开 exon/CDS 语义。
- **调用**：`biohub run convert-gemoma-gff3 -i <gemoma.gff> [-o <converted.gff3>]`
- **输入/输出**：gene/mRNA 重新生成 `Plants<SEQ>geneNNNNN` 风格 ID；其他 feature 每行输出一个 exon 和一个保留原 feature 类型的 CDS-like 行。
- **依赖/注意**：依赖 gene→mRNA→子记录顺序；会替换原属性与 ID。当前 Parent 的历史命名组合需在下游加载前用 GFF validator 和人工样例核对。

### `convert-gene-annotation-contigs2chr-PASA`

- **用途**：按 PASA 历史 placement 背景表将 contig GFF 坐标映射到染色体。
- **调用**：`biohub run convert-gene-annotation-contigs2chr-PASA -gff <input.gff3> -b <placement.tsv> [-o <mapped.gff3>]`
- **输入/输出**：背景表至少 5 列，原第 4 列为 contig key；方向 `0` 使用正向偏移，方向 `1` 反向换算并翻转 `+/-`。未映射记录原样保留。
- **依赖/注意**：背景列位置为历史契约；数值解析失败会回退 `0`。必须对正反链边界做独立核对。

### `convert-gene-annotation-contigs2chr`

- **用途**：`convert-gene-annotation-contigs2chr-PASA` 的兼容别名。
- **调用**：`biohub run convert-gene-annotation-contigs2chr -gff <input.gff3> -b <placement.tsv> [-o <mapped.gff3>]`
- **输入/输出**：与 PASA 主命令完全一致。
- **依赖/注意**：新流程建议使用带 `-PASA` 的明确 ID，以记录背景表契约。

### `convert-gene-annotation-scaffold2chr-nextgenomics`

- **用途**：按 NextGenomics 历史布局表转换 scaffold 注释，并记录分割 scaffold、越界和未映射项目。
- **调用**：`biohub run convert-gene-annotation-scaffold2chr-nextgenomics -gff <input.gff3> -b <layout.tsv> [-o <mapped.gff3>]`
- **输入/输出**：默认主输出 `Change-annot.gff3`；当前目录另写 `change-gene-on-splitSca.txt`、`change-gene-on-splitSca1.txt`、`change-annot.log`、`change-errors-scaffolds.txt`。
- **依赖/注意**：布局表至少 7 个制表字段；固定辅助文件可覆盖。建议为每次运行新建目录，并对 split scaffold 诊断逐条处理。

### `convert-gene-annotation-legacy-alias`

- **用途**：NextGenomics 坐标转换的兼容别名。
- **调用**：`biohub run convert-gene-annotation-legacy-alias -gff <input.gff3> -b <layout.tsv> [-o <mapped.gff3>]`
- **输入/输出**：与 `convert-gene-annotation-scaffold2chr-nextgenomics` 相同，包括固定辅助文件。
- **依赖/注意**：仅保留旧流程可调用性；新分析应使用主 ID。

### `filter-gemoma-as`

- **用途**：按每个 gene 内累计 CDS 长度选择最长 GeMoMa mRNA 块。
- **调用**：`biohub run filter-gemoma-as -i <gemoma.gff3> [-o <longest.gff3>]`
- **输入/输出**：保留注释行和 gene 行，每个 gene 输出选中 mRNA 及其随后的记录。
- **依赖/注意**：输入必须按 gene/mRNA/CDS 连续排序。CDS 长度按 `end-start+1`；并列保留先出现转录本。

### `filter-gemoma-as2`

- **用途**：另一版 GeMoMa 最长 CDS isoform 选择器。
- **调用**：`biohub run filter-gemoma-as2 -i <gemoma.gff3> [-o <longest.gff3>]`
- **输入/输出**：同样按 gene 块输出一个 mRNA；保留注释行。
- **依赖/注意**：本版本 CDS 累计使用 `abs(end-start)`，少计每段 1 bp，与 `filter-gemoma-as` 不完全等价。用于发表前应选定一个实现并通过 golden output 固定。

### `annotation-vcf`

- **用途**：按 GFF 区间给 VCF 变异标记基因、转录本和 CDS/exon/intron/intergenic 区域。
- **调用**：`biohub run annotation-vcf -r <reference.fa> -g <genes.gff3> -v <calls.vcf> [-f tsv|json|pickle] [-o <output>]`
- **输入/输出**：默认 `annotation_vcf.txt`；TSV 含 9 列。JSON/pickle 请求写 JSON 数组，某些后缀组合还会保留主文件。
- **依赖/注意**：仅取多等位 ALT 第一项；不校验 REF，不进行密码子/蛋白效应预测；按 1-based 闭区间判断。GFF 支持 `key=value` 和 `key "value"` 属性的有限解析。

### 6.6 表达量与基因家族

### `merge-fpkm-file`

- **用途**：把目录内多个二列表达文件合并为矩阵，并计算每个 gene 的跨文件总和。
- **调用**：`biohub run merge-fpkm-file -i <fpkm-dir> -oF <matrix.tsv> -oP <profile.tsv>`
- **输入/输出**：每文件前两列为 gene、value；文件按名称排序。矩阵表头记录完整文件路径；profile 为 gene、求和值。
- **依赖/注意**：缺失 gene 不自动补零，可能造成不齐矩阵；目录内所有普通文件均参与。输入前应统一 gene 集并校验列数。

### `extract-gene-family-info`

- **用途**：按成员长度相对 parent 的覆盖比例，把家族表达记录拆为 full-length 和 partial 文件。
- **调用**：`biohub run extract-gene-family-info -l <gene-length.tsv> -e <expression.tsv> -p <parent.tsv> -c <ratio> [-o <out-dir>]`
- **输入/输出**：默认输出当前目录；每家族可生成 `<family>.partial.geneExpression.txt` 与 `<family>.fullenth.geneExpression.txt`。表达表首行跳过。
- **依赖/注意**：长度表使用第 1、3 列；parent 表使用第 1、3 列；gene 的第一个点号前缀作为 family。会排除以 `evm`、`gene`、`DMRT` 开头的表达 ID。`fullenth` 拼写为兼容固定名称。

### `extract-gene-family-matrix`

- **用途**：为指定家族计算 partial 成员和 full-length 成员的平均表达摘要。
- **调用**：`biohub run extract-gene-family-matrix -l <gene-length.tsv> -e <expression.tsv> -f <family-list.txt> -g <gene-family.tsv> -c <ratio> [-o <matrix.tsv>]`
- **输入/输出**：默认 `final.matrix.txt`；每行 `family partial_mean full_mean`。表达每个 gene 的多个数值先取均值；每家族最长成员作为 parent。
- **依赖/注意**：长度表这里使用第 1、2 列，与上一命令不同；gene-family 表为 `family gene`。缺失表达按 0，输出顺序不保证。

### 6.7 BAM 与 FASTQ

### `merge-two-end-bam`

- **用途**：按 read name 取两个 BAM 的交集，分别输出成对 R1/R2 BAM。
- **调用**：`biohub run merge-two-end-bam -i <R1.bam> -j <R2.bam> [-o1 <R1.out.bam>] [-o2 <R2.out.bam>]`
- **输入/输出**：默认文件名为 `test13.h1.R1.outReads.bam` 和 `test13.h1.R2.outReads.bam`；输出使用 R1 header。
- **依赖/注意**：依赖 samtools。R2 全部记录载入内存，同名记录只保留最后一条；不排序、不建索引。大 BAM 可能占用大量内存。

### `merge-two-end-bam1`

- **用途**：与 `merge-two-end-bam` 相同，使用更通用的历史默认文件名。
- **调用**：`biohub run merge-two-end-bam1 -i <R1.bam> -j <R2.bam> [-o1 <R1.out.bam>] [-o2 <R2.out.bam>]`
- **输入/输出**：默认 `R1.outReads.bam`、`R2.outReads.bam`。
- **依赖/注意**：依赖 samtools；匹配规则、内存和排序限制同上。

### `merge-two-end-bam-forMGI`

- **用途**：按 MGI 历史 read name 规则匹配双端 BAM。
- **调用**：`biohub run merge-two-end-bam-forMGI -i <R1.bam> -j <R2.bam> [-o1 <R1.out.bam>] [-o2 <R2.out.bam>]`
- **输入/输出**：默认输出同 `merge-two-end-bam1`。
- **依赖/注意**：依赖 samtools；匹配前无条件删除 read name 末尾 2 个字符。仅适用于确实以两字符 mate 后缀结尾的数据。

### 6.8 密码子比对与四倍简并位点

### `orthofinder-to-pal2nal`

- **用途**：验证蛋白/CDS 对应关系，运行 MAFFT 蛋白比对，再用 PAL2NAL 生成密码子比对。
- **调用**：`biohub run orthofinder-to-pal2nal -p <protein-groups-dir> -o <out-dir> -n <all.cds.fa>`
- **输入/输出**：每个非隐藏普通文件作为一个蛋白组；输出 `validation_summary.tsv`、`README.txt` 和每组的蛋白比对、CDS FASTA、PAML 密码子比对。
- **依赖/注意**：依赖 mafft、`pal2nal.pl`。CDS ID 必须唯一、存在、非空且长度为 3 的倍数；任何组失败或跳过使整体返回非零。Migut ID 有专用前三段映射兼容规则。

### `orthofiner-to-pal2nal`

- **用途**：`orthofinder-to-pal2nal` 的历史拼写兼容入口。
- **调用**：参数、输出和失败行为与规范命令完全相同。
- **兼容期**：保留到 v1.x；新脚本应使用正确拼写。

### `get-diff-sites-from-orthology`

- **用途**：扫描目录内密码子 FASTA，比对所有样本都属于四倍简并密码子的位点，并输出各密码子位置及拼接序列。
- **调用**：`biohub run get-diff-sites-from-orthology -i <codon-fasta-dir> [-o <out-dir>]`
- **输入/输出**：默认原地输出；生成 `allSingleCopyGeneSeq.fasta`、`TheFirstSite.fasta`、`TheSecondSite.fasta`、`TheThirdSite.fasta`、`fourDegenerateSite.fasta`、`fourDegenerateCodenFile.fasta`、`filesImportScripts.txt`。
- **依赖/注意**：样本名截取 FASTA 标题前 9 个字符；当前扫描步长基于首批长度，要求各序列同框同长。目录内所有普通文件都会尝试按 FASTA 读取，输出目录不要与输入目录混用。

### `compare-orthology`

- **用途**：`get-diff-sites-from-orthology` 的兼容别名。
- **调用**：`biohub run compare-orthology -i <codon-fasta-dir> [-o <out-dir>]`
- **输入/输出**：产物和默认原地写入行为与主命令相同。
- **依赖/注意**：新流程使用主 ID，并为输出指定独立目录。

### `get-four-degenerate-sites`

- **用途**：以每个比对首条序列为参考，提取参考为四倍简并密码子且该密码子在组内完全一致的第三位。
- **调用**：`biohub run get-four-degenerate-sites -i <codon-fasta-dir> [-o <out-dir>]`
- **输入/输出**：生成 `fourDegenerateSite.fasta` 和 `fourDegenerateSite.stat`，后者记录总提取位点数和参考密码子总数。
- **依赖/注意**：样本名同样截取标题前 9 字符；各组按最短序列长度扫描。该定义要求整个密码子完全一致，比通常“第三位可变的四倍简并位点”更严格，发表方法部分必须准确描述。

## 7. 直接命令组

这些入口早于统一目录，仍用于稳定夹具和历史工作流。它们不全部显示在 57 个命令目录中。

### 7.1 `rename`

```bash
biohub rename hjjn-genes -i <two-column.tsv> [-o <output>]
biohub rename scaffolds -i <table.tsv> -l <map.tsv> [-o <output>]
biohub rename fasta-scaffolds -i <input.fa> -l <map.tsv> [-o <output.fa>]
```

- `hjjn-genes` 处理前两列，把第二列最后一个 `gene` 后的纯数字补齐为 5 位，并保留点号扩展。
- `scaffolds` 的映射为 `新ID 旧ID`，只替换输入第一列，未映射保留。
- `fasta-scaffolds` 的映射为 `旧ID 新ID`，保留未映射记录，但 FASTA 标题只保留 ID，原描述被删除。

### 7.2 `blast reciprocal`

```bash
biohub blast reciprocal -i <A-to-B.tsv> -r <B-to-A.tsv> [-o <rbh.tsv>]
```

读取每行前两列；每个 query 保留首次命中，再输出互惠对。与目录命令 `compare-two-blast` 的“重复 query 保留最后一条”不同。

### 7.3 `gff`

```bash
biohub gff filter-ncbi --gff <input.gff3> [-o <output.gff3>]
biohub gff filter-gemoma --gff <input.gff3> [-o <output.gff3>]
biohub gff convert-ty1-hjjn --gff <input.gff3> --bed <placement.bed> [-o <output.gff3>]
```

- `filter-ncbi` 默认输出 `TA-filtered.gff3`，重建 gene/mRNA/exon/CDS ID，只处理按块排列的 gene、mRNA、CDS。
- `filter-gemoma` 默认输出 `gemoma-longest.gff3`，按累计 CDS 长度为每个 gene 保留一个转录本。
- `convert-ty1-hjjn` 默认输出 `Results.gff3`，按 BED-like 布局处理正反向坐标；同目录还会写 `<stem>-not-in-bed.txt`、`<stem>-errors-scaffolds.txt`、`change-annot.log` 和 split-scaffold 日志。建议空目录运行。

### 7.4 `fasta longest-transcript`

```bash
biohub fasta longest-transcript -f <input.fa> [-o <longest.fa>]
```

识别标题中 `gene=` 或 `gene_name=`；失败时使用第 4 个空白字段，再回退到首字段。按 gene 保留最长记录，等长保留首次，按 gene 首次出现顺序输出。这是当前推荐的通用最长转录本入口。

### 7.5 `stats`

```bash
biohub stats coverage-ratio -i <observed.tsv> -r <total.tsv> [-o <ratio.tsv>]
biohub stats wgcna-weight -i <header.tsv> [-o <identity.tsv>]
biohub stats hic-matrix-reindex -b <bins.bed> -m <matrix.tsv> -p <group-file-or-dir> [-o <reindexed.tsv>]
```

- `coverage-ratio` 使用两文件前两列按 key 连接，输出 `observed/total`，保留 6 位小数；total 为 0 时输出 0。
- `wgcna-weight` 只读取输入首行名称，生成带行名的单位矩阵，不读取后续权重数值。
- `hic-matrix-reindex` 用 BED 第 1、4 列建立 contig→位置；group 可为单文件或包含 `group*` 文件的目录，每行使用第 2 列 contig、第 3 列方向。方向 `1` 反转位置顺序，再替换矩阵前两列。

### 7.6 `psmc merge`

```bash
biohub psmc merge -d <input-dir> [-p <suffix>] [-o <merged.tsv>]
```

默认后缀 `.0.txt`，默认输出 `merge.psmc.0.txt`。文件按名称排序；样本名取文件名第一个点号前内容；输入前两列作为 Time、Ne，输出表头为 `Sample Time Ne`。

## 8. 验证状态

### 8.1 自动化覆盖

当前仓库测试验证：

- CLI 帮助、未知命令返回码、目录 JSON 后端与依赖元数据；
- HJJN gene 名规范化；
- scaffold 重命名的未映射保留；
- `fasta longest-transcript` 选择逻辑；
- PAL2NAL 前 CDS/蛋白 ID、缺失、重复和三联体长度检查；
- `annotation-vcf` 的 TSV、JSON、边界位置和变异类型；
- 本说明书对目录中全部脚本 ID 的章节覆盖；
- PAF 和 MUMmer coords 点阵图合成数据 smoke test（需 Rscript 环境）。
- 13 个 recipe 的 schema、目录契约、Python/R 语法和 Snakemake DAG dry-run；
- family de novo、功能富集、microbiome RDA、DESeq2 四个低风险 recipe 的合成数据实跑；
- 装配覆盖率完整分母、SyRI 第 11 列类型、MCMCTree 空白分隔 chain、PLINK2 hybrid
  结果与群体 pi 汇总的回归测试。

### 8.2 尚需真实数据验证的范围

以下类别已有实现，但仍应视为“迁移后待领域验证”：

- PASA、GeMoMa、NextGenomics 等来源的多种 GFF 属性和记录排序；
- 大规模 BAM 的内存、顺序、同名 reads 和 header 一致性；
- 历史 family/GO/BUSCO 表的非标准字段；
- LASTZ→JCVI 伪基因切分与历史结果逐项一致性；
- 多物种 orthology 与四倍简并位点的生物学定义；
- 真实超大基因组下的内存和运行时间。

发表前最低验证建议：

1. 保存精确命令、BioHub 版本、输入校验和、依赖版本和容器 digest。
2. 为每个关键命令选择去标识化小型代表数据及人工批准输出。
3. 对坐标转换检查正链、负链、首尾边界、未映射和 split scaffold。
4. 对同源关系检查并列 score、重复 query、无命中和物种数不齐。
5. 对密码子流程检查 reading frame、stop codon、gap、ID 映射和样本名截断。
6. 将结果与历史脚本或独立工具比较，并解释任何差异。

## 9. 故障排查

### 9.1 `output exists` 或 `File exists`

`biohub run` 默认保护主输出。确认旧文件可替换后传 `--force`，或改用新路径：

```bash
biohub run dotplot --input align.paf --output plot.pdf --force
```

多文件固定输出命令不一定受保护，仍建议新目录运行。

### 9.2 `R backend script not found`

从仓库根目录/`biohub-rs` 目录运行，使用完整发布归档，或设置：

```bash
export BIOHUB_R_DIR=/path/to/biohub/r
```

该目录必须包含 `dotplot.R` 和 `psmc_plot.R`。

### 9.3 `missing external dependency`

运行 `biohub doctor --strict`。确认程序名称准确位于 `PATH`；PAL2NAL 可执行文件必须名为 `pal2nal.pl`。容器镜像提供主要依赖。

### 9.4 输出为空但返回成功

常见原因：字段数不足、ID/染色体名不匹配、标题不符合历史规则、长度阈值过高、目录混入非目标文件。先用小样本检查分隔符和首几行：

```bash
sed -n '1,5p' input.tsv
biohub catalog --format json
```

再检查 stderr 和输出行数。不要把“空输出”自动解释为“没有生物学信号”。

### 9.5 输出顺序变化

部分实现使用哈希聚合，跨运行或平台不保证行顺序。若顺序不具生物学意义，在比较前使用稳定业务键排序；若顺序有意义，请提交 issue 并提供去标识化最小复现。

## 10. 可重复性、数据安全与引用

- 不提交私有基因组、未发表样本 ID、凭据或原始患者数据。
- 推荐保存 `biohub --version`、`biohub doctor --json`、命令行、标准错误、输入/输出 SHA256 和运行环境。
- 使用 `cargo build --locked` 或带 digest 的容器镜像。
- 发布工作应引用仓库 `CITATION.cff`；获得 Zenodo DOI 后优先引用 DOI 对应版本。
- BioHub 使用 MIT License；第三方版权与来源见 `NOTICE`。
- 问题报告应包含版本、精确命令、预期/实际结果和最小去标识化数据，见仓库 Issue 模板与 `SECURITY.md`。

## 11. 维护者文档

- `docs/RECIPES.zh-CN.md`：13 个科研 workflow 的详细使用、结果解释和发表复现检查；
- `docs/COMMANDS.md`：命令支持等级与晋级规则；
- `docs/SCRIPT_MIGRATION.md` 与矩阵：91 项脚本审计决策、验证状态和延期原因；
- `docs/RELEASE.md`：发布检查清单；
- `r/README.md`：R 后端契约；
- `CONTRIBUTING.md`：贡献与测试要求；
- `CHANGELOG.md`：版本变化；
- `docs/AI_USAGE.md`：AI 使用披露。

目录是机器可读事实源，本说明书是面向用户的行为契约。新增、删除或重命名脚本 ID 时，文档覆盖测试会要求同步更新本文件。
