# BioHub v0.4 Recipe 使用与发表复现手册

> 文档版本：0.4.0<br>
> 最后更新：2026-08-14<br>
> 适用对象：运行 BioHub 多步骤科研流程、复核结果或准备论文方法与补充材料的分析人员

本手册覆盖 BioHub `0.4.0` 内置的 13 个 Snakemake recipe。命令级工具见
[BioHub v0.4 功能说明书](USER_GUIDE.zh-CN.md)。Recipe 状态均为
`experimental`：软件会检查配置、输入一致性、运行状态和可重复性元数据，但不会替研究者
决定样本设计、物种树、校准点、前景分支、过滤阈值或因果解释。

示例中的尖括号表示必须替换的占位值。表中范围来自提交的 JSON Schema；没有提供“推荐值”
时，应依据研究设计、数据质量和领域评审填写。不要把合成示例数值直接用于发表分析。

## 1. 通用运行模型

### 1.1 安装与发现

从源码构建核心程序：

```bash
git clone https://github.com/Wzhennan-icmm/biohub-github-release.git
cd biohub-github-release/biohub-rs
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"
biohub --version
biohub recipe list
```

先查看目标流程声明，不猜依赖：

```bash
biohub recipe describe selection-branch-site
biohub recipe describe selection-branch-site --format json
biohub doctor --recipe selection-branch-site --strict
```

领域容器按目录声明分组：

| 镜像 | Recipe |
| --- | --- |
| `ghcr.io/wzhennan-icmm/biohub-comparative:0.4.0` | orthology/codon、CAFE、branch-site、MCMCTree、SyRI |
| `ghcr.io/wzhennan-icmm/biohub-assembly:0.4.0` | T2T assembly 评估 |
| `ghcr.io/wzhennan-icmm/biohub-population:0.4.0` | SNP GWAS、群体选择 |
| `ghcr.io/wzhennan-icmm/biohub-variant:0.4.0` | family de novo rate |
| `ghcr.io/wzhennan-icmm/biohub-omics:0.4.0` | DESeq2、功能富集、microbiome RDA |

`kmer-gwas` 依赖已停止维护的 Python 2 和用户提供的 legacy 脚本，因此不提供官方领域镜像。
发表运行应保存实际镜像 digest；在容器内设置：

```bash
export BIOHUB_CONTAINER_DIGEST='sha256:<实际镜像摘要>'
```

非标准安装可设置 `BIOHUB_RECIPE_DIR` 指向包含 recipe 子目录和 `_lib/` 的目录。

### 1.2 初始化、填写与路径

`init` 要求目标目录不存在或为空，生成 `config.yaml`、`config.schema.yaml` 和
`README.txt`：

```bash
biohub recipe init <recipe-id> --workdir configs/<recipe-id>
```

编辑时保留模板中的全部键。`additionalProperties: false`，拼错或多余键会使 schema
校验失败；`REQUIRED` 与不允许 `null` 的模板值必须替换。由于 Snakemake 以 run 目录作为
工作目录，输入文件推荐写绝对路径。若使用相对路径，必须先确认它们能从目标 run 目录正确
解析；正式执行后再核对 `inputs.manifest.tsv` 中记录的最终路径和校验和。

### 1.3 预检、验证与正式运行

推荐顺序：

```bash
biohub doctor --recipe <recipe-id> --strict
biohub recipe validate <recipe-id> --config configs/<recipe-id>/config.yaml
biohub recipe run <recipe-id> \
  --config configs/<recipe-id>/config.yaml \
  --workdir runs/<recipe-id>-<run-label> \
  --cores <N>
biohub recipe report --workdir runs/<recipe-id>-<run-label>
```

`validate` 使用临时目录执行 Snakemake dry-run，不保留结果。需要检查实际 run 目录布局时，
使用 `recipe run ... --dry-run`。Dry-run 只验证 DAG、schema 和可见输入，不验证真实领域软件
能否完成，也不证明结果具备生物学有效性。

### 1.4 工作目录、覆盖与恢复

BioHub 管理以下 run 路径：

| 路径 | 含义 |
| --- | --- |
| `config.resolved.yaml`、`config.sha256`、`recipe.id` | 固定运行配置及身份 |
| `command.sh` | 实际 Snakemake 调用，可用于审计 |
| `run.json` | `running`、`failed`、`finalization_failed` 或 `complete` 状态 |
| `versions.tsv` | BioHub、Snakemake 和全部声明依赖的版本探针结果 |
| `provenance.json` | workflow 哈希、profile、容器提示和可选 runtime digest |
| `recipe.sources.sha256` | Snakefile、schema、脚本和共享 provenance helper 校验和 |
| `inputs.manifest.tsv` | 逻辑输入、解析路径、字节数和 SHA256 |
| `logs/`、`results/`、`report/` | 日志、主结果、图表和归档 |
| `checksums.sha256` | 稳定 bundle 文件校验和；故意排除可变的 `run.json` |

已有 `run.json` 的目录默认拒绝重跑。失败原因修复后，仅在 recipe ID 和配置 SHA256 完全
一致时使用：

```bash
biohub recipe run <recipe-id> --config <同一配置> \
  --workdir <原run目录> --cores <N> --resume
```

修改任何配置后新建 run 目录。`recipe report` 默认保护已有 `report/report.html`；确认替换
后使用 `--force`。不要手工修改 `config.resolved.yaml`、哈希或运行状态来绕过保护。

### 1.5 Local、Slurm 与自定义 profile

Local 为默认 profile：

```bash
biohub recipe run <recipe-id> --config <config> --workdir <run> \
  --profile local --cores 8
```

内置 Slurm profile 依赖 Snakemake 8.6+ 和 `snakemake-executor-plugin-slurm`，默认只提交一个
作业：

```bash
biohub recipe run <recipe-id> --config <config> --workdir <run> \
  --profile slurm --cores 8
```

生产集群应复制 profile，再填写获批准的账户、分区、内存、时间和并发策略；自定义 profile
直接传目录路径。BioHub 不自动推断集群资源。

### 1.6 通用发表复现清单

每次正式分析至少归档：

1. Git commit、BioHub 版本、recipe ID、`config.resolved.yaml` 和 `command.sh`。
2. `versions.tsv`、`provenance.json`、容器 digest、输入与 recipe source 校验和。
3. 全部日志、验证汇总、主结果、生成报告和 recipe 压缩包。
4. 样本/物种纳入排除规则、阈值来源、随机种子、重复数、背景集或前景分支定义。
5. 去标识化代表数据上的 golden output、独立工具对照和领域专家审核记录。

`complete` 仅表示工作流成功结束和 bundle 完成，不表示统计假设、研究设计或生物学结论已通过审查。

## `comparative-orthology-codon`

### 适用场景与边界

验证 orthogroup 蛋白与统一 CDS FASTA 的 ID、数量和读框，调用 MAFFT 与 PAL2NAL 生成
gap-free PAML 密码子比对。依赖 `snakemake`、`python3`、`mafft`、`pal2nal.pl`。它不负责
推断 orthogroup，也不自动修复错误 CDS、内部终止密码子或物种命名。

### 输入准备

- `protein_groups_dir`：每个普通文件代表一个 orthogroup 的蛋白 FASTA；不要混入日志或结果。
- `cds_fasta`：所有蛋白 ID 对应的 CDS；ID 必须唯一，CDS 非空且长度为 3 的倍数。
- 每组应包含 `expected_taxa` 条蛋白和同数 CDS。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `protein_groups_dir` | 非空目录路径 |
| `cds_fasta` | 非空 CDS FASTA 路径 |
| `expected_taxa` | 整数，至少 2；每个 orthogroup 预期 taxa 数 |
| `biohub_executable` | 可执行 BioHub 名称或路径，通常为 `biohub` |

### 运行步骤

```bash
biohub recipe init comparative-orthology-codon --workdir configs/orthology-codon
biohub doctor --recipe comparative-orthology-codon --strict
biohub recipe validate comparative-orthology-codon --config configs/orthology-codon/config.yaml
biohub recipe run comparative-orthology-codon \
  --config configs/orthology-codon/config.yaml \
  --workdir runs/orthology-codon-01 --cores 8
```

### 产物与解读

- `results/codon_alignments/`：蛋白比对、选中 CDS、PAML 密码子比对及原始
  `validation_summary.tsv`。
- `results/summary.tsv`：总组数、完成数、失败/跳过数和 taxa 数不匹配数。
- `logs/validation.log`：逐组失败原因；任何失败或 taxa 数不符使流程失败。

### 失败、恢复与发表核对

先处理重复/缺失 ID、CDS 非三联体长度、组内数量不齐和 MAFFT/PAL2NAL 错误，再用原配置
`--resume`。发表时记录 orthogroup 来源、蛋白过滤规则、预期物种数、MAFFT/PAL2NAL 参数，
并验证每个最终比对物种齐全、长度同框、PAML header 正确。

## `gene-family-cafe`

### 适用场景与边界

准备基因家族 count matrix 和超度量二叉树，运行多个 CAFE5 replicate，按似然与参数稳定性
选择收敛簇，再提取 expansion/contraction。依赖 `snakemake`、`python3`、`cafe5`。流程不会
替用户修正物种树、估计合理误差模型或判断扩张收缩的功能意义。

### 输入准备

- `gene_counts`：`orthofinder` 格式首列必须为 `Orthogroup`；`cafe` 格式需含
  `Family ID` 或 `FamilyID`。至少两个唯一物种列，count 为非负整数。
- `tree`：带正分支长度的 Newick 二叉树；叶名称必须与 matrix 物种列完全一致。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `gene_counts`、`tree` | 非空输入路径 |
| `input_format` | `orthofinder` 或 `cafe` |
| `max_family_size` | 正整数；任一物种超过该 count 的家族被排除 |
| `ultrametric_tolerance` | 非负数；允许 root-to-tip 最大偏差 |
| `cafe_executable` | CAFE5 可执行名称或路径 |
| `cores`、`replicates`、`model_k` | 正整数 |
| `root_distribution` | 正数，对应 CAFE `-p` |
| `error_model` | 非负数、非空模型路径/字符串或 `null` |
| `lambda_value` | 正数或 `null`；固定 lambda 时填写 |
| `alpha_value` | 正数或 `null`；`model_k > 1` 时必须填写 |
| `family_pvalue` | 0–1，对应 CAFE `-P`；来源必须在方法中说明 |
| `minimum_converged_replicates` | 正整数；收敛簇最低 replicate 数 |
| `likelihood_tolerance` | 非负数；收敛簇允许的 `-lnL` 范围 |
| `parameter_cv_tolerance` | 非负数；lambda/alpha 变异系数阈值 |

### 运行步骤

```bash
biohub recipe init gene-family-cafe --workdir configs/cafe
biohub doctor --recipe gene-family-cafe --strict
biohub recipe validate gene-family-cafe --config configs/cafe/config.yaml
biohub recipe run gene-family-cafe --config configs/cafe/config.yaml \
  --workdir runs/cafe-01 --cores 8
```

### 产物与解读

- `results/family_filter_manifest.tsv`：每个家族 included/excluded、原因、最大和总 count。
- `results/tree_qc.tsv`、`results/input_summary.tsv`：超度量深度和过滤统计。
- `results/model_fit_replicates.tsv`：各 replicate `-lnL`、lambda、alpha、数值警告及收敛簇标记。
- `results/selected_run.txt`：通过最低收敛数后选择的 run。
- `results/expanded_families.tsv`、`contracted_families.tsv`、`summary_by_node.tsv`：非零 branch
  change 的机械分类，不等于显著功能改变。

### 失败、恢复与发表核对

常见失败：taxa 不一致、非二叉/非超度量树、全部家族被过滤、CAFE 数值警告、收敛 replicate
不足。不要修改配置后 resume；调整模型必须新建 run。发表时报告 count 来源、过滤规则、树及
时间单位、误差模型、lambda/alpha、replicate 数、收敛判据和 branch change 定义。

## `selection-branch-site`

### 适用场景与边界

为每个 test 运行 codeml branch-site Model A alternative/null，计算 50:50 mixture LRT p 值，
对 manifest 中全部成功 test 做全局 BH，并仅从 alternative model 提取 BEB 位点。依赖
`snakemake`、`python3`、`codeml`。流程不选择前景分支，也不证明正选择结论。

### 输入准备

`tests_manifest` 必须为制表文件：

```text
test_id	alignment	marked_tree	foreground
<安全唯一ID>	<多序列PAML密码子比对>	<仅含一个#1的Newick>	<前景标签>
```

比对首行 taxa 数必须等于 `expected_taxa`，位点数必须为正且是 3 的倍数。树必须括号平衡、
分号结尾且恰有一个 `#1`。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `tests_manifest` | 非空 manifest 路径 |
| `expected_taxa` | 整数，至少 2 |
| `codeml_executable` | codeml 可执行名称或路径 |
| `codon_frequency` | 0–3，对应 codeml `CodonFreq` |
| `initial_kappa`、`alternative_initial_omega` | 正数；模型初值，不是结果 |
| `clean_data`、`optimizer_method` | 仅允许 0 或 1 |
| `beb_thresholds` | 至少一个唯一 0–1 阈值 |
| `bh_family` | 固定为 `all_manifest_tests` |

### 运行步骤

```bash
biohub recipe init selection-branch-site --workdir configs/branch-site
biohub doctor --recipe selection-branch-site --strict
biohub recipe validate selection-branch-site --config configs/branch-site/config.yaml
biohub recipe run selection-branch-site --config configs/branch-site/config.yaml \
  --workdir runs/branch-site-01 --cores 8
```

### 产物与解读

- `results/tests.normalized.tsv` 和 `input_summary.tsv`：规范路径、taxa 和 codon sites 审计。
- `results/runs/<test_id>/{alternative,null}/`：复制的输入、`codeml.ctl` 和 `mlc`。
- `results/branch_site_lrt.tsv`：两模型 lnL、LRT、mixture p、`bh_qvalue_global`。
- `results/beb_sites.tsv`：test、前景、位置、氨基酸、posterior、阈值和 PAML stars。
- `logs/validation.log`：未完成或解析失败 test；存在失败时流程非零。

### 失败、恢复与发表核对

排查未标记/多标记树、taxa 不齐、比对非 codon、codeml 非零或 `mlc` 缺少 lnL。发表时明确
前景分支、固定树、模型设置、mixture 分布、BH family、显著阈值和 BEB 阈值；保留每个
alternative/null `mlc`。不要把 BEB 位点或 calibration-bound node 当作独立验证。

## `dating-mcmctree`

### 适用场景与边界

为多个 MCMCTree replicate 执行 `usedata=3` 到 `usedata=2 in.BV` 两阶段流程，解析
`mcmc.txt` 年龄列并计算均值、分位数和简化 ESS。依赖 `snakemake`、`python3`、
`mcmctree`。流程不生成校准、不决定时间单位，也不替代 Tracer 等独立收敛诊断。

### 输入准备

`runs_manifest`：

```text
run_id	stage1_ctl	stage2_ctl
<唯一ID>	<usedata=3控制文件>	<usedata=2控制文件>
```

两个 ctl 必须含 `seqfile`、`treefile`、`outfile`、`ndata`、`usedata`；`ndata` 必须等于
`expected_loci`。控制文件相对输入路径按 ctl 所在目录解析。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `runs_manifest` | 非空 replicate manifest |
| `mcmctree_executable` | MCMCTree 可执行名称或路径 |
| `expected_loci`、`expected_internal_nodes` | 正整数 |
| `age_column_regex` | 非空正则，用于选择 `mcmc.txt` 年龄列 |
| `burnin_samples` | 非负整数；按样本行丢弃 |
| `minimum_ess` | 正数；每个 run/node 的最低 ESS |

### 运行步骤

```bash
biohub recipe init dating-mcmctree --workdir configs/mcmctree
biohub doctor --recipe dating-mcmctree --strict
biohub recipe validate dating-mcmctree --config configs/mcmctree/config.yaml
biohub recipe run dating-mcmctree --config configs/mcmctree/config.yaml \
  --workdir runs/mcmctree-01 --cores 4
```

### 产物与解读

- `results/runs/<run_id>/`：重写后的两阶段 ctl、`out.BV`/`in.BV`、`mcmctree.out`、
  `mcmc.txt`。
- `results/node_age_summary.tsv`：run、node、样本数、mean、SD、median、min/max、95% 分位
  区间、ESS 和 pass 标记。
- `results/summary.tsv`：请求/失败 run、node-run 记录数和低 ESS 数。
- `logs/<run_id>.stage{1,2}.log` 与 `logs/validation.log`：运行和链解析错误。

### 失败、恢复与发表核对

常见失败：`ndata` 不符、usedata 顺序错误、stage1 无 `out.BV`、stage2 无 `mcmc.txt`、年龄
正则错配、burn-in 后样本不足、node 数不符或 ESS 未达标。发表时记录序列组织方式、所有校准
上下界、时间单位、ctl 参数、burn-in、采样数、重复策略和独立收敛检查。

## `assembly-t2t-evaluate`

### 适用场景与边界

统计 assembly FASTA、端粒 motif 末端证据，并用 minimap2 对 reference 生成 PAF 覆盖和
identity 汇总。依赖 `snakemake`、`python3`、`minimap2`。它不宣称仅凭序列数或 motif 即可
判定 T2T，也不替代结构准确性和人工染色体审核。

### 输入准备

`assemblies_manifest`：

```text
assembly_id	fasta
<唯一ID>	<assembly FASTA>
```

FASTA ID 必须唯一、序列非空。`reference_fasta` 提供 alignment denominator；PAF 中所有名称、
长度和坐标必须与两份 FASTA 一致。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `assemblies_manifest`、`reference_fasta` | 非空路径 |
| `expected_chromosomes` | 正整数；用于 sequence count 检查 |
| `minimum_contig_length` | 正整数；long sequence 计数阈值 |
| `telomere_motif` | 仅 A/C/G/T/N；同时搜索反向互补 |
| `telomere_window_bp`、`minimum_telomere_hits` | 正整数；末端窗口和最低命中数 |
| `minimap2_executable` | minimap2 名称或路径 |
| `minimap_preset` | `asm5`、`asm10` 或 `asm20` |
| `threads` | 正整数 |

### 运行步骤

```bash
biohub recipe init assembly-t2t-evaluate --workdir configs/t2t
biohub doctor --recipe assembly-t2t-evaluate --strict
biohub recipe validate assembly-t2t-evaluate --config configs/t2t/config.yaml
biohub recipe run assembly-t2t-evaluate --config configs/t2t/config.yaml \
  --workdir runs/t2t-01 --cores 8
```

### 产物与解读

- `results/stats/<assembly_id>.sequences.tsv`：每条序列长度、N、长度 pass、左右 motif 命中和
  双端标记。
- `results/paf/<assembly_id>.paf` 与 `paf_summary/`：alignment、完整 FASTA query/target
  denominator 的 union coverage、加权 identity。
- `results/assembly_summary.tsv`：sequence count、total bp、N50、long sequence、双端端粒、N、
  coverage/identity 和 expected chromosome 检查。

### 失败、恢复与发表核对

空/重复 FASTA ID、PAF 长度或坐标不一致、无 alignment 会失败；sequence count 不符写 warning，
不自动宣告失败。发表时记录 reference、minimap2 preset、motif、窗口、命中阈值和完整 denominator，
并结合 BUSCO、QV、结构一致性、缺口及人工染色体端点证据。

## `synteny-sv`

### 适用场景与边界

对 assembly pair 运行 minimap2、可选 PAF 过滤及 SyRI，按官方输出第 11 列 annotation type
汇总记录数。依赖 `snakemake`、`python3`、`minimap2`、`syri`。结果是调用记录汇总，不是独立
验证过的真值 SV 集。

### 输入准备

`pairs_manifest`：

```text
pair_id	reference_fasta	query_fasta
<唯一ID>	<reference FASTA>	<query FASTA>
```

两份 FASTA 的 ID 必须唯一且序列非空。启用 `require_matching_sequence_ids` 时，ID 集必须相同；
顺序差异会记录但不失败。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `pairs_manifest` | 非空 pair manifest |
| `require_matching_sequence_ids` | 布尔；是否强制 reference/query ID 集一致 |
| `minimap2_executable`、`syri_executable` | 可执行名称或路径 |
| `minimap_preset` | `asm5`、`asm10` 或 `asm20` |
| `threads`、`minimum_alignment_length`、`syri_cores` | 正整数 |
| `minimum_mapq` | 0–255 |
| `syri_use_filtered_paf` | 布尔；SyRI 使用 raw 或过滤 PAF |
| `syri_include_cigar`、`syri_include_snps` | 布尔；控制 SyRI `--cigar` 与 SNP 输出 |

### 运行步骤

```bash
biohub recipe init synteny-sv --workdir configs/syri
biohub doctor --recipe synteny-sv --strict
biohub recipe validate synteny-sv --config configs/syri/config.yaml
biohub recipe run synteny-sv --config configs/syri/config.yaml \
  --workdir runs/syri-01 --cores 8
```

### 产物与解读

- `results/qc/<pair_id>.tsv`：序列数、ID 集与顺序检查。
- `results/paf/<pair_id>.paf` 和 `.filtered.paf`：raw 与按 block length/MAPQ 筛选结果。
- `results/syri/<pair_id>/pair.syri.out`：SyRI 正式输出。
- `results/syri_summary.tsv`：pair、有效记录数、按第 11 列统计的 type/count。

### 失败、恢复与发表核对

ID 集不符、PAF 全部被过滤、SyRI 成功但无非空输出、少于 11 列或空 type 会失败。发表时记录
assembly 方向、minimap2 preset、是否使用 filtered PAF、长度/MAPQ 阈值、CIGAR/SNP 选项和
SyRI 版本；使用独立比对、read support 或人工检查验证关键 SV。

## `population-gwas`

### 适用场景与边界

验证 VCF/PGEN/BED 与 phenotype/covariate 样本，按 trait 调用 PLINK2 `--glm`，汇总 ADD test
有效 p 值和每个 trait 的机械最小 p。依赖 `snakemake`、`python3`、`plink2`。流程不提供群体
结构建模策略、阈值校正选择、位点独立验证或因果推断。

### 输入准备

`traits_manifest`：

```text
trait_id	phenotype_file	phenotype_column
<唯一ID>	<含IID列的TSV>	<表型列名>
```

VCF 从 `#CHROM` 后读取样本；PGEN prefix 需要 `.pgen/.pvar/.psam`；BED prefix 需要
`.bed/.bim/.fam`。phenotype 和 covariate 表必须含唯一非空 `IID`。`require_complete_samples`
决定是否允许 genotype 中存在未进入 phenotype 的样本。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `genotype_kind` | `vcf`、`pfile` 或 `bfile` |
| `genotype`、`traits_manifest` | 非空路径；pfile/bfile 使用 prefix |
| `covariates_file` | 含 `IID` 的 TSV 或 `null` |
| `covariate_columns` | 唯一列名数组，可为空 |
| `require_complete_samples` | 布尔；要求 genotype/phenotype 样本完整一致 |
| `plink2_executable` | PLINK2 名称或路径 |
| `minor_allele_frequency` | 0–0.5，对应 `--maf` |
| `maximum_variant_missingness`、`maximum_sample_missingness` | 0–1，对应 `--geno`、`--mind` |
| `hardy_weinberg_pvalue`、`significance_threshold` | 大于 0 且不超过 1；来源必须说明 |

### 运行步骤

```bash
biohub recipe init population-gwas --workdir configs/gwas
biohub doctor --recipe population-gwas --strict
biohub recipe validate population-gwas --config configs/gwas/config.yaml
biohub recipe run population-gwas --config configs/gwas/config.yaml \
  --workdir runs/gwas-01 --cores 8
```

### 产物与解读

- `results/traits.normalized.tsv`、`input_summary.tsv`：trait 路径和样本审计。
- `results/gwas/<trait_id>/`：PLINK2 `.glm.*`，包括 logistic hybrid 文件。
- `results/trait_status.tsv`：结果文件数、有效 ADD rows、阈值内 rows。
- `results/lead_associations.tsv`：每 trait 最小有限 p 记录；`passes_threshold` 仅按配置机械判断。

### 失败、恢复与发表核对

缺少 genotype 组件、重复 IID、样本集合不符、表型列缺失、无 `.glm` 或无有限 ADD p 会失败。
发表时记录表型编码、样本排除、QC 阈值、协变量、群体结构、模型类型、multiple-testing 方法和
独立验证；不要把 lead table 当作校正后显著或因果证据。

## `population-selection`

### 适用场景与边界

按 comparison 调用 VCFtools，计算 windowed Weir FST 和两群体 windowed nucleotide diversity，
按配置阈值提取候选 FST window。依赖 `snakemake`、`python3`、`vcftools`。它不定义群体、
中性模型、重组率或经验显著性。

### 输入准备

`comparisons_manifest`：

```text
comparison_id	population1_samples	population2_samples
<唯一ID>	<每行一个VCF样本ID的文件>	<另一不重叠样本文件>
```

VCF 支持未压缩、`.gz` 和 `.bgz`。样本 ID 必须唯一；每个 comparison 的两群体不得重叠，且
全部样本必须存在于 VCF。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `vcf`、`comparisons_manifest` | 非空路径 |
| `vcftools_executable` | VCFtools 名称或路径 |
| `window_size_bp`、`window_step_bp` | 正整数 |
| `minor_allele_frequency` | 0–0.5 |
| `minimum_site_call_rate` | 0–1，对应 VCFtools `--max-missing` |
| `candidate_fst_threshold` | 数值；候选 window 的机械筛选阈值 |

### 运行步骤

```bash
biohub recipe init population-selection --workdir configs/selection
biohub doctor --recipe population-selection --strict
biohub recipe validate population-selection --config configs/selection/config.yaml
biohub recipe run population-selection --config configs/selection/config.yaml \
  --workdir runs/selection-01 --cores 8
```

### 产物与解读

- `results/comparisons.normalized.tsv`、`input_summary.tsv`：绝对样本文件路径与群体大小。
- `results/comparisons/<comparison_id>/`：原始 `.windowed.weir.fst` 和两个 `.windowed.pi`。
- `results/candidate_fst_windows.tsv`：达到配置阈值的 chromosome/window/FST。
- `results/comparison_summary.tsv`：有效和候选 window 数。
- `results/nucleotide_diversity_summary.tsv`：每群体有限 pi window 数、均值和范围。

### 失败、恢复与发表核对

重复/缺失/重叠样本、VCF 无 header、VCFtools 失败、pi 文件缺失或负/非有限 pi 会失败。
发表时记录样本量、群体定义、VCF QC、窗口/步长、MAF、call rate、阈值来源，并结合多样性、
单倍型、重组和人口史模型交叉验证候选区域。

## `kmer-gwas`

### 适用场景与边界

调用用户提供的 legacy Python 2 k-mer GWAS 脚本，按 trait 收集其
`pass_threshold_5per|10per` 结果并提取可识别 DNA k-mer。依赖 `snakemake`、`python3`、
`python2` 和外部脚本。无官方容器；接口强依赖目标 legacy 脚本行为。

### 输入准备

`traits_manifest`：

```text
trait_id	phenotype_file
<唯一ID>	<legacy脚本接受的表型文件>
```

`kmers_table_prefix` 必须与外部脚本预生成 table prefix 一致。BioHub 不验证 phenotype 和
k-mer table 样本集合；正式运行前应单独审计。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `traits_manifest`、`kmers_table_prefix`、`kmers_gwas_script` | 非空路径/prefix |
| `python2_executable` | Python 2 可执行名称或隔离环境路径 |
| `kmer_length`、`threads`、`minimum_data_points`、`minor_allele_count`、`kmers_number`、`permutations` | 正整数 |
| `minor_allele_frequency` | 0–0.5 |
| `permutation_tail_percent` | 仅 5 或 10；决定读取哪个 threshold 文件 |

### 运行步骤

```bash
biohub recipe init kmer-gwas --workdir configs/kmer-gwas
biohub doctor --recipe kmer-gwas --strict
biohub recipe validate kmer-gwas --config configs/kmer-gwas/config.yaml
biohub recipe run kmer-gwas --config configs/kmer-gwas/config.yaml \
  --workdir runs/kmer-gwas-01 --cores 8
```

### 产物与解读

- `results/associations/<trait_id>/`：外部脚本原始目录。
- `results/trait_status.tsv`：每 trait significant k-mer 行数和源文件。
- `results/significant_kmers.tsv`：保留 trait、源行号和完整 raw record。
- `results/significant_kmers.fasta`：每行中首个仅含 A/C/G/T/N 字段；无法识别者只保留 TSV。

### 失败、恢复与发表核对

Python 2/脚本不可用、外部返回非零或缺少预期 threshold 文件会失败。必须固定 legacy script
checksum、Python 2 环境、k-mer table 生成流程、样本映射、permutation、MAF/MAC 和阈值；
候选 k-mer 需要回贴参考基因组、控制群体结构并独立验证。Python 2 环境应隔离，禁止联网安装
不受信代码。

## `family-denovo-rate`

### 适用场景与边界

审计 family/pair candidate、callable BED、tier 和 evidence class，以 callable bp × ploidy 为
机会数，报告 pair 和 combined de novo rate 及 Garwood exact Poisson 区间。依赖
`snakemake`、`python3` 标准库。它不执行 variant calling、亲缘验证、污染检测或父母来源判定。

### 输入准备

`pairs_manifest`：

```text
family_id	pair_id	candidates_tsv	callable_bed
<family>	<唯一pair>	<candidate表>	<0-based half-open BED>
```

candidate 表必须含：

```text
chrom	position	ref	alt	tier	evidence_class
```

`position` 为 1-based；同一 candidate 文件内 `chrom/position/ref/alt` 必须唯一。callable BED
区间必须有效，重叠区间自动合并。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `pairs_manifest` | 非空 manifest |
| `included_tiers` | 至少一个唯一非空 tier |
| `included_evidence_classes` | `candidate` 和/或 `experimentally_confirmed` |
| `ploidy` | 正整数；机会数乘数 |
| `confidence_level` | 0–1 开区间 |
| `require_all_candidates_callable` | 布尔；任一不可调用 candidate 是否使流程失败 |

### 运行步骤

```bash
biohub recipe init family-denovo-rate --workdir configs/denovo
biohub doctor --recipe family-denovo-rate --strict
biohub recipe validate family-denovo-rate --config configs/denovo/config.yaml
biohub recipe run family-denovo-rate --config configs/denovo/config.yaml \
  --workdir runs/denovo-01 --cores 2
```

### 产物与解读

- `results/candidate_audit.tsv`：callable、tier/class 纳入状态和是否进入 numerator。
- `results/pair_rates.tsv`：candidate count、callable bp、ploidy、机会数、rate 和置信区间。
- `results/combined_rate.tsv`：所有 pair 合并 numerator/opportunity 和 Poisson 区间。
- `logs/validation.log`：逐个不可调用 candidate；strict callable 模式下导致失败。

### 失败、恢复与发表核对

空 BED、非法坐标、重复 candidate、字段缺失或 strict callable 违规会失败。发表时记录 callable
定义、过滤和验证 tier、每个 pair 机会数、ploidy、代数假设、置信水平及 combined 策略；零
candidate 仍需报告上限。使用独立家系/测序 QC 验证每个纳入事件。

## `rnaseq-deseq2`

### 适用场景与边界

读取 raw integer count、sample design 和 contrasts，过滤低 count gene，验证 design 满秩，运行
DESeq2 并输出 normalized count、差异表和 PCA。依赖 `snakemake`、`python3`、`Rscript`、
`DESeq2`。不接受 TPM/FPKM 作为 count，也不替代批次设计和样本 QC。

### 输入准备

- `counts_matrix`：首列必须为 `gene_id`，至少两个样本列；gene 唯一，值为有限非负整数。
- `samples_tsv`：含唯一 `sample_id` 及 design 使用的全部列；样本集与 count 列完全一致。
- `contrasts_manifest`：

```text
contrast_id	factor	numerator	denominator
<唯一安全ID>	<design中的分类变量>	<水平A>	<水平B>
```

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `counts_matrix`、`samples_tsv`、`contrasts_manifest` | 非空路径 |
| `design` | 以 `~` 开头的 R formula；变量必须存在 |
| `minimum_count` | 非负整数 |
| `minimum_samples` | 正整数；达到 minimum count 的最低样本数 |
| `alpha` | 大于 0 且不超过 1 |
| `minimum_absolute_log2_fold_change` | 非负数；significant 标记阈值 |
| `independent_filtering` | 布尔，传给 DESeq2 results |

### 运行步骤

```bash
biohub recipe init rnaseq-deseq2 --workdir configs/deseq2
biohub doctor --recipe rnaseq-deseq2 --strict
biohub recipe validate rnaseq-deseq2 --config configs/deseq2/config.yaml
biohub recipe run rnaseq-deseq2 --config configs/deseq2/config.yaml \
  --workdir runs/deseq2-01 --cores 4
```

### 产物与解读

- `results/differential_expression.tsv`：contrast、gene、baseMean、log2FC、SE、stat、p、padj、
  significant。
- `results/contrast_summary.tsv`：tested/adjusted/significant/up/down gene 数。
- `results/normalized_counts.tsv`：DESeq2 size-factor normalized count。
- `results/sample_qc.tsv` 和 `report/sample_pca.pdf`：library size、detected gene 和 PCA。
- `logs/deseq2.log`：DESeq2 版本、过滤后规模、design、alpha 和 LFC 阈值。

### 失败、恢复与发表核对

样本不一致、非整数 count、design 缺列/缺失/非满秩、contrast factor 非分类或水平不存在、过滤
后无 gene 会失败。发表时记录 count 生成方法、design formula、reference level、contrast 方向、
过滤阈值、alpha、LFC 条件、independent filtering 和 DESeq2 版本；PCA 只作 QC，不替代批次诊断。

## `functional-enrichment`

### 适用场景与边界

对多个 foreground set 和显式 background/association 做单侧超几何 over-representation test，
按配置 scope 执行 BH。依赖 `snakemake`、`python3`、`Rscript`，统计部分使用 base R。流程不
生成功能注释，不允许用“只注释 foreground”代替背景。

### 输入准备

- `foreground_tsv`：`set_id<TAB>gene_id`；set ID 安全，gene 非空。
- `background_genes`：无 header，每行一个可被选择的 background gene。
- `associations_tsv`：`gene_id`、`term_id`、`source`，可选 `term_name`；同一 source/term 的名称
  不得冲突。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `foreground_tsv`、`background_genes`、`associations_tsv` | 非空路径 |
| `sources` | 至少一个唯一非空 source，如项目实际使用的 GO/KEGG 标签 |
| `minimum_term_size`、`maximum_term_size`、`minimum_overlap` | 正整数；minimum 不得大于 maximum |
| `fdr` | 大于 0 且不超过 1 |
| `adjustment_scope` | `set_source`、`set` 或 `global` |
| `require_foreground_in_background` | 布尔；严格要求 foreground 属于 background |
| `plot_top_terms` | 正整数；每个 set 图中最多显示数 |

### 运行步骤

```bash
biohub recipe init functional-enrichment --workdir configs/enrichment
biohub doctor --recipe functional-enrichment --strict
biohub recipe validate functional-enrichment --config configs/enrichment/config.yaml
biohub recipe run functional-enrichment --config configs/enrichment/config.yaml \
  --workdir runs/enrichment-01 --cores 2
```

### 产物与解读

- `results/enrichment_all.tsv`：set/source/term、背景与前景大小、overlap、fold enrichment、p、
  padj、significant 和 overlap genes。
- `results/enrichment_significant.tsv`：`padj <= fdr` 的机械子集。
- `results/set_summary.tsv`：输入、背景内、背景外、已注释 gene 及 tested/significant term 数。
- `report/enrichment.pdf`：每个 set 按 padj 排序的 top term；空结果也会生成说明页。

### 失败、恢复与发表核对

背景为空、foreground 越界、source 无背景 association、term name 冲突或参数不合法会失败。
发表时明确 foreground、可检验 background、注释数据库/版本、source、term size、overlap、BH
scope 和 FDR；同时报告覆盖率与完整结果，不只展示显著条目。

## `microbiome-rda`

### 适用场景与边界

过滤 feature、转换 abundance，按显式 constraint 和可选 `Condition()` covariate 运行 vegan
RDA，输出 overall/term/axis permutation tests 与 scores。依赖 `snakemake`、`python3`、
`Rscript`、`vegan`。流程不决定 compositional 方法是否合适，也不作因果解释。

### 输入准备

- `feature_table`：首列必须为 `feature_id`，至少两个样本列；feature 唯一，值有限且非负。
- `metadata_tsv`：含唯一 `sample_id`、全部 constraint 和 condition 列；样本集必须完全一致。
- 字符/逻辑变量转为 factor；数值变量保持连续。

### 配置字段

| 字段 | 约束与含义 |
| --- | --- |
| `feature_table`、`metadata_tsv` | 非空路径 |
| `constraints` | 至少一个唯一合法列名 |
| `condition_covariates` | 唯一列名数组，可为空；不得与 constraints 重叠 |
| `transform` | `hellinger`、`relative`、`log1p` 或 `none` |
| `minimum_prevalence` | 0–1 |
| `minimum_total_abundance` | 非负数 |
| `drop_incomplete_samples` | 布尔；否则任一模型缺失值使流程失败 |
| `permutations` | 正整数 |
| `random_seed` | 非负整数 |
| `scaling` | 1 或 2 |

### 运行步骤

```bash
biohub recipe init microbiome-rda --workdir configs/rda
biohub doctor --recipe microbiome-rda --strict
biohub recipe validate microbiome-rda --config configs/rda/config.yaml
biohub recipe run microbiome-rda --config configs/rda/config.yaml \
  --workdir runs/rda-01 --cores 2
```

### 产物与解读

- `results/features_audit.tsv`：prevalence、total abundance 和 retained 标记。
- `results/model_summary.tsv`：输入/分析样本、输入/保留 feature、总/约束 inertia、比例和轴数。
- `results/permutation_overall.tsv`、`permutation_terms.tsv`、`permutation_axes.tsv`：固定随机种子的
  permutation test。
- `results/site_scores.tsv`、`feature_scores.tsv`、`biplot_scores.tsv`：RDA1/RDA2 scores；rank<2
  时第二轴可为 NA。
- `report/rda.pdf` 与 `logs/rda.log`：ordination 图、vegan 版本、formula、转换、seed 和 eigenvalues。

### 失败、恢复与发表核对

样本不一致、负/非有限 abundance、constraint 缺列或重叠、缺失数据策略不符、过滤后少于两个
feature、有效样本少于三个、零 sample total 或无 constrained axis 会失败。发表时记录 feature
定义、过滤、转换、constraint/Condition formula、缺失样本处理、permutation、seed、scaling 和
vegan 版本；使用诊断和敏感性分析评估 compositional 与混杂影响。

## 结果交付与论文材料

正式交付目录建议保留 BioHub 原始 run，不抽取后删除 provenance。向合作者提供时可额外建立：

```text
delivery/
├── README.txt
├── metadata/       # config、版本、命令、校验和、研究设计说明
├── tables/         # 审核后的主表和字段说明
├── figures/        # 审核后的图及图注
├── logs/           # 原始验证和运行日志
└── <recipe-id>.tar.gz
```

Methods 至少陈述：BioHub 版本/commit、recipe ID、输入来源和过滤、外部软件版本、全部影响统计
结果的参数、重复与随机种子、显著性与多重校正范围、失败/排除数量。Results 和 Discussion 才解释
生物学含义；不要把工作流成功状态写成科学验证。

## 问题报告

提交 GitHub issue 时提供：

- `biohub --version` 与 `biohub recipe describe <id> --format json`；
- 去标识化后的 `config.resolved.yaml`、`run.json`、相关日志和最小输入；
- 精确预期、实际结果和返回码；
- 是否容器运行、镜像 digest、操作系统和架构；
- 禁止上传私有基因组、未发表样本 ID、凭据或受控数据。
