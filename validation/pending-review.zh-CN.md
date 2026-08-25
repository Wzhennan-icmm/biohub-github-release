# BioHub v0.4 待人工批准清单

自动验证不能代替科研领域判断。正式 `v0.4.0` 发布前，以下14项必须由具名审核人
使用去标识化代表数据复核，并按 `review-template.md` 保存记录。

公开 CC0-1.0 合成数据、参数、参考基线和验收容差已分为四个验证包：

```bash
python3 tools/validation_review.py build --pack all
python3 tools/validation_review.py verify --pack all
```

每个包生成 `validation/evidence/<pack>/review.md`、输入/输出 SHA256、软件
版本、完整命令和自动差异。脚本不提供自动批准入口。

## GFF 与坐标转换

- `014 convert-gemoma-gff3`：用 GFF3 validator 检查 ID/Parent、feature 类型、phase 和层级。
- `015 convert-gene-annotation-contigs2chr-PASA`：复核正反链首尾坐标、1-based 闭区间和未映射记录。
- `016 convert-gene-annotation-scaffold2chr-nextgenomics`：复核 split scaffold、越界记录和辅助日志。

## 同源与密码子

- `030 get-best-hit-from-blast`：确认物种数、重复命中和互惠集合语义符合实际项目。
- `033 get-diff-sites-from-orthology`：与独立脚本对照各密码子位置及四倍简并位点。
- `034 get-four-degenerate-sites`：确认“所有样本密码子完全相同”严格规则符合论文方法。
- `043 comparative-orthology-codon`：用真实小型 orthogroup 检查蛋白/CDS配对、MAFFT和PAL2NAL结果。

## 可视化

- `045 plot-depth-pandepth`
- `046 plot-depth-pandepth2`
- `047 plot-mosdepth-point`

三项已有整文件 fingerprint 回归。审核人仍需查看标题、轴、范围、颜色、点位遮挡和
出版尺寸，并注明接受的渲染环境与容差。

## 统计 Recipe

- `049/050 microbiome-rda`：核对 Hellinger 转换、条件变量、置换检验和解释率。
- `053 functional-enrichment`：核对 foreground/background、超几何检验和 BH 校正。
- `072 rnaseq-deseq2`：核对设计矩阵、contrast、归一化、独立过滤和多重检验。

## 批准步骤

1. 固定 BioHub commit、输入 SHA256、参考软件版本和完整参数。
2. 运行 BioHub 与独立参考实现，保存输出 SHA256 和差异说明。
3. 审核人签署姓名、单位、日期与 `approved/rejected` 结论。
4. 将记录加入 `validation/`；把 `reviews.tsv` 状态改为 `approved`。
5. 同一提交移除迁移矩阵对应 `pending` 标记。
6. 运行 `python3 tools/validate_release.py --release --tag v0.4.0`；必须通过。

建议逐包批准口令：`批准 annotation-coordinates`、`批准 orthology-codon`、
`批准 visualization`、`批准 statistics`。每个口令只覆盖对应 inventory IDs；
未明确批准的包继续保持 `pending`。
