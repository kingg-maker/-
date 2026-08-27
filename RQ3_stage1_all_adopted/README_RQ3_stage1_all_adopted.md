# RQ3 Stage 1 全部采纳版（对照 RQ3改 826.docx）

这一版按队友文档逐条全部采纳，先跑出一套“照文档执行”的基准结果，后续再与推荐修正版对比。

## 实现内容与文档条目对应

| # | 文档要求 | 本版实现 | 状态 |
|---|---|---|---|
| 1 | RQ3 按物质/技能/社会参与三维拆解 DVI 诊断异质性 | 三个维度分别建模，不再只用 DVI 总分 | 采纳 |
| 2 | 主分析样本仅保留三个维度均完整观测者 | 过滤 Dmat/Dskill/Dsocial 非缺失；默认慢性基期主样本（`USE_CHRONIC_ONLY=TRUE`） | 采纳 |
| 3 | 物质维度 = 无电脑、无手机均值（0/0.5/1） | `Dmat=round(dvi_mat,1)` | 采纳 |
| 4 | 技能维度 = 低教育、无大专学历子女均值（0/0.5/1） | `Dskill=round(dvi_skill,1)` | 采纳 |
| 5 | 社会参与维度 = DA056 1/2/4/5 + 6/8/9 代理 | `Dsocial=round(dvi_mot,1)`（与文档逐项一致） | 采纳 |
| 6 | 结果变量 BHCI 与功能健康 | `bhci_fixed`、`func_cap` | 采纳 |
| 7 | 九个维度×状态子样本分别跑 | 核心表逐格建模 | 采纳 |
| 8 | 按文档公式跑 TWFE 基准 | `treat + 控制变量 | ID_num + wave`，城市聚类 | 采纳 |
| 9 | 与 RQ1 一致的执行口径 | 同时输出堆叠 CS-DID（cohort×unit + cohort×wave，城市聚类） | 采纳 |
| 10 | 控制变量含“是否与子女同住” | `live_with_child`：优先 2011 家庭名册，缺失回退 `cb053==1/2`；该控制会带来约 49% 缺失 | 采纳 |
| 11 | 核心表：N/Treated N/Control N/ATT/SE/95%CI/p | `rq3_stage1_core.csv` | 采纳 |
| 12 | ATT(1)-ATT(0) 及 SE/CI/p | `rq3_stage1_contrasts.csv`（交互项估计，非手动相减） | 采纳 |
| 13 | Global heterogeneity test p-value | `rq3_stage1_global_test.csv`：每维度 2df 联合 Wald + 全维 6df 全局检验 | 采纳 |

## 运行方式

前置：R 已安装 `tidyverse`、`fixest`、`haven`。

```r
source("RQ3_stage1_all_adopted.R", encoding = "UTF-8")
```

## 输入

- `rq2_final_data.csv`（08-15 outputs 的 `rq2_final_v1`，已含 DVI 三维与全部控制变量）
- `CHARLS.csv`、`hhmember.dta`（构建“是否与子女同住”，路径在脚本顶部可改；缺失时回退 `family_information.dta` 的 `cb053==1/2`）

## 输出

- `rq3_stage1_core.csv`：九个子样本 × TWFE/CSDID × BHCI/功能健康
- `rq3_stage1_contrasts.csv`：每维度 ATT(1)-ATT(0)
- `rq3_stage1_global_test.csv`：每维度联合 Wald 与全局异质性检验
- `rq3_stage1_sample_cells.csv`：各档基期人数与人年 N

## 口径说明与已知代价（全部采纳带来的）

- TWFE 与 CS-DID 结果并存；文档第 8 条按原文 TWFE 执行，第 9 条按 RQ1 主口径 CS-DID 执行。
- “是否与子女同住”作为控制变量会丢弃名册未覆盖的样本（覆盖率约 51%），TWFE 的 N 会明显小于 CS-DID 的 N。
- `Treated_N_pw`/`Control_N_pw` 为人年行数；CSDID 中从未处理组跨两个队列重复计入，故 N 比真实人数大。
- 极小档（部分维度 state=0 或 0.5）可能估计失败，脚本会保留 NA 行并打印原因，不中断。
