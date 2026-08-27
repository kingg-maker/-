# RQ3 Stage 1 推荐版结果说明（2026-08-27）

主口径：慢性基期 + 2011-2018 + 三维完整样本（4922 人，16662 人年）；堆叠 CS-DID 为主，TWFE 仅作稳健性（不含子女同住控制）。

## 文件
- rq3_recommended_core.csv：九格维度×状态 ATT（CSDID 含预趋势 p 与 BH-FDR；TWFE 基准）
- rq3_recommended_joint.csv：联合交互模型（Dmat/Dskill/Dsocial 同时进入）
- rq3_recommended_composition.csv：同 DVI 总分、不同构成的组内 ATT
- rq3_recommended_composition_contrasts.csv：构成间差异（含 FDR）
- rq3_recommended_sample_cells.csv / rq3_recommended_dvi3_counts.csv：样本量
- rq3_recommended_mechanism_desc.csv：轻量渠道描述表

## 关键结果速览
1. 高脆弱组（Dmat=1、Dskill=1、Dsocial=1）BHCI 与功能健康 ATT 均显著，且预趋势通过、FDR 后仍存（BHCI Dmat1/Dsocial1 q<0.01，Dskill1 q=0.019；功能 Dskill1/Dsocial1 q<0.04）。
2. 联合模型（CSDID）：BHCI 只有 Dsocial 梯度边际显著（0.022，p=0.035，q=0.104）；功能健康只有 Dskill 负梯度边际显著（-0.023，p=0.021，q=0.062）；Dmat 梯度被 Dsocial 吸收。
3. 同总分构成对比：唯一 FDR 后显著的对比是 DVI3=0.6667 组内 (0.5_0.5_1) vs (0.5_1_0.5)（BHCI 差 0.036，q=0.040），方向为“社会参与缺口 > 技能缺口”。
4. 两个非主格预趋势不过：Dmat=0（BHCI，p=0.003）、Dskill=0.5（功能，p=0.0009），论文中作为稳健性说明。
5. “控制总量后维度是否仍有解释力”的回归因 dvi3 与三维均值完全共线而不可用，已从脚本移除，不报告。
