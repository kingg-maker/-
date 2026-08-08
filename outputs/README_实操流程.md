# 宽带中国 x CHARLS 实操流水线（修订版方案落地）

## 一、数据审计结果（2026-08-04）

| 数据 | 结果 |
| --- | --- |
| CHARLS.csv | 176 列，96,629 行，5 波次（17,708 / 18,612 / 21,097 / 19,816 / 19,395） |
| 城市覆盖 | 126 个省份|城市单元 |
| 宏观 GDP | 296 个地级市，2011/2013/2015/2018/2020 |
| 宏观床位 | 297 个地级市，同上 |
| 城市匹配 | 116 个直接匹配，1 个手工映射（巢湖->合肥 340100），9 个自治州/盟/地区无宏观数据 |

### 关键缺口（必须先解决）
1. **宽带中国试点名单**：需要你提供 broadband_pilot.csv（province, city, batch, announce_year, announce_month）。脚本会因缺此文件停止。
2. **互联网相关变量**：现有 CHARLS.csv 中没有 ha011/ha012/cb001/互联网使用/未满足医疗需求。DVI 已用代理版（高龄、文盲、农村、低消费、无配偶）；若能从原始 CHARLS 提取 ha011/ha012/cb001，可升级为六维度原版 DVI。
3. **空间数据**：空间 SLX 需要城市坐标/邻接矩阵，脚本已留骨架。

## 二、执行顺序

1. 把宽带试点名单整理成 broadband_pilot.csv 放入 C:/Users/26301/Desktop/CHARLS_project/
2. 运行 01_data_prep.R：城市匹配、样本筛选、BHCI/CHE/DVI 构建、宏观合并
3. 运行 02_main_analysis.R：RQ1 C&S、RQ2a 堆叠交互+因果森林、RQ2b CI-DID
4. 运行 03_mechanism_robustness.R：RQ3 渠道、IPW、Lee Bounds、Bacon、空间骨架

## 三、输出
- analysis_df.rds / analysis_df.csv：最终分析面板
- RQ1-RQ3 结果表与事件研究图

## 四、审稿防守要点
- 处理编码已按时间线校准（第二批/第三批 2015 波次=0）
- 交互项只用基期变量（DVI、z 均基于 2011）
- 因果森林不做手动残差化（grf 内置 DML）
- 生存者偏差用 IPW + Lee Bounds 双层防线
