# RQ1 + RQ2 unified rerun (2026-08-29)

Frozen spec: main sample = baseline chronic elderly (chronic_base==1),
window 2011-2018, never-treated controls, CS-DID (did::att_gt), city clustering.
Unified controls = age, gender, rural, edu, marry, hhcperc, city GDP, smoke,
drink, insurance, children, chronic count.

RQ1 chronic outputs:
- rq1_chronic_att_main.csv / rq1_chronic_att_ext.csv: 8 outcomes, CS-DID
- rq1_chronic_event_study_{func_cap,bhci_fixed,cog_cap}.csv/png
- rq1_chronic_pre_trend_*.csv
- rq1_chronic_baseline_table1.csv (2011, treated vs control, SMD)

RQ2 (3-dim DVI: motivation + material + skill) outputs:
- rq2_dvi3_distribution.csv
- rq2_dvi3_value_att_csdid.csv: ATT by each actual DVI3 value
- rq2_dvi3_value_table_word.csv: wide version for Word Table3
- rq2_dvi3_value_att_stacked.csv: complete by-value table (stacked DID),
  used for DVI3=1 where CS-DID cannot converge
- rq2_dvi3_tercile_att_csdid.csv + high-low contrast
- rq2_dvi3_linear_interactions.csv: TWFE / stacked robustness

Notes:
- DVI3 terciles are value-based: Low = 0-0.5, Mid = 0.6667/0.75,
  High = 0.8333/1 (quantile cuts split identical DVI values and made
  the High cell contain only DVI=1, which CS-DID cannot estimate).
- cov_used = ext: extended controls; main_fallback: only age/gender/rural/edu;
  failed: estimator did not converge for that cell.
- est_method/control_group: DR + never-treated preferred; ipw/reg or
  not-yet-treated controls are used only when the preferred estimator
  fails (extreme DVI cells with collinear or tiny control groups).
- n in att tables is person-wave rows after complete-case controls.
- baseline_n is unique baseline persons with nonmissing outcome in that cell.
- edu_base 1-4 coded as none / primary / junior+ for Table 1 (verify with
  CHARLS codebook before final submission).
