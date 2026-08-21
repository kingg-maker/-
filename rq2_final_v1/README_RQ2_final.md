# RQ2 final run (2026-08-21)

Main sample: baseline chronic elderly (chronic_base==1), 2011-2018.
Main DVI: material + skill equal-weight; four-dim versions are robustness.
Main estimator: stacked DID (cohort x unit + cohort x wave), city clustering.

Outputs:
- rq2_final_data.csv: merged analysis data with DVI versions
- rq2a_stacked_models.csv: linear interaction models
- rq2a_dvi_tercile_att.csv: DVI tercile ATTs
- rq2a_high_low_test.csv: High-Low Wald test
- rq2a_causal_forest.rds / importance csv: exploratory causal forest
- rq2b_ci_did.csv / erreygers csv: income-dimension results

Data supplements used: none beyond analysis_df_v2.csv and revised_data.csv.
