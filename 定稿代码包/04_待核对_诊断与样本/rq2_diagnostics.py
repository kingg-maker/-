# RQ2 diagnostics: control missingness, GDP balance, DVI=1 support.
# Read-only analysis over existing data; writes only new diagnostic CSVs.

import pandas as pd
from pathlib import Path

root = Path("C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq2_final_v1")
out_dir = root / "diagnostics"
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    root / "rq2_final_data.csv",
    dtype={"ID": str, "city_code": str},
    low_memory=False,
)
df["wave"] = pd.to_numeric(df["wave"], errors="coerce")
df["gvar"] = pd.to_numeric(df["gvar"], errors="coerce")

# ---------- control missingness ----------
main_ctrl = ["age_base", "gender_base", "rural_base", "edu_base"]
extra3 = ["marry_base", "hhcperc_base", "chronic_count_base"]
extra_gdp = ["marry_base", "hhcperc_base", "gdp_pc_log_base", "chronic_count_base"]
ext_full = [
    "marry_base", "hhcperc_base", "gdp_pc_log_base",
    "smoke_base", "drink_base", "ins_base", "children_base",
    "chronic_count_base",
]
all_ctrl = main_ctrl + ["marry_base", "hhcperc_base", "gdp_pc_log_base",
                        "smoke_base", "drink_base", "ins_base",
                        "children_base", "chronic_count_base"]

all_main = df[df["wave"].isin([1, 2, 3, 4])]
chronic_main = all_main[all_main["chronic_base"] == 1]

missing_rows = []
for sample_name, d in [("all_60plus_2011_2018", all_main),
                       ("chronic_2011_2018", chronic_main)]:
    for var in all_ctrl:
        miss = d[var].isna().sum()
        missing_rows.append({
            "sample": sample_name,
            "variable": var,
            "n_rows": len(d),
            "n_missing": int(miss),
            "missing_pct": round(100 * miss / len(d), 2),
        })
missing_tab = pd.DataFrame(missing_rows)
missing_tab.to_csv(out_dir / "rq2_diag_control_missing.csv", index=False)

sample_rows = []
for sample_name, d in [("all_60plus_2011_2018", all_main),
                       ("chronic_2011_2018", chronic_main)]:
    for set_name, vars_ in [
        ("main4", main_ctrl),
        ("main4_plus_3", main_ctrl + extra3),
        ("main4_plus_gdp", main_ctrl + extra_gdp),
        ("extended12", main_ctrl + ext_full),
    ]:
        sub = d.dropna(subset=vars_)
        sample_rows.append({
            "sample": sample_name,
            "control_set": set_name,
            "n_rows": len(sub),
            "n_ids": sub["ID"].nunique(),
            "row_retention_pct": round(100 * len(sub) / len(d), 2),
        })
sample_tab = pd.DataFrame(sample_rows)
sample_tab.to_csv(out_dir / "rq2_diag_sample_by_controls.csv", index=False)

# ---------- GDP balance at baseline ----------
base = df[df["wave"] == 1].copy()
base["ever_treat"] = (base["gvar"] > 0).astype(int)

def bal(v):
    x0 = base.loc[base["ever_treat"] == 0, v].dropna()
    x1 = base.loc[base["ever_treat"] == 1, v].dropna()
    m0, m1 = x0.mean(), x1.mean()
    s0, s1 = x0.std(ddof=1), x1.std(ddof=1)
    sp = ((x0.count() - 1) * s0**2 + (x1.count() - 1) * s1**2) / max(
        x0.count() + x1.count() - 2, 1
    )
    sp = sp**0.5
    return {
        "variable": v,
        "mean_never": round(m0, 4),
        "mean_ever": round(m1, 4),
        "std_diff": round((m1 - m0) / sp, 4) if sp else None,
        "n_never": int(x0.count()),
        "n_ever": int(x1.count()),
    }

bal_tab = pd.DataFrame([
    bal("gdp_pc_log_base"),
    bal("hhcperc_base"),
    bal("age_base"),
    bal("chronic_count_base"),
])
bal_tab.to_csv(out_dir / "rq2_diag_gdp_balance.csv", index=False)

# ---------- DVI=1 support ----------
df["dvi3"] = df[["dvi_mot", "dvi_mat", "dvi_skill"]].mean(axis=1, skipna=True).round(4)
df["dvi1_3dim"] = (df["dvi3"] == 1).astype(int)
df["dvi1_main"] = (df["dvi_main"].round(4) == 1).astype(int)

support_rows = []
for grp_name, grp_col in [("dvi3==1", "dvi1_3dim"), ("dvi_main==1", "dvi1_main")]:
    sub = df[df[grp_col] == 1]
    ever = sub[sub["gvar"] > 0]
    nev = sub[sub["gvar"] == 0]
    support_rows.append({
        "group": grp_name,
        "stat": "all_baseline_ids",
        "n": sub[sub["wave"] == 1]["ID"].nunique(),
    })
    support_rows.append({
        "group": grp_name,
        "stat": "all_person_wave_rows",
        "n": len(sub[sub["wave"].isin([1, 2, 3, 4])]),
    })
    support_rows.append({
        "group": grp_name,
        "stat": "all_cities",
        "n": sub["city_code"].nunique(),
    })
    for g in [3, 4]:
        gd = ever[ever["gvar"] == g]
        support_rows.append({
            "group": grp_name,
            "stat": f"gvar{g}_treated_baseline_ids",
            "n": gd[gd["wave"] == 1]["ID"].nunique(),
        })
        support_rows.append({
            "group": grp_name,
            "stat": f"gvar{g}_treated_rows",
            "n": len(gd[gd["wave"].isin([1, 2, 3, 4])]),
        })
        support_rows.append({
            "group": grp_name,
            "stat": f"gvar{g}_treated_cities",
            "n": gd["city_code"].nunique(),
        })
    support_rows.append({
        "group": grp_name,
        "stat": "never_baseline_ids",
        "n": nev[nev["wave"] == 1]["ID"].nunique(),
    })
    support_rows.append({
        "group": grp_name,
        "stat": "never_rows",
        "n": len(nev[nev["wave"].isin([1, 2, 3, 4])]),
    })
    support_rows.append({
        "group": grp_name,
        "stat": "never_cities",
        "n": nev["city_code"].nunique(),
    })
    for g in [3, 4]:
        gd = ever[ever["gvar"] == g]
        for w in [1, 2, 3, 4]:
            support_rows.append({
                "group": grp_name,
                "stat": f"gvar{g}_wave{w}_rows",
                "n": len(gd[gd["wave"] == w]),
            })

support_tab = pd.DataFrame(support_rows)
support_tab.to_csv(out_dir / "rq2_diag_dvi1_support.csv", index=False)

print("Missingness:")
print(missing_tab.to_string(index=False))
print("\nSample by control set:")
print(sample_tab.to_string(index=False))
print("\nGDP balance:")
print(bal_tab.to_string(index=False))
print("\nDVI=1 support:")
print(support_tab.to_string(index=False))
