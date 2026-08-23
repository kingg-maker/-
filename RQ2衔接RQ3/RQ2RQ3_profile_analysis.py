# -*- coding: utf-8 -*-
"""
RQ2 -> RQ3 bridging descriptive tables.

Table 1: vulnerability shares by the nine values of the three-dimension DVI.
Table 2: baseline (2011) profile of the high-vulnerability tail
         (DVI3 in {0.8333, 1.0}) versus the rest, chronic-base main sample.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Users\26301\Documents\Codex\2026-08-23\c-users-26301-documents-codex-2026\outputs"
DF_PATH = r"C:\Users\26301\Documents\Codex\2026-08-15\new-chat\outputs\rq2_final_v1\rq2_final_data.csv"
FAM_PATH = r"C:\Users\26301\Documents\Codex\2026-08-07\w\work\raw\新增数据\2011\family_information.dta"
HH_PATH = r"C:\Users\26301\Documents\Codex\2026-08-07\w\work\raw\新增数据\2011\hhmember.dta"
CHARLS_PATH = r"C:\Users\26301\Documents\Codex\2026-08-07\w\work\raw\CHARLS.csv"


def norm_id_2011(value: object) -> object:
    """Insert a zero before the final digit (matches process_raw.py)."""
    if pd.isna(value):
        return value
    s = str(value).strip()
    if len(s) >= 2:
        return s[:-1] + "0" + s[-1]
    return s


def chi2_pvalue(obs: np.ndarray) -> float:
    """Chi-square survival probability without scipy."""
    exp = np.outer(obs.sum(axis=1), obs.sum(axis=0)) / obs.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.nansum(np.where(exp > 0, (obs - exp) ** 2 / exp, 0.0))
    df = (obs.shape[0] - 1) * (obs.shape[1] - 1)
    if df <= 0 or stat <= 0:
        return 1.0
    if stat > 300:
        return 0.0
    a = df / 2.0
    x = stat / 2.0
    # Regularized lower incomplete gamma P(a, x) via a log-space series.
    log_term = 0.0
    log_total = 0.0
    k = 0
    while k < 500:
        log_term += math.log(x) - math.log(a + k)
        m = max(log_total, log_term)
        log_total = m + math.log(math.exp(log_total - m) + math.exp(log_term - m))
        if log_term < log_total - 50:
            break
        k += 1
    p_low = math.exp(a * math.log(x) - x - math.lgamma(a) + log_total)
    return min(max(1.0 - p_low, 0.0), 1.0)


def prop_pvalue(a1: int, n1: int, a0: int, n0: int) -> float:
    """2x2 chi-square with Yates continuity correction."""
    tab = np.array([[a1, n1 - a1], [a0, n0 - a0]], dtype=float)
    exp = np.outer(tab.sum(axis=1), tab.sum(axis=0)) / tab.sum()
    stat = np.sum((np.abs(tab - exp) - 0.5) ** 2 / exp)
    return _chi2_1df(stat)


def _chi2_1df(stat: float) -> float:
    if stat <= 0:
        return 1.0
    return math.erfc(math.sqrt(stat / 2.0))


def welch_pvalue(x1: np.ndarray, x0: np.ndarray) -> float:
    x1 = x1[~np.isnan(x1)]
    x0 = x0[~np.isnan(x0)]
    if len(x1) < 2 or len(x0) < 2:
        return np.nan
    m1, m0 = x1.mean(), x0.mean()
    v1, v0 = x1.var(ddof=1), x0.var(ddof=1)
    se = math.sqrt(v1 / len(x1) + v0 / len(x0))
    if se == 0:
        return 1.0
    t = (m1 - m0) / se
    return 2 * math.erfc(abs(t) / math.sqrt(2.0))


def q_stats(s: pd.Series) -> dict:
    v = s.dropna().astype(float)
    if len(v) == 0:
        return {"n": 0, "mean": np.nan, "sd": np.nan, "p25": np.nan,
                "median": np.nan, "p75": np.nan}
    return {
        "n": len(v),
        "mean": v.mean(),
        "sd": v.std(ddof=1),
        "p25": v.quantile(0.25),
        "median": v.median(),
        "p75": v.quantile(0.75),
    }


def pct_str(n: int, d: int) -> str:
    if d == 0:
        return "NA"
    return f"{100 * n / d:.1f}%"


def num(v: float, digits: int = 2) -> str:
    if pd.isna(v):
        return "NA"
    return f"{v:.{digits}f}"


def main() -> None:
    df = pd.read_csv(DF_PATH, low_memory=False, dtype={"ID": "string"})
    df["dvi3"] = df[["dvi_mot", "dvi_mat", "dvi_skill"]].mean(axis=1, skipna=True)

    # Chronic-base main sample, baseline wave only.
    base = df[(df["wave"] == 1) & (df["chronic_base"] == 1)].copy()
    base["dvi3_round"] = base["dvi3"].round(4)

    # Family module (2011 raw IDs need normalization).
    fam_cols = (
        ["ID"]
        + [f"cb053_{i}_" for i in range(1, 15)]
        + [f"cd004_{i}_" for i in range(1, 15)]
    )
    fam = pd.read_stata(FAM_PATH, convert_categoricals=False, columns=fam_cols)
    fam["ID"] = fam["ID"].map(norm_id_2011).astype("string")
    fam = fam.drop_duplicates(subset=["ID"])

    ch = pd.read_csv(
        CHARLS_PATH,
        low_memory=False,
        usecols=["ID", "wave", "family_size", "hchild", "fcamt", "adlab_c", "iadl"],
        dtype={"ID": "string"},
    )
    ch = ch[ch["wave"] == 1].drop_duplicates(subset=["ID"])

    hh = pd.read_stata(HH_PATH, convert_categoricals=False, columns=["ID", "a006"])
    hh["ID"] = hh["ID"].map(norm_id_2011).astype("string")
    hh["live_with_child"] = hh["a006"].isin([7, 8])
    hh_live = hh.groupby("ID")["live_with_child"].max().reset_index()

    base = (
        base.merge(ch, on="ID", how="left")
        .merge(fam, on="ID", how="left")
        .merge(hh_live, on="ID", how="left")
    )
    base["live_with_child"] = base["live_with_child"].fillna(0)

    # Derived family / support variables.
    cb = base[[f"cb053_{i}_" for i in range(1, 15)]]
    cd = base[[f"cd004_{i}_" for i in range(1, 15)]]
    noncores_cond = cb.isin([3, 4, 5, 6, 7, 8, 9, 10]).to_numpy()
    contact = cd.where(noncores_cond)
    base["any_child_near_cb"] = cb.isin([1, 2]).any(axis=1)
    base["college_child_near"] = (
        (base["has_college_child_2011"] == 1) & (base["live_with_child"] == 1)
    )
    best = contact.apply(lambda r: r[r.between(1, 9)].min(), axis=1)
    base["contact_almost_daily"] = best.eq(1)
    base["contact_weekly_plus"] = best.between(2, 3)
    base["contact_monthly_plus"] = best.between(4, 5)
    base["contact_less"] = best.between(6, 9)
    base["has_noncores_child"] = cd.notna().any(axis=1)
    base["has_any_child"] = base["hchild"].fillna(0) > 0

    base["female"] = (base["gender_base"] == 0).astype(float)
    base["rural_bin"] = (base["rural_base"] == 1).astype(float)
    base["married"] = (base["marry_base"] == 1).astype(float)
    base["live_alone"] = (base["family_size"] == 1).astype(float)
    base["fcamt_pos"] = (base["fcamt"] > 0).astype(float)
    base["adl_limited"] = (base["adlab_c"] > 0).astype(float)
    base["iadl_limited"] = (base["iadl"] > 0).astype(float)
    base["cesd_total"] = (1 - base["psych_cap"]) * 30
    base["log_hhcperc"] = np.log1p(base["hhcperc_base"])
    base["age60_69"] = base["age_base"].between(60, 69).astype(float)
    base["age70_79"] = base["age_base"].between(70, 79).astype(float)
    base["age80p"] = (base["age_base"] >= 80).astype(float)
    base["edu_none"] = (base["edu_base"] == 1).astype(float)
    base["edu_primary"] = (base["edu_base"] == 2).astype(float)
    base["edu_junior_plus"] = (base["edu_base"] >= 3).astype(float)
    base["srh_bad"] = (base["srh_cap_fixed"] <= 0.25).astype(float)
    base["srh_fair"] = (base["srh_cap_fixed"] == 0.5).astype(float)
    base["srh_good"] = (base["srh_cap_fixed"] >= 0.75).astype(float)

    # ---------------- Table 1: nine DVI values ----------------
    order = [0.0, 0.1667, 0.25, 0.3333, 0.5, 0.6667, 0.75, 0.8333, 1.0]
    rows1 = []
    for val in order:
        g = base[base["dvi3_round"] == val]
        if len(g) == 0:
            continue
        rows1.append({
            "DVI值": f"{val:.4f}",
            "N": len(g),
            "无电脑(%)": pct_str(int(g["no_computer"].sum()), int(g["no_computer"].notna().sum())),
            "无手机(%)": pct_str(int(g["no_mobile"].sum()), int(g["no_mobile"].notna().sum())),
            "低教育(%)": pct_str(int(g["low_edu"].sum()), int(g["low_edu"].notna().sum())),
            "无高学历子女(%)": pct_str(int(g["no_college"].sum()), int(g["no_college"].notna().sum())),
            "无社会信息活动(%)": pct_str(int((g["social_info"] == 0).sum()), int(g["social_info"].notna().sum())),
            "无学习活动(%)": pct_str(int((g["learning"] == 0).sum()), int(g["learning"].notna().sum())),
            "动机维度全脆弱(%)": pct_str(int((g["dvi_mot"] == 1).sum()), int(g["dvi_mot"].notna().sum())),
            "物质维度均值": num(g["dvi_mat"].mean()),
            "技能维度均值": num(g["dvi_skill"].mean()),
            "动机维度均值": num(g["dvi_mot"].mean()),
        })
    tab1 = pd.DataFrame(rows1)

    # ---------------- Table 2: high tail vs others ----------------
    high = base["dvi3_round"].isin([0.8333, 1.0])
    base["group"] = np.where(high, "High", "Non-high")
    hi = base[high]
    nh = base[~high]

    cont_vars = [
        ("年龄", "age_base"),
        ("log家庭人均消费", "log_hhcperc"),
        ("认知健康(0-31)", "cog_cap"),
        ("CES-D总分(0-30)", "cesd_total"),
        ("自评健康(0-1)", "srh_cap_fixed"),
    ]
    rows2 = []
    for label, col in cont_vars:
        a = q_stats(hi[col])
        b = q_stats(nh[col])
        p = welch_pvalue(hi[col].to_numpy(dtype=float), nh[col].to_numpy(dtype=float))
        rows2.append({
            "类别": "连续变量", "变量": label, "指标": "Mean",
            "High": num(a["mean"]), "Non-high": num(b["mean"]), "P值": num(p, 3),
        })
        rows2.append({
            "类别": "连续变量", "变量": label, "指标": "SD",
            "High": num(a["sd"]), "Non-high": num(b["sd"]), "P值": "",
        })
        rows2.append({
            "类别": "连续变量", "变量": label, "指标": "Median",
            "High": num(a["median"]), "Non-high": num(b["median"]), "P值": "",
        })
        rows2.append({
            "类别": "连续变量", "变量": label, "指标": "P25/P75",
            "High": f"{num(a['p25'])}/{num(a['p75'])}",
            "Non-high": f"{num(b['p25'])}/{num(b['p75'])}",
            "P值": "",
        })

    bin_vars = [
        ("女性(%)", "female"),
        ("农村(%)", "rural_bin"),
        ("已婚/有配偶(%)", "married"),
        ("独居(%)", "live_alone"),
        ("子女或子女配偶同住(%)", "live_with_child"),
        ("子女同住或同院相邻(cb053)(%)", "any_child_near_cb"),
        ("高学历子女在身边(%)", "college_child_near"),
        ("收到子女经济支持(%)", "fcamt_pos"),
        ("ADL至少一项受限(%)", "adl_limited"),
        ("IADL至少一项受限(%)", "iadl_limited"),
        ("自评健康差(%)", "srh_bad"),
        ("自评健康一般(%)", "srh_fair"),
        ("自评健康好(%)", "srh_good"),
        ("年龄60-69(%)", "age60_69"),
        ("年龄70-79(%)", "age70_79"),
        ("年龄80+(%)", "age80p"),
        ("未受教育(%)", "edu_none"),
        ("小学(%)", "edu_primary"),
        ("初中及以上(%)", "edu_junior_plus"),
    ]
    for label, col in bin_vars:
        n1 = int(hi[col].sum())
        n0 = int(nh[col].sum())
        p = prop_pvalue(n1, len(hi), n0, len(nh))
        rows2.append({
            "类别": "二元变量", "变量": label, "指标": "比例",
            "High": pct_str(n1, len(hi)), "Non-high": pct_str(n0, len(nh)),
            "P值": num(p, 3),
        })

    # Multi-category rows.
    mult = [
        ("年龄组(60-69/70-79/80+)", ["age60_69", "age70_79", "age80p"]),
        ("教育(未受/小学/初中及以上)", ["edu_none", "edu_primary", "edu_junior_plus"]),
        ("自评健康(差/一般/好)", ["srh_bad", "srh_fair", "srh_good"]),
    ]
    for label, cols in mult:
        tab = np.array([[hi[c].sum() for c in cols], [nh[c].sum() for c in cols]])
        p = chi2_pvalue(tab)
        rows2.append({
            "类别": "分组变量", "变量": label, "指标": "组间卡方P值",
            "High": "", "Non-high": "", "P值": num(p, 3),
        })

    # Contact-frequency categories.
    freq_cols = ["contact_almost_daily", "contact_weekly_plus",
                 "contact_monthly_plus", "contact_less"]
    tab_c = np.array([[hi[c].sum() for c in freq_cols],
                      [nh[c].sum() for c in freq_cols]])
    p_c = chi2_pvalue(tab_c)
    rows2.append({
        "类别": "分组变量", "变量": "子女联系频率(非共居子女中最频繁)",
        "指标": "组间卡方P值", "High": "", "Non-high": "", "P值": num(p_c, 3),
    })
    for label, col in zip(
        ["几乎每天(%)", "每周1-3次(%)", "每月1-2次(%)", "更少或几乎不(%)"],
        freq_cols,
    ):
        rows2.append({
            "类别": "子女联系频率", "变量": label, "指标": "比例",
            "High": pct_str(int(hi[col].sum()), len(hi)),
            "Non-high": pct_str(int(nh[col].sum()), len(nh)),
            "P值": "",
        })
    rows2.append({
        "类别": "子女联系频率", "变量": "全部子女同住/相邻(无非同住子女)(%)",
        "指标": "比例",
        "High": pct_str(
            int((hi["has_any_child"] & ~hi["has_noncores_child"]).sum()),
            int(hi["has_any_child"].sum()),
        ),
        "Non-high": pct_str(
            int((nh["has_any_child"] & ~nh["has_noncores_child"]).sum()),
            int(nh["has_any_child"].sum()),
        ),
        "P值": "",
    })

    # Economic support amount (yuan), conditional on any support.
    rows2.append({
        "类别": "连续变量", "变量": "子女经济支持金额(全部样本, 元)", "指标": "Mean",
        "High": num(hi["fcamt"].mean()), "Non-high": num(nh["fcamt"].mean()),
        "P值": num(welch_pvalue(hi["fcamt"].to_numpy(float), nh["fcamt"].to_numpy(float)), 3),
    })
    rows2.append({
        "类别": "连续变量", "变量": "子女经济支持金额(全部样本, 元)", "指标": "Median",
        "High": num(hi["fcamt"].median()), "Non-high": num(nh["fcamt"].median()), "P值": "",
    })
    rows2.append({
        "类别": "连续变量", "变量": "子女经济支持金额(有支持者, 元)", "指标": "Mean",
        "High": num(hi.loc[hi["fcamt_pos"] == 1, "fcamt"].mean()),
        "Non-high": num(nh.loc[nh["fcamt_pos"] == 1, "fcamt"].mean()), "P值": "",
    })

    tab2 = pd.DataFrame(rows2)

    tab1.to_csv(f"{OUT_DIR}/RQ2RQ3_table1_dvi9.csv", index=False, encoding="utf-8-sig")
    tab2.to_csv(f"{OUT_DIR}/RQ2RQ3_profile_high_vs_nonhigh.csv", index=False, encoding="utf-8-sig")

    n_hi = len(hi)
    n_nh = len(nh)
    md = []
    md.append("# RQ2→RQ3 衔接：DVI 取值描述与高脆弱尾部画像")
    md.append("")
    md.append(f"样本：2011 年基期 60 岁及以上、城市可匹配、`chronic_base=1` 的慢病人群（RQ2 主样本），共 **{len(base)}** 人；")
    md.append(f"高脆弱尾部（DVI3=0.8333、1.0）**{n_hi}** 人，其余 **{n_nh}** 人。")
    md.append("")
    md.append("## 表 1：按三维 DVI 九个取值分组的维度弱势占比")
    md.append("")
    md.append("| " + " | ".join(tab1.columns) + " |")
    md.append("|" + "|".join(["---"] * len(tab1.columns)) + "|")
    for _, r in tab1.iterrows():
        md.append("| " + " | ".join(str(r[c]) for c in tab1.columns) + " |")
    md.append("")
    md.append("## 表 2：高脆弱尾部 vs 其他人群（2011 基期画像）")
    md.append("")
    md.append("| 类别 | 变量 | 指标 | High | Non-high | P值 |")
    md.append("|---|---|---|---|---|---|")
    for _, r in tab2.iterrows():
        md.append(f"| {r['类别']} | {r['变量']} | {r['指标']} | {r['High']} | {r['Non-high']} | {r['P值']} |")
    md.append("")
    md.append("## 变量定义与口径")
    md.append("")
    md.append("- DVI3 = 动机 + 物质 + 技能 三维等权（与 RQ2 三维主口径一致）；数值为 0、0.1667、0.25、0.3333、0.5、0.6667、0.75、0.8333、1.0 共 9 档。")
    md.append("- 物质：无电脑、无手机；技能：低教育（文盲/未脱盲）、无大专及以上学历子女；动机：无社会信息活动（DA056 1/2/4/5）、无学习活动（DA056 6/8/9）代理。")
    md.append("- 与子女同住：2011 年家庭名册中任意成员与受访者关系为子女或子女配偶（`hhmember.a006`=7/8）；另附 `cb053`=1/2（本户或同院相邻）口径作为对照。")
    md.append("- 子女联系频率：对非同住子女的 `cd004`（`cb053`>=3 的非共居子女），取最频繁者（代码 1=几乎每天，2-3=每周，4-5=每月，6-9=更少/几乎不）。")
    md.append("- 子女经济支持：`fcamt`（子女对父母经济支持，元）。")
    md.append("- ADL 受限 = `adlab_c>0`；IADL 受限 = `iadl>0`；CES-D 总分 = `(1-psych_cap)×30`；自评健康：差≤0.25、一般=0.5、好≥0.75。")
    md.append("- log 家庭人均消费采用 `log(1+hhcperc_base)`，避免基期消费为 0 时无定义。")
    md.append("- 连续变量组间差异为 Welch t 近似检验；比例差异为 2×2 卡方（Yates 校正）；多分类为 r×2 卡方。画像为描述性统计，P 值仅作参考。")
    md.append("")
    with open(f"{OUT_DIR}/RQ2RQ3_画像表.md", "w", encoding="utf-8") as f:
        f.write("\n".join(md))

    print("N high:", n_hi, "N non-high:", n_nh)
    print(tab1.to_string(index=False))
    print(tab2.to_string(index=False))
    print("Saved:", f"{OUT_DIR}/RQ2RQ3_画像表.md")


if __name__ == "__main__":
    main()
