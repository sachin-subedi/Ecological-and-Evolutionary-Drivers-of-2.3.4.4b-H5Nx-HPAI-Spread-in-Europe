import time
from io import StringIO
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from dateutil import parser
import baltic as bt

ANCHOR_DEC_YEAR = 2025.2520547945205

GEOCLUSTERS = [
    "HC_Cluster_1_Alpine",
    "HC_Cluster_1_Atlantic",
    "HC_Cluster_1_Continental",
    "HC_Cluster_2_Alpine",
    "HC_Cluster_2_Continental",
    "HC_Cluster_2_Mediterranean",
    "HC_Cluster_2_Pannonian",
    "HC_Cluster_3_Alpine",
    "HC_Cluster_3_Boreal",
    "HC_Cluster_4_Mediterranean",
]

GEOCLUSTER_COLORS = {
    "HC_Cluster_1_Alpine"        : "#CC79A7",
    "HC_Cluster_1_Atlantic"      : "#CD5C5C",
    "HC_Cluster_1_Continental"   : "#BCBD22",
    "HC_Cluster_2_Alpine"        : "#26A69A",
    "HC_Cluster_2_Continental"   : "#0072B2",
    "HC_Cluster_2_Mediterranean" : "#708090",
    "HC_Cluster_2_Pannonian"     : "#6B4C1B",
    "HC_Cluster_3_Alpine"        : "#7570B3",
    "HC_Cluster_3_Boreal"        : "#006D2C",
    "HC_Cluster_4_Mediterranean" : "#56B4E9",
}

GEOCLUSTER_LABELS = {
    "HC_Cluster_1_Alpine"        : "Central Alpine",
    "HC_Cluster_1_Atlantic"      : "Atlantic",
    "HC_Cluster_1_Continental"   : "Western Continental",
    "HC_Cluster_2_Alpine"        : "Eastern Alpine",
    "HC_Cluster_2_Continental"   : "Eastern Continental",
    "HC_Cluster_2_Mediterranean" : "Southeast Mediterranean",
    "HC_Cluster_2_Pannonian"     : "Pannonian",
    "HC_Cluster_3_Alpine"        : "Scandinavian Highlands",
    "HC_Cluster_3_Boreal"        : "Boreal Baltic",
    "HC_Cluster_4_Mediterranean" : "Iberian",
}
GEOCLUSTER_LABELS_ORDERED = [GEOCLUSTER_LABELS[g] for g in GEOCLUSTERS]

START      = pd.to_datetime("2016-01-01")
END        = pd.to_datetime("2025-04-03")
ALL_MONTHS = pd.date_range(START, END, freq="MS")
YEAR_TICKS = pd.date_range(START, END, freq="YS")
YEARS      = list(range(START.year, END.year + 1))

DATASETS = [
    {
        "label"  : "Equal_rep1",
        "trees"  : (
            "GeoCluster_history_equal_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_equal_rep1",
    },
    {
        "label"  : "Equal_rep2",
        "trees"  : (
            "GeoCluster_history_equal_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_equal_rep2",
    },
    {
        "label"  : "Equal_rep3",
        "trees"  : (
            "GeoCluster_history_equal_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_equal_rep3",
    },
    {
        "label"  : "Proportional_rep1",
        "trees"  : (
            "GeoCluster_history_proportional_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_proportional_rep1",
    },
    {
        "label"  : "Proportional_rep2",
        "trees"  : (
            "GeoCluster_history_proportional_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_proportional_rep2",
    },
    {
        "label"  : "Proportional_rep3",
        "trees"  : (
            "GeoCluster_history_proportional_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_proportional_rep3",
    },
    {
        "label"  : "Stratified_rep1",
        "trees"  : (
            "GeoCluster_history_stratified_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_stratified_rep1",
    },
    {
        "label"  : "Stratified_rep2",
        "trees"  : (
            "GeoCluster_history_stratified_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_stratified_rep2",
    },
    {
        "label"  : "Stratified_rep3",
        "trees"  : (
            "GeoCluster_history_stratified_combined.trees"
        ),
        "out_prefix": "GeoCluster_Reward_stratified_rep3",
    },
]

def decyear_to_month(dec):
    yr   = int(dec)
    frac = dec - yr
    dt   = datetime(yr, 1, 1) + pd.Timedelta(days=frac * 365.25)
    return dt.strftime("%Y-%m")

def hpd_interval(values, cred_mass=0.95):
    s = np.sort(np.asarray(values, dtype=float))
    s = s[~np.isnan(s)]
    n = len(s)
    if n == 0:
        return (np.nan, np.nan)
    if n == 1:
        return (s[0], s[0])
    interval_idx = int(np.floor(cred_mass * n))
    if interval_idx == 0:
        return (s[0], s[-1])
    n_windows = n - interval_idx
    widths    = s[interval_idx:] - s[:n_windows]
    min_i     = int(np.argmin(widths))
    return (s[min_i], s[min_i + interval_idx])

def eti_interval(values, cred_mass=0.95):
    s = np.asarray(values, dtype=float)
    s = s[~np.isnan(s)]
    if s.size == 0:
        return (np.nan, np.nan)
    lo_q = (1.0 - cred_mass) / 2.0 * 100.0
    hi_q = (1.0 + cred_mass) / 2.0 * 100.0
    return (np.percentile(s, lo_q), np.percentile(s, hi_q))

def parse_trees(trees_file):
    if not Path(trees_file).exists():
        raise FileNotFoundError(f"Cannot find: {trees_file}")

    start                 = time.time()
    tree_counter          = 0
    posterior_prop        = {}
    all_months_set        = set()
    per_tree_rewards      = []
    per_tree_year_rewards = []

    with open(trees_file) as fh:
        for ln in fh:
            if not ln.lower().startswith("tree state_"):
                continue
            tree_counter += 1
            try:
                newick = ln.split("= [&R] ")[1]
                tree   = bt.loadNewick(StringIO(newick), absoluteTime=False)
                tree.setAbsoluteTime(ANCHOR_DEC_YEAR)
            except Exception:
                continue

            month_cluster_raw  = {}
            month_totals_raw   = {}
            tree_reward_totals = {cl: 0.0 for cl in GEOCLUSTERS}
            tree_year_reward   = {} 

            for node in tree.Objects:
                if not node.traits:
                    continue
                month = decyear_to_month(node.absoluteTime)
                yr    = int(node.absoluteTime)
                all_months_set.add(month)

                node_rewards = {}
                for cl in GEOCLUSTERS:
                    val = node.traits.get(f"{cl}_reward")
                    if val is None:
                        continue
                    val = float(val)
                    node_rewards[cl] = val
                    tree_reward_totals[cl] += val
                    tree_year_reward[(yr, cl)] = tree_year_reward.get((yr, cl), 0.0) + val
                    pair_key = (month, cl)
                    month_cluster_raw[pair_key] = (
                        month_cluster_raw.get(pair_key, 0.0) + val
                    )
                node_total = float(node.traits.get("GeoCluster.count", 0.0))
                if node_total == 0:
                    node_total = sum(node_rewards.values())
                month_totals_raw[month] = month_totals_raw.get(month, 0.0) + node_total

            for (month, cl), raw_val in month_cluster_raw.items():
                total = month_totals_raw.get(month, 0.0)
                if total == 0:
                    continue
                prop = raw_val / total
                posterior_prop[(month, cl)] = posterior_prop.get((month, cl), 0.0) + prop

            per_tree_rewards.append(tree_reward_totals)
            per_tree_year_rewards.append(tree_year_reward)

    print(f"  ✓ {tree_counter:,} trees in {time.time()-start:.1f}s")
    return tree_counter, posterior_prop, per_tree_rewards, per_tree_year_rewards

def make_reward_proportions(tree_counter, posterior_prop):
    rows = []
    for (mm, gc), summed in posterior_prop.items():
        rows.append({
            "YearMonth" : pd.to_datetime(mm + "-01"),
            "GeoCluster": gc,
            "prop"      : summed / tree_counter,
        })

    if not rows:
        return pd.DataFrame(index=ALL_MONTHS, columns=GEOCLUSTERS).fillna(0.0)

    df  = pd.DataFrame(rows)
    mat = (
        df.pivot_table(
            index="YearMonth", columns="GeoCluster",
            values="prop", fill_value=0.0,
        )
        .reindex(ALL_MONTHS)
        .fillna(0.0)
    )
    for gc in GEOCLUSTERS:
        if gc not in mat.columns:
            mat[gc] = 0.0
    mat = mat[GEOCLUSTERS]

    row_sums = mat.sum(axis=1)
    mat = mat.div(row_sums.replace(0, pd.NA), axis=0).fillna(0.0)

    zero_rows = mat.sum(axis=1) == 0
    mat.loc[zero_rows] = pd.NA
    mat = mat.ffill().fillna(0.0)

    return mat


def make_reward_table(per_tree_rewards, label):
    if not per_tree_rewards:
        return pd.DataFrame()

    rdf            = pd.DataFrame(per_tree_rewards).reindex(columns=GEOCLUSTERS).fillna(0.0)
    total_per_tree = rdf.sum(axis=1)
    safe_total     = total_per_tree.replace(0, np.nan)

    rows = []
    for gc in GEOCLUSTERS:
        vals = rdf[gc].to_numpy(dtype=float)
        prop = (rdf[gc] / safe_total).to_numpy(dtype=float)
        r_hlo, r_hhi = hpd_interval(vals)
        r_clo, r_chi = eti_interval(vals)
        p_hlo, p_hhi = hpd_interval(prop)
        p_clo, p_chi = eti_interval(prop)
        rows.append({
            "Dataset"             : label,
            "GeoCluster"          : gc,
            "reward_mean_yr"      : np.mean(vals),
            "reward_median_yr"    : np.median(vals),
            "reward_HPD95_lo_yr"  : r_hlo,
            "reward_HPD95_hi_yr"  : r_hhi,
            "reward_CrI95_lo_yr"  : r_clo,
            "reward_CrI95_hi_yr"  : r_chi,
            "prop_mean"           : np.nanmean(prop),
            "prop_median"         : np.nanmedian(prop),
            "prop_HPD95_lo"       : p_hlo,
            "prop_HPD95_hi"       : p_hhi,
            "prop_CrI95_lo"       : p_clo,
            "prop_CrI95_hi"       : p_chi,
        })
    out = pd.DataFrame(rows)
    out["GeoCluster"] = out["GeoCluster"].map(GEOCLUSTER_LABELS)
    return out

def make_yearly_reward_table(per_tree_year_rewards, label):
    if not per_tree_year_rewards:
        return pd.DataFrame()

    dwell_samples = {(y, cl): [] for y in YEARS for cl in GEOCLUSTERS}
    prop_samples  = {(y, cl): [] for y in YEARS for cl in GEOCLUSTERS}

    for tree_yr in per_tree_year_rewards:
        year_tot = {y: 0.0 for y in YEARS}
        for (y, cl), v in tree_yr.items():
            if y in year_tot:
                year_tot[y] += v
        for y in YEARS:
            tot = year_tot[y]
            for cl in GEOCLUSTERS:
                v = tree_yr.get((y, cl), 0.0)
                dwell_samples[(y, cl)].append(v)
                prop_samples[(y, cl)].append(v / tot if tot > 0 else np.nan)

    rows = []
    for y in YEARS:
        for cl in GEOCLUSTERS:
            dv = np.asarray(dwell_samples[(y, cl)], dtype=float)
            pv = np.asarray(prop_samples[(y, cl)], dtype=float)
            r_hlo, r_hhi = hpd_interval(dv)
            r_clo, r_chi = eti_interval(dv)
            p_hlo, p_hhi = hpd_interval(pv)
            p_clo, p_chi = eti_interval(pv)
            rows.append({
                "Dataset"            : label,
                "Year"               : y,
                "GeoCluster"         : cl,
                "reward_mean_yr"     : np.nanmean(dv) if dv.size else np.nan,
                "reward_median_yr"   : np.nanmedian(dv) if dv.size else np.nan,
                "reward_HPD95_lo_yr" : r_hlo,
                "reward_HPD95_hi_yr" : r_hhi,
                "reward_CrI95_lo_yr" : r_clo,
                "reward_CrI95_hi_yr" : r_chi,
                "prop_pct_mean"      : np.nanmean(pv) * 100,
                "prop_pct_median"    : np.nanmedian(pv) * 100,
                "prop_pct_HPD95_lo"  : p_hlo * 100,
                "prop_pct_HPD95_hi"  : p_hhi * 100,
                "prop_pct_CrI95_lo"  : p_clo * 100,
                "prop_pct_CrI95_hi"  : p_chi * 100,
            })
    out = pd.DataFrame(rows)
    out["GeoCluster"] = out["GeoCluster"].map(GEOCLUSTER_LABELS)
    return out

def save_yearly_wide(yearly_table, out_prefix):
    if yearly_table.empty:
        return None
    prop_wide = (
        yearly_table.pivot(index="GeoCluster", columns="Year", values="prop_pct_mean")
        .reindex(index=GEOCLUSTER_LABELS_ORDERED, columns=YEARS)
    )
    dwell_wide = (
        yearly_table.pivot(index="GeoCluster", columns="Year", values="reward_mean_yr")
        .reindex(index=GEOCLUSTER_LABELS_ORDERED, columns=YEARS)
    )
    prop_wide.to_csv(f"{out_prefix}_reward_by_year_prop_pct.csv")
    dwell_wide.to_csv(f"{out_prefix}_reward_by_year_dwell.csv")
    return prop_wide

def plot_reward(reward_prop_df, save_prefix):
    colors = [GEOCLUSTER_COLORS[g] for g in GEOCLUSTERS]

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("white")

    ax.stackplot(ALL_MONTHS, [reward_prop_df[g].values for g in GEOCLUSTERS],
                 colors=colors, linewidth=0)

    ax.set_facecolor("white")
    ax.set_ylim(0, 1)
    ax.set_xlim(START, END)
    ax.set_ylabel("Reward proportion", fontsize=18)
    ax.set_xlabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    ax.tick_params(axis="y", labelsize=18)

    ax.set_xticks(YEAR_TICKS)
    ax.set_xticklabels([d.year for d in YEAR_TICKS], fontsize=18, rotation=45, ha="right")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")  # force override pandas default

    plt.tight_layout()
    fig.savefig(f"{save_prefix}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{save_prefix}.pdf", bbox_inches="tight")
    plt.show()
    plt.close(fig)
    print(f"  ✓ saved → {save_prefix}.png / .pdf")

def plot_yearly_heatmap(prop_wide, save_prefix):
    if prop_wide is None or prop_wide.empty:
        return
    mat = prop_wide.reindex(index=GEOCLUSTER_LABELS_ORDERED, columns=YEARS).to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.patch.set_facecolor("white")

    im = ax.imshow(mat, aspect="auto", cmap="magma_r", vmin=0,
                   vmax=np.nanmax(mat) if np.isfinite(np.nanmax(mat)) else 1)

    ax.set_xticks(range(len(YEARS)))
    ax.set_xticklabels(YEARS, fontsize=13)
    ax.set_yticks(range(len(GEOCLUSTER_LABELS_ORDERED)))
    ax.set_yticklabels(GEOCLUSTER_LABELS_ORDERED, fontsize=12)

    thr = np.nanmax(mat) * 0.55 if np.isfinite(np.nanmax(mat)) else 50
    for i in range(len(GEOCLUSTERS)):
        for j in range(len(YEARS)):
            v = mat[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    fontsize=10, color="white" if v > thr else "black")
    for j in range(len(YEARS)):
        col = mat[:, j]
        if np.all(np.isnan(col)):
            continue
        i_max = int(np.nanargmax(col))
        ax.add_patch(plt.Rectangle((j - 0.5, i_max - 0.5), 1, 1,
                     fill=False, edgecolor="#111111", lw=2.0))

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Within-year reward (%)", fontsize=13)
    ax.set_xlabel(""); ax.set_ylabel("")
    plt.tight_layout()
    fig.savefig(f"{save_prefix}_year_heatmap.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{save_prefix}_year_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ saved → {save_prefix}_year_heatmap.png / .pdf")


def round3(df):
    if df is None or df.empty:
        return df
    out  = df.copy()
    cols = out.select_dtypes(include=[np.number]).columns
    out[cols] = out[cols].round(3)
    return out


def build_combined_wide(all_yearly, value_col, out_name):
    frames = []
    for y in all_yearly:
        if y is None or y.empty:
            continue
        ds   = y["Dataset"].iloc[0]
        wide = (
            y.pivot(index="GeoCluster", columns="Year", values=value_col)
             .reindex(index=GEOCLUSTER_LABELS_ORDERED, columns=YEARS)
             .reset_index()
        )
        wide.insert(0, "Dataset", ds)
        frames.append(wide)
    if not frames:
        return
    combined = pd.concat(frames, ignore_index=True)
    combined.columns = ["Dataset", "GeoCluster"] + [str(c) for c in YEARS]
    round3(combined).to_csv(out_name, index=False)
    print(f"  ✓ combined wide table → {out_name}")


def main():
    all_overall, all_yearly = [], []

    for ds in DATASETS:
        print(f"\n==== {ds['label']} ====")
        tree_counter, posterior_prop, per_tree_rewards, per_tree_year_rewards = \
            parse_trees(ds["trees"])

        reward_df = make_reward_proportions(tree_counter, posterior_prop)
        plot_reward(reward_df, save_prefix=ds["out_prefix"])

        overall = round3(make_reward_table(per_tree_rewards, ds["label"]))
        overall.to_csv(f"{ds['out_prefix']}_reward_table.csv", index=False)
        all_overall.append(overall)

        yearly = round3(make_yearly_reward_table(per_tree_year_rewards, ds["label"]))
        yearly.to_csv(f"{ds['out_prefix']}_reward_by_year_long.csv", index=False)
        prop_wide = save_yearly_wide(yearly, ds["out_prefix"])
        plot_yearly_heatmap(prop_wide, ds["out_prefix"])
        all_yearly.append(yearly)

        print(f"  ✓ tables → {ds['out_prefix']}_reward_table.csv, "
              f"_reward_by_year_long.csv, _reward_by_year_prop_pct.csv, "
              f"_reward_by_year_dwell.csv")

    if all_overall:
        round3(pd.concat(all_overall, ignore_index=True)).to_csv(
            "GeoCluster_Reward_table_all_datasets.csv", index=False)
    if all_yearly:
        round3(pd.concat(all_yearly, ignore_index=True)).to_csv(
            "GeoCluster_Reward_by_year_all_datasets.csv", index=False)
        build_combined_wide(all_yearly, "reward_mean_yr",
                            "GeoCluster_Reward_by_year_dwell_all_datasets.csv")
        build_combined_wide(all_yearly, "prop_pct_mean",
                            "GeoCluster_Reward_by_year_prop_pct_all_datasets.csv")
    print("\n  ✓ combined tables written for all datasets")


if __name__ == "__main__":
    main()
