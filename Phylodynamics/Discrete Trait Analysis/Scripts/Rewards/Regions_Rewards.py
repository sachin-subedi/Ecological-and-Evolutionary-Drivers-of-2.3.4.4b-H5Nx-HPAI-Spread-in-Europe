import time
from io import StringIO
from datetime import datetime
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from dateutil import parser
import baltic as bt

ANCHOR_DEC_YEAR = 2025.2520547945205

RegionS = [
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

Region_COLORS = {
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

START      = pd.to_datetime("2016-01-01")
END        = pd.to_datetime("2025-04-03")
ALL_MONTHS = pd.date_range(START, END, freq="MS")
YEAR_TICKS = pd.date_range(START, END, freq="YS")

DATASETS = [
    {
        "label"  : "Equal",
        "trees"  : (
            "Region_history_equal_combined.trees"
        ),
        "out_png": "Region_Reward_equal_no_legend",
    },
    {
        "label"  : "Proportional",
        "trees"  : (
            "Region_history_proportional_combined.trees"
        ),
        "out_png": "Region_Reward_proportional_no_legend",
    },
    {
        "label"  : "Stratified",
        "trees"  : (
            "Region_history_stratified_combined.trees"
        ),
        "out_png": "Region_Reward_stratified_no_legend",
    },
]


def decyear_to_month(dec):
    yr   = int(dec)
    frac = dec - yr
    dt   = datetime(yr, 1, 1) + pd.Timedelta(days=frac * 365.25)
    return dt.strftime("%Y-%m")


def parse_trees(trees_file):
    if not Path(trees_file).exists():
        raise FileNotFoundError(f"Cannot find: {trees_file}")

    start          = time.time()
    tree_counter   = 0
    posterior_prop = {}
    all_months_set = set()

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

            month_cluster_raw = {}
            month_totals_raw  = {}

            for node in tree.Objects:
                if not node.traits:
                    continue
                month = decyear_to_month(node.absoluteTime)
                all_months_set.add(month)

                node_rewards = {}
                for cl in RegionS:
                    val = node.traits.get(f"{cl}_reward")
                    if val is None:
                        continue
                    node_rewards[cl] = float(val)
                    pair_key = (month, cl)
                    month_cluster_raw[pair_key] = (
                        month_cluster_raw.get(pair_key, 0.0) + float(val)
                    )
                node_total = float(node.traits.get("Region.count", 0.0))
                if node_total == 0:
                    node_total = sum(node_rewards.values())
                month_totals_raw[month] = month_totals_raw.get(month, 0.0) + node_total

            for (month, cl), raw_val in month_cluster_raw.items():
                total = month_totals_raw.get(month, 0.0)
                if total == 0:
                    continue
                prop = raw_val / total
                posterior_prop[(month, cl)] = posterior_prop.get((month, cl), 0.0) + prop

    print(f"  ✓ {tree_counter:,} trees in {time.time()-start:.1f}s")
    return tree_counter, posterior_prop

# ── Build reward df ───────────────────────────────────────────────────────────

def make_reward_proportions(tree_counter, posterior_prop):
    rows = []
    for (mm, gc), summed in posterior_prop.items():
        rows.append({
            "YearMonth" : pd.to_datetime(mm + "-01"),
            "Region": gc,
            "prop"      : summed / tree_counter,
        })

    if not rows:
        return pd.DataFrame(index=ALL_MONTHS, columns=RegionS).fillna(0.0)

    df  = pd.DataFrame(rows)
    mat = (
        df.pivot_table(
            index="YearMonth", columns="Region",
            values="prop", fill_value=0.0,
        )
        .reindex(ALL_MONTHS)
        .fillna(0.0)
    )
    for gc in RegionS:
        if gc not in mat.columns:
            mat[gc] = 0.0
    mat = mat[RegionS]

    row_sums = mat.sum(axis=1)
    mat = mat.div(row_sums.replace(0, pd.NA), axis=0).fillna(0.0)

    zero_rows = mat.sum(axis=1) == 0
    mat.loc[zero_rows] = pd.NA
    mat = mat.ffill().fillna(0.0)

    return mat


def plot_reward(reward_prop_df, save_prefix):
    colors = [Region_COLORS[g] for g in RegionS]

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("white")

    ax.stackplot(ALL_MONTHS, [reward_prop_df[g].values for g in RegionS],
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
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(f"{save_prefix}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{save_prefix}.pdf", bbox_inches="tight")
    plt.show()
    plt.close(fig)
    print(f"  ✓ saved → {save_prefix}.png / .pdf")


def main():
    for ds in DATASETS:
        print(f"\n==== {ds['label']} ====")
        tree_counter, posterior_prop = parse_trees(ds["trees"])
        reward_df = make_reward_proportions(tree_counter, posterior_prop)
        plot_reward(reward_df, save_prefix=ds["out_png"])


if __name__ == "__main__":
    main()
