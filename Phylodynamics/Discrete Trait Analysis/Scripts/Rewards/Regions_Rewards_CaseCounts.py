import time
from io import StringIO
from datetime import datetime
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from dateutil import parser
import baltic as bt

ANCHOR_DEC_YEAR = 2025.2520547945205
RAW_CSV         = "Latest Reported Events.xlsx"

RegionsS = [
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

Regions_COLORS = {
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

CLUSTER_LOOKUP = {
    "Albania"               : "HC_Cluster_2_Mediterranean",
    "Austria"               : "HC_Cluster_2_Alpine",
    "Belarus"               : "HC_Cluster_2_Alpine",
    "Belgium"               : "HC_Cluster_1_Atlantic",
    "Bosnia and Herz."      : "HC_Cluster_2_Alpine",
    "Bosnia and Herzegovina": "HC_Cluster_2_Alpine",
    "Bulgaria"              : "HC_Cluster_2_Continental",
    "Croatia"               : "HC_Cluster_2_Continental",
    "Cyprus"                : "HC_Cluster_2_Mediterranean",
    "Czechia"               : "HC_Cluster_2_Continental",
    "Denmark"               : "HC_Cluster_2_Continental",
    "Estonia"               : "HC_Cluster_3_Boreal",
    "Finland"               : "HC_Cluster_3_Boreal",
    "France"                : "HC_Cluster_1_Atlantic",
    "Germany"               : "HC_Cluster_1_Continental",
    "Greece"                : "HC_Cluster_2_Mediterranean",
    "Hungary"               : "HC_Cluster_2_Pannonian",
    "Iceland"               : "HC_Cluster_1_Atlantic",
    "Ireland"               : "HC_Cluster_1_Atlantic",
    "Italy"                 : "HC_Cluster_2_Mediterranean",
    "Kosovo"                : "HC_Cluster_2_Continental",
    "Latvia"                : "HC_Cluster_3_Boreal",
    "Lithuania"             : "HC_Cluster_3_Boreal",
    "Luxembourg"            : "HC_Cluster_1_Continental",
    "Moldova"               : "HC_Cluster_2_Continental",
    "Montenegro"            : "HC_Cluster_2_Continental",
    "Netherlands"           : "HC_Cluster_1_Atlantic",
    "North Macedonia"       : "HC_Cluster_2_Continental",
    "Norway"                : "HC_Cluster_3_Alpine",
    "Poland"                : "HC_Cluster_2_Continental",
    "Portugal"              : "HC_Cluster_4_Mediterranean",
    "Romania"               : "HC_Cluster_2_Continental",
    "Serbia"                : "HC_Cluster_2_Continental",
    "Republic of Serbia"    : "HC_Cluster_2_Continental",
    "Slovakia"              : "HC_Cluster_2_Alpine",
    "Slovenia"              : "HC_Cluster_2_Continental",
    "Spain"                 : "HC_Cluster_4_Mediterranean",
    "Sweden"                : "HC_Cluster_3_Boreal",
    "Switzerland"           : "HC_Cluster_1_Alpine",
    "Ukraine"               : "HC_Cluster_2_Continental",
    "United Kingdom"        : "HC_Cluster_1_Atlantic",
}

START      = pd.to_datetime("2016-01-01")
END        = pd.to_datetime("2025-04-03")
ALL_MONTHS = pd.date_range(START, END, freq="MS")
YEAR_TICKS = pd.date_range(START, END, freq="YS")

DATASETS = [
    {
        "label"  : "Equal",
        "trees"  : (
            "Regions_history_equal_combined.trees"
        ),
        "out_png": "Regions_CaseCounts_Reward_equal_no_legend",
    },
    {
        "label"  : "Proportional",
        "trees"  : (
            "Regions_history_proportional_combined.trees"
        ),
        "out_png": "Regions_CaseCounts_Reward_proportional_no_legend",
    },
    {
        "label"  : "Stratified",
        "trees"  : (
            "Regions_history_stratified_combined.trees"
        ),
        "out_png": "Regions_CaseCounts_Reward_stratified_no_legend",
    },
]


def parse_date_safe(x):
    if pd.isna(x):
        return pd.NaT
    for fmt in ("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y"):
        try:
            return datetime.strptime(str(x), fmt)
        except Exception:
            pass
    try:
        return parser.parse(str(x), dayfirst=False, fuzzy=True)
    except Exception:
        return pd.NaT


def to_month_start(s):
    s = pd.to_datetime(s, errors="coerce")
    if getattr(s.dt, "tz", None) is not None:
        s = s.dt.tz_localize(None)
    return s.dt.to_period("M").dt.to_timestamp()


def decyear_to_month(dec):
    yr   = int(dec)
    frac = dec - yr
    dt   = datetime(yr, 1, 1) + pd.Timedelta(days=frac * 365.25)
    return dt.strftime("%Y-%m")

def make_case_counts():
    cases_raw = pd.read_excel(RAW_CSV)
    cases_raw["observation_date"] = cases_raw["observation date"].apply(parse_date_safe)
    cases = cases_raw.dropna(subset=["observation_date"]).copy()
    cases["Country"]    = cases["Country"].astype(str)
    cases["Regions"] = cases["Country"].map(CLUSTER_LOOKUP)
    cases = cases.dropna(subset=["Regions"]).copy()
    cases["YearMonth"]  = to_month_start(cases["observation_date"])

    monthly_counts = (
        cases.groupby(["YearMonth", "Regions"])
             .size()
             .reset_index(name="case_counts")
             .pivot_table(
                 index="YearMonth",
                 columns="Regions",
                 values="case_counts",
                 fill_value=0,
             )
             .reindex(ALL_MONTHS)
             .fillna(0)
    )
    for gc in RegionsS:
        if gc not in monthly_counts.columns:
            monthly_counts[gc] = 0
    return monthly_counts[RegionsS]


def parse_trees(trees_file):
    if not Path(trees_file).exists():
        raise FileNotFoundError(f"Cannot find: {trees_file}")

    start        = time.time()
    tree_counter = 0
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
                for cl in RegionsS:
                    val = node.traits.get(f"{cl}_reward")
                    if val is None:
                        continue
                    node_rewards[cl] = float(val)
                    pair_key = (month, cl)
                    month_cluster_raw[pair_key] = (
                        month_cluster_raw.get(pair_key, 0.0) + float(val)
                    )
                node_total = float(node.traits.get("Regions.count", 0.0))
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
    return tree_counter, posterior_prop, all_months_set


def make_reward_proportions(tree_counter, posterior_prop):
    rows = []
    for (mm, gc), summed in posterior_prop.items():
        rows.append({
            "YearMonth" : pd.to_datetime(mm + "-01"),
            "Regions": gc,
            "prop"      : summed / tree_counter,
        })

    if not rows:
        return pd.DataFrame(index=ALL_MONTHS, columns=RegionsS).fillna(0.0)

    df  = pd.DataFrame(rows)
    mat = (
        df.pivot_table(
            index="YearMonth", columns="Regions",
            values="prop", fill_value=0.0,
        )
        .reindex(ALL_MONTHS)
        .fillna(0.0)
    )
    for gc in RegionsS:
        if gc not in mat.columns:
            mat[gc] = 0.0
    mat = mat[RegionsS]

    row_sums = mat.sum(axis=1)
    mat = mat.div(row_sums.replace(0, pd.NA), axis=0).fillna(0.0)

    zero_rows = mat.sum(axis=1) == 0
    mat.loc[zero_rows] = pd.NA
    mat = mat.ffill().fillna(0.0)

    return mat

def plot_combined(case_counts_df, reward_prop_df, save_prefix):
    colors = [Regions_COLORS[g] for g in RegionsS]

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(14, 9), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0]},
    )
    fig.patch.set_facecolor("white")

    ax_top.stackplot(ALL_MONTHS, [case_counts_df[g].values for g in RegionsS],
                     colors=colors, linewidth=0)
    ax_top.set_ylabel("Case counts", fontsize=22)

    ax_bot.stackplot(ALL_MONTHS, [reward_prop_df[g].values for g in RegionsS],
                     colors=colors, linewidth=0)
    ax_bot.set_ylabel("Reward proportion", fontsize=22)
    ax_bot.set_ylim(0, 1)

    for ax in (ax_top, ax_bot):
        ax.set_facecolor("white")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.2)
        ax.spines["left"].set_linewidth(1.2)
        ax.tick_params(axis="y", labelsize=20)

    ax_bot.set_xlim(START, END)
    ax_bot.set_xticks(YEAR_TICKS)
    ax_bot.set_xticklabels([d.year for d in YEAR_TICKS], fontsize=22)
    ax_bot.set_xlabel("")

    plt.tight_layout()
    fig.savefig(f"{save_prefix}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{save_prefix}.pdf", bbox_inches="tight")
    plt.show()
    plt.close(fig)
    print(f"  ✓ saved → {save_prefix}.png / .pdf")


def main():
    case_df = make_case_counts()

    for ds in DATASETS:
        print(f"\n==== {ds['label']} ====")
        tree_counter, posterior_prop, _ = parse_trees(ds["trees"])
        reward_df = make_reward_proportions(tree_counter, posterior_prop)
        plot_combined(case_df, reward_df, save_prefix=ds["out_png"])


if __name__ == "__main__":
    main()
