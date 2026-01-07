import time, re
from io import StringIO
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import baltic as bt

ANCHOR_DEC_YEAR  = 2025.2520547945205

GEOCLUSTERS = [
    "GeoCluster_One",
    "GeoCluster_Two",
    "GeoCluster_Three",
    "GeoCluster_Four",
    "GeoCluster_Five"
]

GEOCLUSTER_COLORS = {
    "GeoCluster_One"   : "#4DBBD5",
    "GeoCluster_Two"   : "#F39B7F",
    "GeoCluster_Three" : "#3C5488",
    "GeoCluster_Four"  : "#BCBD22",
    "GeoCluster_Five"  : "#1B9E77",
}

def lighten_color(color, amount=0.5):
    r, g, b = mcolors.to_rgb(color)
    r = r + (1.0 - r) * amount
    g = g + (1.0 - g) * amount
    b = b + (1.0 - b) * amount
    return mcolors.to_hex((r, g, b))

def decyear_to_month(dec):
    yr   = int(dec)
    frac = dec - yr
    dt   = datetime(yr, 1, 1) + pd.Timedelta(days=frac * 365.25)
    return dt.strftime("%Y-%m")

def make_geocluster_reward_plot(trees_file, out_png):
    if not Path(trees_file).exists():
        raise FileNotFoundError(f"Cannot find trees file: {trees_file}")

    start = time.time()
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
            except Exception as e:
                print(f"Skipped tree {tree_counter} (parse error: {e})")
                continue

            month_cluster_raw = {}
            month_totals_raw  = {}

            for node in tree.Objects:
                if not node.traits:
                    continue

                month = decyear_to_month(node.absoluteTime)
                all_months_set.add(month)

                # gather this node’s GeoCluster rewards
                node_rewards = {}
                for cl in GEOCLUSTERS:
                    key = f"{cl}_reward"
                    val = node.traits.get(key)
                    if val is None:
                        continue
                    node_rewards[cl] = float(val)
                    pair_key = (month, cl)
                    month_cluster_raw[pair_key] = month_cluster_raw.get(pair_key, 0.0) + float(val)

                # denominator: use GeoCluster.count if >0, else sum of cluster rewards
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

    elapsed = time.time() - start
    print(f"✓ processed {tree_counter:,} trees from {trees_file} in {elapsed:.1f}s")

    mean_rows = []
    for (month, cl), summed_prop in posterior_prop.items():
        mean_rows.append({
            "year_month": month,
            "cluster": cl,
            "prop": summed_prop / tree_counter
        })

    mean_df = pd.DataFrame(mean_rows)
    full_months = sorted(all_months_set)
    full_index  = pd.MultiIndex.from_product(
        [full_months, GEOCLUSTERS], names=["year_month", "cluster"]
    )
    mean_df = (
        mean_df
        .set_index(["year_month", "cluster"])
        .reindex(full_index, fill_value=0.0)
        .reset_index()
    )

    plot_df = (
        mean_df
        .pivot_table(index="year_month", columns="cluster", values="prop", fill_value=0)
        .sort_index()
    )

    plot_df = plot_df.loc[(plot_df != 0).any(axis=1)]

    row_sums = plot_df.sum(axis=1)
    plot_df  = plot_df.div(row_sums, axis=0).fillna(0)

    plot_df.index = pd.to_datetime(plot_df.index + "-01")


    fig, ax = plt.subplots(figsize=(10, 6))

    plot_df[GEOCLUSTERS].plot.area(
        ax=ax,
        linewidth=0,
        color=[lighten_color(GEOCLUSTER_COLORS[c], amount=0.25) for c in GEOCLUSTERS],
        legend=False
    )

    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    ax.set_xlim(plot_df.index.min(), plot_df.index.max())
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)

    year_ticks = pd.date_range(plot_df.index.min(), plot_df.index.max(), freq="YS")
    ax.set_xticks(year_ticks)
    ax.set_xticklabels([d.year for d in year_ticks], fontsize=18)
    ax.set_xlabel("")
    ax.set_ylabel("Reward proportion", fontsize=18)
    ax.tick_params(axis="y", labelsize=18)

    # extra safety: remove any legend if created
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()

    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    print(f"✓ saved plot to {out_png}")
    plt.close(fig)


if __name__ == "__main__":
    subsamples = [
        dict(
            tag="Subsample1",
            trees_file=(
                "GeoCluster_history_Subsample1_combined.trees"
            ),
        ),
        dict(
            tag="Subsample2",
            trees_file=(
                "GeoCluster_history_Subsample2_combined.trees"
            ),
        ),
        dict(
            tag="Subsample3",
            trees_file=(
                "GeoCluster_history_Subsample3_combined.trees"
            ),
        ),
    ]

    for cfg in subsamples:
        tag        = cfg["tag"]
        trees_file = cfg["trees_file"]
        out_png    = f"GeoCluster_RewardProportion_{tag}.png"

        print(f"\n==== Processing {tag} ====")
        make_geocluster_reward_plot(trees_file, out_png)
