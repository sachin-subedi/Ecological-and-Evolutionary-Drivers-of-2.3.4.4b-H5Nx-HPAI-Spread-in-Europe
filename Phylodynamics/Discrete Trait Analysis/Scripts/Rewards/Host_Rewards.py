import time, re
from io import StringIO
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import baltic as bt

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial"]

ANCHOR_DEC_YEAR  = 2025.1068493150685

HOST = [
    "Domestic_Bird", "Domestic_Mammal",
    "Human", "Wild_Bird", "Wild_Mammal"
]

HOST_COLORS = {
  "Domestic_Bird"   : "#CD5C5C",
  "Domestic_Mammal" : "#6A3D9A",
  "Human"           : "#1B9E77",
  "Wild_Bird"       : "#3C5488",
  "Wild_Mammal"     : "#BCBD22",
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


def make_host_reward_plot(trees_file, out_png):
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

            # ---- load tree --------------------------------------------------
            try:
                newick = ln.split("= [&R] ")[1]
                tree   = bt.loadNewick(StringIO(newick), absoluteTime=False)
                tree.setAbsoluteTime(ANCHOR_DEC_YEAR)
            except Exception as e:
                print(f"Skipped tree {tree_counter} (parse error: {e})")
                continue

            month_host_raw   = {}
            month_totals_raw = {}

            for node in tree.Objects:
                if not node.traits:
                    continue

                month = decyear_to_month(node.absoluteTime)
                all_months_set.add(month)
                
                node_rewards = {}
                for h in HOST:
                    key = f"{h}_reward"
                    val = node.traits.get(key)
                    if val is None:
                        continue
                    node_rewards[h] = float(val)
                    pair_key = (month, h)
                    month_host_raw[pair_key] = month_host_raw.get(pair_key, 0.0) + float(val)

                node_total = float(node.traits.get("Host.count", 0.0))
                if node_total == 0:
                    node_total = sum(node_rewards.values())

                month_totals_raw[month] = month_totals_raw.get(month, 0.0) + node_total

            for (month, h), raw_val in month_host_raw.items():
                total = month_totals_raw.get(month, 0.0)
                if total == 0:
                    continue
                prop = raw_val / total
                posterior_prop[(month, h)] = posterior_prop.get((month, h), 0.0) + prop

    elapsed = time.time() - start
    print(f"✓ processed {tree_counter:,} trees from {trees_file} in {elapsed:.1f}s")

    mean_rows = []
    for (month, h), summed_prop in posterior_prop.items():
        mean_rows.append({
            "year_month": month,
            "host": h,
            "prop": summed_prop / tree_counter
        })

    mean_df = pd.DataFrame(mean_rows)
    full_months = sorted(all_months_set)
    full_index  = pd.MultiIndex.from_product(
                     [full_months, HOST], names=["year_month", "host"])
    mean_df = (
        mean_df
        .set_index(["year_month", "host"])
        .reindex(full_index, fill_value=0.0)
        .reset_index()
    )

    plot_df = (
        mean_df
        .pivot_table(
            index="year_month",
            columns="host",
            values="prop",
            fill_value=0
        )
        .sort_index()
    )
    plot_df = plot_df.loc[(plot_df != 0).any(axis=1)]


    row_sums = plot_df.sum(axis=1)
    plot_df  = plot_df.div(row_sums, axis=0).fillna(0)

    plot_df.index = pd.to_datetime(plot_df.index + "-01")

    fig, ax = plt.subplots(figsize=(10, 6))

    plot_df[HOST].plot.area(
        ax=ax,
        linewidth=0,
        color=[lighten_color(HOST_COLORS[h], amount=0.25) for h in HOST]
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
    ax.set_xlabel("", fontsize=18)
    ax.set_ylabel("Reward proportion", fontsize=18)
    ax.tick_params(axis="y", labelsize=18)

    legend_labels = [
        "Domestic Bird", "Domestic Mammal",
        "Human", "Wild Bird", "Wild Mammal"
    ]

    fig.subplots_adjust(top=0.80, bottom=0.12, left=0.10, right=0.98)

    handles, _ = ax.get_legend_handles_labels()
    ax.legend().remove()

    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        fontsize=12
    )

    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    print(f"✓ saved plot to {out_png}")
    plt.close(fig)

if __name__ == "__main__":
    subsamples = [
        dict(
            tag="Subsample1",
            trees_file=(
                "Host_history_Subsample1_combined.trees"
            ),
        ),
        dict(
            tag="Subsample2",
            trees_file=(
                "Host_history_Subsample2_combined.trees"
            ),
        ),
        dict(
            tag="Subsample3",
            trees_file=(
                "Host_history_Subsample3_combined.trees"
            ),
        ),
    ]

    for cfg in subsamples:
        tag        = cfg["tag"]
        trees_file = cfg["trees_file"]
        out_png    = f"Host_RewardProportion_{tag}.png"

        print(f"\n==== Processing {tag} ====")
        make_host_reward_plot(trees_file, out_png)
