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

ANCHOR_DEC_YEAR = 2025.2520547945205

HABITATS = [
    "Coastal", "Farm", "Forest", "Grassland", "Human_Modified",
    "Marine", "Shrubland", "Urban", "Woodland", "Rock", "Wetland"
]

HABITAT_COLORS = {
    "Coastal":        "#0072B2",
    "Farm":           "#CD5C5C",
    "Forest":         "#BCBD22",
    "Grassland":      "#26A69A",
    "Human_Modified": "#56B4E9",
    "Marine":         "#CC79A7",
    "Shrubland":      "lightslategray",
    "Urban":          "#E6842A",
    "Woodland":       "#6B4C1B",
    "Rock":           "#7570B3",
    "Wetland":        "#006D2C",
}

DATASETS = [
    {
        "trees_file": (
            "Habitat_history_equal_combined.trees"
        ),
        "out_png": "Habitat_Reward_equal_NL.png",
    },
    {
        "trees_file": (
            "Habitat_history_proportional_combined.trees"
        ),
        "out_png": "Habitat_Reward_proportional_NL.png",
    },
    {
        "trees_file": (
            "Habitat_history_stratified_combined.trees"
        ),
        "out_png": "Habitat_Reward_stratified_NL.png",
    },
]


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

def make_habitat_reward_plot(trees_file, out_png):
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

            month_habitat_raw = {}
            month_totals_raw  = {}

            for node in tree.Objects:
                if not node.traits:
                    continue

                month = decyear_to_month(node.absoluteTime)
                all_months_set.add(month)

                node_rewards = {}
                for hab in HABITATS:
                    key = f"{hab}_reward"
                    val = node.traits.get(key)
                    if val is None:
                        continue
                    node_rewards[hab] = float(val)
                    pair_key = (month, hab)
                    month_habitat_raw[pair_key] = month_habitat_raw.get(pair_key, 0.0) + float(val)

                node_total = float(node.traits.get("Habitat.count", 0.0))
                if node_total == 0:
                    node_total = sum(node_rewards.values())

                month_totals_raw[month] = month_totals_raw.get(month, 0.0) + node_total

            for (month, hab), raw_val in month_habitat_raw.items():
                total = month_totals_raw.get(month, 0.0)
                if total == 0:
                    continue
                prop = raw_val / total
                posterior_prop[(month, hab)] = posterior_prop.get((month, hab), 0.0) + prop

    elapsed = time.time() - start
    print(f"✓ processed {tree_counter:,} trees from {trees_file} in {elapsed:.1f}s")

    mean_rows = []
    for (month, hab), summed_prop in posterior_prop.items():
        mean_rows.append({
            "year_month": month,
            "habitat":    hab,
            "prop":       summed_prop / tree_counter
        })

    mean_df = pd.DataFrame(mean_rows)
    full_months = sorted(all_months_set)
    full_index  = pd.MultiIndex.from_product(
                     [full_months, HABITATS], names=["year_month", "habitat"])
    mean_df = (
        mean_df
        .set_index(["year_month", "habitat"])
        .reindex(full_index, fill_value=0.0)
        .reset_index()
    )

    plot_df = (
        mean_df
        .pivot_table(
            index="year_month",
            columns="habitat",
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

    plot_df[HABITATS].plot.area(
        ax=ax,
        linewidth=0,
        color=[lighten_color(HABITAT_COLORS[h], amount=0.25) for h in HABITATS]
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
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_xlabel("", fontsize=18)
    ax.set_ylabel("Reward proportion", fontsize=18)
    ax.tick_params(axis="y", labelsize=18)

    ax.legend().remove()

    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    print(f"✓ saved plot to {out_png}")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    for ds in DATASETS:
        make_habitat_reward_plot(ds["trees_file"], ds["out_png"])
