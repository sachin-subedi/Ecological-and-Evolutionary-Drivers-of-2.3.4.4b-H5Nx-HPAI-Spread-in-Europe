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

YEARS = list(range(2016, 2026))

DATASETS = [
    {
        "label"     : "Equal_rep1",
        "trees_file": (
            "Habitat_history_equal_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_equal_rep1",
    },
    {
        "label"     : "Equal_rep2",
        "trees_file": (
            "Habitat_history_equal_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_equal_rep2",
    },
    {
        "label"     : "Equal_rep3",
        "trees_file": (
            "Habitat_history_equal_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_equal_rep3",
    },
    {
        "label"     : "Proportional_rep1",
        "trees_file": (
            "Habitat_history_proportional_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_proportional_rep1",
    },
    {
        "label"     : "Proportional_rep2",
        "trees_file": (
            "Habitat_history_proportional_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_proportional_rep2",
    },
    {
        "label"     : "Proportional_rep3",
        "trees_file": (
            "Habitat_history_proportional_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_proportional_rep3",
    },
    {
        "label"     : "Stratified_rep1",
        "trees_file": (
            "Habitat_history_stratified_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_stratified_rep1",
    },
    {
        "label"     : "Stratified_rep2",
        "trees_file": (
            "Habitat_history_stratified_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_stratified_rep2",
    },
    {
        "label"     : "Stratified_rep3",
        "trees_file": (
            "Habitat_history_stratified_combined.trees"
        ),
        "out_prefix": "Habitat_Reward_stratified_rep3",
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

def round3(df):
    if df is None or df.empty:
        return df
    out  = df.copy()
    cols = out.select_dtypes(include=[np.number]).columns
    out[cols] = out[cols].round(3)
    return out

def make_reward_table(per_tree_rewards, label):
    if not per_tree_rewards:
        return pd.DataFrame()

    rdf            = pd.DataFrame(per_tree_rewards).reindex(columns=HABITATS).fillna(0.0)
    total_per_tree = rdf.sum(axis=1)
    safe_total     = total_per_tree.replace(0, np.nan)

    rows = []
    for hab in HABITATS:
        vals = rdf[hab].to_numpy(dtype=float)
        prop = (rdf[hab] / safe_total).to_numpy(dtype=float)
        r_hlo, r_hhi = hpd_interval(vals)
        r_clo, r_chi = eti_interval(vals)
        p_hlo, p_hhi = hpd_interval(prop)
        p_clo, p_chi = eti_interval(prop)
        rows.append({
            "Dataset"             : label,
            "Habitat"             : hab,
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
    return pd.DataFrame(rows)

def make_yearly_reward_table(per_tree_year_rewards, label):
    if not per_tree_year_rewards:
        return pd.DataFrame()

    dwell_samples = {(y, h): [] for y in YEARS for h in HABITATS}
    prop_samples  = {(y, h): [] for y in YEARS for h in HABITATS}

    for tree_yr in per_tree_year_rewards:
        year_tot = {y: 0.0 for y in YEARS}
        for (y, h), v in tree_yr.items():
            if y in year_tot:
                year_tot[y] += v
        for y in YEARS:
            tot = year_tot[y]
            for h in HABITATS:
                v = tree_yr.get((y, h), 0.0)
                dwell_samples[(y, h)].append(v)
                prop_samples[(y, h)].append(v / tot if tot > 0 else np.nan)

    rows = []
    for y in YEARS:
        for h in HABITATS:
            dv = np.asarray(dwell_samples[(y, h)], dtype=float)
            pv = np.asarray(prop_samples[(y, h)], dtype=float)
            r_hlo, r_hhi = hpd_interval(dv)
            r_clo, r_chi = eti_interval(dv)
            p_hlo, p_hhi = hpd_interval(pv)
            p_clo, p_chi = eti_interval(pv)
            rows.append({
                "Dataset"            : label,
                "Year"               : y,
                "Habitat"            : h,
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
    return pd.DataFrame(rows)

def save_yearly_wide(yearly_table, out_prefix):
    if yearly_table.empty:
        return None
    prop_wide = (
        yearly_table.pivot(index="Habitat", columns="Year", values="prop_pct_mean")
        .reindex(index=HABITATS, columns=YEARS)
    )
    dwell_wide = (
        yearly_table.pivot(index="Habitat", columns="Year", values="reward_mean_yr")
        .reindex(index=HABITATS, columns=YEARS)
    )
    prop_wide.round(3).to_csv(f"{out_prefix}_reward_by_year_prop_pct.csv")
    dwell_wide.round(3).to_csv(f"{out_prefix}_reward_by_year_dwell.csv")
    return prop_wide

def plot_yearly_heatmap(prop_wide, out_prefix):
    if prop_wide is None or prop_wide.empty:
        return
    mat = prop_wide.reindex(index=HABITATS, columns=YEARS).to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(11, 7))
    fig.patch.set_facecolor("white")

    vmax = np.nanmax(mat) if np.isfinite(np.nanmax(mat)) else 1
    im = ax.imshow(mat, aspect="auto", cmap="magma_r", vmin=0, vmax=vmax)

    ax.set_xticks(range(len(YEARS)))
    ax.set_xticklabels(YEARS, fontsize=13)
    ax.set_yticks(range(len(HABITATS)))
    ax.set_yticklabels([h.replace("_", " ") for h in HABITATS], fontsize=12)

    thr = vmax * 0.55
    for i in range(len(HABITATS)):
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
    fig.savefig(f"{out_prefix}_year_heatmap.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_year_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"✓ saved heatmap → {out_prefix}_year_heatmap.png / .pdf")

def make_habitat_reward_plot(trees_file, out_prefix):
    if not Path(trees_file).exists():
        raise FileNotFoundError(f"Cannot find trees file: {trees_file}")

    start = time.time()
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
            except Exception as e:
                print(f"Skipped tree {tree_counter} (parse error: {e})")
                continue

            month_habitat_raw  = {}
            month_totals_raw   = {}
            tree_reward_totals = {hab: 0.0 for hab in HABITATS}
            tree_year_reward   = {}

            for node in tree.Objects:
                if not node.traits:
                    continue

                month = decyear_to_month(node.absoluteTime)
                yr    = int(node.absoluteTime)
                all_months_set.add(month)

                node_rewards = {}
                for hab in HABITATS:
                    key = f"{hab}_reward"
                    val = node.traits.get(key)
                    if val is None:
                        continue
                    val = float(val)
                    node_rewards[hab] = val
                    tree_reward_totals[hab] += val
                    tree_year_reward[(yr, hab)] = tree_year_reward.get((yr, hab), 0.0) + val
                    pair_key = (month, hab)
                    month_habitat_raw[pair_key] = month_habitat_raw.get(pair_key, 0.0) + val

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

            per_tree_rewards.append(tree_reward_totals)
            per_tree_year_rewards.append(tree_year_reward)

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
    ax.set_xlabel("", fontsize=18)
    ax.set_ylabel("Reward proportion", fontsize=18)
    ax.tick_params(axis="y", labelsize=18)

    legend_labels = [
        "Coastal", "Farm", "Forest", "Grassland", "Human Modified",
        "Marine", "Shrubland", "Urban", "Woodland", "Rock", "Wetland"
    ]

    fig.subplots_adjust(top=0.80, bottom=0.12, left=0.10, right=0.98)

    handles, _ = ax.get_legend_handles_labels()
    ax.legend().remove()

    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        fontsize=12
    )

    fig.savefig(f"{out_prefix}.png", dpi=600, bbox_inches="tight")
    print(f"✓ saved plot to {out_prefix}.png")
    plt.show()
    plt.close(fig)

    return per_tree_rewards, per_tree_year_rewards

def build_combined_wide(all_yearly, value_col, out_name):
    frames = []
    for y in all_yearly:
        if y is None or y.empty:
            continue
        ds   = y["Dataset"].iloc[0]
        wide = (
            y.pivot(index="Habitat", columns="Year", values=value_col)
             .reindex(index=HABITATS, columns=YEARS)
             .reset_index()
        )
        wide.insert(0, "Dataset", ds)
        frames.append(wide)
    if not frames:
        return
    combined = pd.concat(frames, ignore_index=True)
    combined.columns = ["Dataset", "Habitat"] + [str(c) for c in YEARS]
    round3(combined).to_csv(out_name, index=False)
    print(f"✓ combined wide table → {out_name}")


def main():
    all_overall, all_yearly = [], []

    for ds in DATASETS:
        print(f"\n==== {ds['label']} ====")
        per_tree_rewards, per_tree_year_rewards = \
            make_habitat_reward_plot(ds["trees_file"], ds["out_prefix"])

        overall = round3(make_reward_table(per_tree_rewards, ds["label"]))
        overall.to_csv(f"{ds['out_prefix']}_reward_table.csv", index=False)
        all_overall.append(overall)

        yearly = round3(make_yearly_reward_table(per_tree_year_rewards, ds["label"]))
        yearly.to_csv(f"{ds['out_prefix']}_reward_by_year_long.csv", index=False)
        prop_wide = save_yearly_wide(yearly, ds["out_prefix"])
        plot_yearly_heatmap(prop_wide, ds["out_prefix"])
        all_yearly.append(yearly)

        print(f"✓ tables → {ds['out_prefix']}_reward_table.csv, "
              f"_reward_by_year_long.csv, _reward_by_year_prop_pct.csv, "
              f"_reward_by_year_dwell.csv")

    if all_overall:
        round3(pd.concat(all_overall, ignore_index=True)).to_csv(
            "Habitat_Reward_table_all_datasets.csv", index=False)
    if all_yearly:
        round3(pd.concat(all_yearly, ignore_index=True)).to_csv(
            "Habitat_Reward_by_year_all_datasets.csv", index=False)
        build_combined_wide(all_yearly, "reward_mean_yr",
                            "Habitat_Reward_by_year_dwell_all_datasets.csv")
        build_combined_wide(all_yearly, "prop_pct_mean",
                            "Habitat_Reward_by_year_prop_pct_all_datasets.csv")
    print("\n✓ combined tables written for all datasets")


if __name__ == "__main__":
    main()
