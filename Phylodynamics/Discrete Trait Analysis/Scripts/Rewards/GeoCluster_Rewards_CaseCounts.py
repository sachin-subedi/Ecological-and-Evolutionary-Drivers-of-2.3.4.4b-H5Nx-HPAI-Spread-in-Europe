from datetime import datetime
from io import StringIO
from pathlib import Path
import time

import pandas as pd
import matplotlib.pyplot as plt
from dateutil import parser
import baltic as bt
import matplotlib.colors as mcolors

# ── Inputs ─────────────────────────────────────────────────────────────────
RAW_CSV         = "Latest Reported Events.csv"
ANCHOR_DEC_YEAR = 2025.2520547945205

GEOCLUSTERS = [
    "GeoCluster_One", "GeoCluster_Two", "GeoCluster_Three",
    "GeoCluster_Four", "GeoCluster_Five"
]

GEOCLUSTER_COLORS = {
    "GeoCluster_One"   : "#4DBBD5",
    "GeoCluster_Two"   : "#F39B7F",
    "GeoCluster_Three" : "#3C5488",
    "GeoCluster_Four"  : "#BCBD22",
    "GeoCluster_Five"  : "#1B9E77",
}

def lighten_color(color, amount=0.2):
    try:
        c = mcolors.cnames[color]
    except KeyError:
        c = color
    r, g, b = mcolors.to_rgb(c)
    r_new = (1 - amount) * r + amount * 1.0
    g_new = (1 - amount) * g + amount * 1.0
    b_new = (1 - amount) * b + amount * 1.0
    return (r_new, g_new, b_new)

GEOCLUSTER_COLORS_LIGHT = {
    k: lighten_color(v, amount=0.2) for k, v in GEOCLUSTER_COLORS.items()
}

CLUSTER_LOOKUP = {
    "Albania": "GeoCluster_One",
    "Austria": "GeoCluster_Five",
    "Belgium": "GeoCluster_Three",
    "Bosnia and Herzegovina": "GeoCluster_Five",
    "Bulgaria": "GeoCluster_One",
    "Croatia": "GeoCluster_Five",
    "Cyprus": "GeoCluster_One",
    "Czechia": "GeoCluster_Five",
    "Denmark": "GeoCluster_Three",
    "Estonia": "GeoCluster_Four",
    "Finland": "GeoCluster_Four",
    "France": "GeoCluster_Three",
    "Germany": "GeoCluster_Three",
    "Greece": "GeoCluster_One",
    "Hungary": "GeoCluster_Two",
    "Iceland": "GeoCluster_Three",
    "Ireland": "GeoCluster_Three",
    "Italy": "GeoCluster_Five",
    "Kosovo": "GeoCluster_One",
    "Latvia": "GeoCluster_Four",
    "Lithuania": "GeoCluster_Four",
    "Luxembourg": "GeoCluster_Three",
    "Moldova": "GeoCluster_Two",
    "Netherlands": "GeoCluster_Three",
    "North Macedonia": "GeoCluster_One",
    "Norway": "GeoCluster_Four",
    "Poland": "GeoCluster_Two",
    "Portugal": "GeoCluster_Three",
    "Romania": "GeoCluster_One",
    "Serbia": "GeoCluster_One",
    "Slovakia": "GeoCluster_Two",
    "Slovenia": "GeoCluster_Five",
    "Spain": "GeoCluster_Three",
    "Sweden": "GeoCluster_Four",
    "Switzerland": "GeoCluster_Three",
    "Ukraine": "GeoCluster_Two",
    "United Kingdom": "GeoCluster_Three",
}

START = pd.to_datetime("2016-01-01")
END   = pd.to_datetime("2025-04-03")
ALL_MONTHS = pd.date_range(START, END, freq="MS")
YEAR_TICKS = pd.date_range(START, END, freq="YS")

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

def to_month_start(series_like):
    s = pd.to_datetime(series_like, errors="coerce")
    if getattr(s.dt, "tz", None) is not None:
        s = s.dt.tz_localize(None)
    return s.dt.to_period("M").dt.to_timestamp()

def decyear_to_yearmonth(dec):
    yr = int(dec)
    frac = float(dec) - yr
    dt = datetime(yr, 1, 1) + pd.Timedelta(days=frac * 365.25)
    return dt.strftime("%Y-%m")

def normalize_rows(df, columns):
    row_sums = df[columns].sum(axis=1)
    scaled = df[columns].div(row_sums.replace(0, pd.NA), axis=0)
    scaled = scaled.fillna(0.0)
    df_out = df.copy()
    df_out[columns] = scaled
    return df_out

def make_case_counts():
    cases_raw = pd.read_csv(RAW_CSV)

    if "observation date" not in cases_raw.columns:
        raise KeyError("Column 'observation date' not found in CSV.")
    if "Country" not in cases_raw.columns:
        raise KeyError("Column 'Country' not found in CSV.")

    cases_raw["observation_date"] = cases_raw["observation date"].apply(parse_date_safe)
    cases = cases_raw.dropna(subset=["observation_date"]).copy()
    cases["Country"] = cases["Country"].astype(str)

    cases["GeoCluster"] = cases["Country"].map(CLUSTER_LOOKUP)
    cases = cases.dropna(subset=["GeoCluster"]).copy()

    cases["YearMonth"] = to_month_start(cases["observation_date"])

    monthly_counts = (
        cases.groupby(["YearMonth", "GeoCluster"])
             .size()
             .reset_index(name="case_counts")
             .pivot_table(
                 index="YearMonth",
                 columns="GeoCluster",
                 values="case_counts",
                 fill_value=0
             )
             .reindex(ALL_MONTHS)
             .fillna(0)
    )

    for gc in GEOCLUSTERS:
        if gc not in monthly_counts.columns:
            monthly_counts[gc] = 0
    monthly_counts = monthly_counts[GEOCLUSTERS]

    return monthly_counts

def make_reward_proportions(trees_file):
    if not Path(trees_file).exists():
        raise FileNotFoundError(f"Cannot find trees file: {trees_file}")

    start_t = time.time()
    tree_counter = 0
    posterior_prop = {}

    with open(trees_file) as fh:
        for ln in fh:
            if not ln.lower().startswith("tree state_"):
                continue

            try:
                newick = ln.split("= [&R] ")[1]
                tree = bt.loadNewick(StringIO(newick), absoluteTime=False)
                tree.setAbsoluteTime(ANCHOR_DEC_YEAR)
            except Exception:
                continue

            tree_counter += 1

            month_cluster_raw = {}
            month_totals_raw  = {}

            for node in tree.Objects:
                if not getattr(node, "traits", None):
                    continue

                month = decyear_to_yearmonth(node.absoluteTime)
                node_rewards = {}
                for gc in GEOCLUSTERS:
                    key = f"{gc}_reward"
                    val = node.traits.get(key)
                    if val is None:
                        continue
                    try:
                        v = float(val)
                    except Exception:
                        continue
                    node_rewards[gc] = v
                    kpair = (month, gc)
                    month_cluster_raw[kpair] = month_cluster_raw.get(kpair, 0.0) + v
                try:
                    node_total = float(node.traits.get("GeoCluster.count", 0.0))
                except Exception:
                    node_total = 0.0
                if node_total == 0.0:
                    node_total = sum(node_rewards.values())

                month_totals_raw[month] = month_totals_raw.get(month, 0.0) + node_total
            for (mm, gc), raw_val in month_cluster_raw.items():
                total = month_totals_raw.get(mm, 0.0)
                if total <= 0.0:
                    continue
                prop = raw_val / total
                posterior_prop[(mm, gc)] = posterior_prop.get((mm, gc), 0.0) + prop

    if tree_counter <= 0:
        raise RuntimeError(f"No valid trees were parsed from {trees_file}")

    elapsed = time.time() - start_t
    print(f"✓ processed {tree_counter:,} trees from {trees_file} in {elapsed:.1f}s")

    rows = []
    for (mm, gc), summed in posterior_prop.items():
        rows.append({
            "YearMonth": pd.to_datetime(mm + "-01"),
            "GeoCluster": gc,
            "prop": summed / tree_counter
        })

    if not rows:
        out = pd.DataFrame(index=ALL_MONTHS, columns=GEOCLUSTERS).fillna(0.0)
        return out

    df = pd.DataFrame(rows)
    mat = (
        df.pivot_table(
            index="YearMonth",
            columns="GeoCluster",
            values="prop",
            fill_value=0.0
        )
        .reindex(ALL_MONTHS)
        .fillna(0.0)
    )

    for gc in GEOCLUSTERS:
        if gc not in mat.columns:
            mat[gc] = 0.0
    mat = mat[GEOCLUSTERS]

    mat = normalize_rows(mat, GEOCLUSTERS)

    return mat

def plot_combined(case_counts_df, reward_prop_df, save_prefix="GeoCluster_Proportions_SharedX"):
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(14, 9), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0]}
    )

    colors_light = [GEOCLUSTER_COLORS_LIGHT[g] for g in GEOCLUSTERS]

    case_counts_df[GEOCLUSTERS].plot.area(
        ax=ax_top,
        linewidth=0,
        color=colors_light
    )
    ax_top.set_ylabel("Case counts", fontsize=22)

    reward_prop_df[GEOCLUSTERS].plot.area(
        ax=ax_bot,
        linewidth=0,
        color=colors_light
    )
    ax_bot.set_ylabel("Reward proportion", fontsize=22)

    for ax in (ax_top, ax_bot):
        ax.set_facecolor("white")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.2)
        ax.spines["left"].set_linewidth(1.2)
        ax.tick_params(axis="y", labelsize=20)
        
    ax_bot.set_ylim(0, 1)


    ax_bot.set_xlim(START, END)
    ax_bot.set_xticks(YEAR_TICKS)
    ax_bot.set_xticklabels([d.year for d in YEAR_TICKS], fontsize=22)
    ax_bot.set_xlabel("")

    if ax_top.get_legend() is not None:
        ax_top.get_legend().remove()
    if ax_bot.get_legend() is not None:
        ax_bot.get_legend().remove()


    plt.tight_layout()
    fig.savefig(f"{save_prefix}.png", dpi=600)
    fig.savefig(f"{save_prefix}.pdf")
    plt.close(fig)

def main():
    case_df = make_case_counts()

    subsamples = [
        dict(
            tag="Subsample1",
            trees_file=(
                "/scratch/ss11645/old_combined/subsampled1/DTA/"
                "Habitat_Host_GC/combined_history_trees/GeoCluster/"
                "GeoCluster_history_Subsample1_combined.trees"
            ),
        ),
        dict(
            tag="Subsample2",
            trees_file=(
                "/scratch/ss11645/old_combined/subsampled2/DTA/"
                "Habitat_Host_GC/combined_history_trees/GeoCluster/"
                "GeoCluster_history_Subsample2_combined.trees"
            ),
        ),
        dict(
            tag="Subsample3",
            trees_file=(
                "/scratch/ss11645/old_combined/subsampled3/DTA/"
                "Habitat_Host_GC/combined_history_trees/GeoCluster/"
                "GeoCluster_history_Subsample3_combined.trees"
            ),
        ),
    ]

    for cfg in subsamples:
        tag        = cfg["tag"]
        trees_file = cfg["trees_file"]
        print(f"\n==== Processing GeoCluster rewards for {tag} ====")
        reward_df = make_reward_proportions(trees_file)
        prefix    = f"GeoCluster_Proportions_SharedX_{tag}"
        plot_combined(case_df, reward_df, save_prefix=prefix)
        print(f"Saved {prefix}.png / {prefix}.pdf")

if __name__ == "__main__":
    main()
