from io       import StringIO
from pathlib  import Path
import numpy  as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as clr
from   matplotlib.ticker import MaxNLocator
import baltic as bt

UNIFORM_FONTSIZE = 22
mpl.rcParams.update({
    "font.family"     : "sans-serif",
    "font.sans-serif" : ["DejaVu Sans"],
    "font.size"       : UNIFORM_FONTSIZE,
    "axes.linewidth"  : 1.2,
})

TRAIT = "Habitat"
COLORS = {
    "Coastal"        : "#0072B2",
    "Farm"           : "#CD5C5C",
    "Forest"         : "#BCBD22",
    "Grassland"      : "#26A69A",
    "Human_Modified" : "#56B4E9",
    "Marine"         : "#CC79A7",
    "Shrubland"      : "lightslategray",
    "Urban"          : "#E6842A",
    "Woodland"       : "#6B4C1B",
    "Rock"           : "#7570B3",
    "Wetland"        : "#006D2C",
}
HABITATS  = list(COLORS.keys())
LABEL_MAP = {k: k.replace("_", " ") for k in HABITATS}

YEAR_MIN, YEAR_MAX = 2015.5, 2026
YEAR_TICKS         = np.arange(2016, 2026, 1)

def set_axis_font(ax):
    ax.tick_params(axis='both', which='major', labelsize=UNIFORM_FONTSIZE)
    ax.xaxis.label.set_size(UNIFORM_FONTSIZE)
    ax.yaxis.label.set_size(UNIFORM_FONTSIZE)

def draw_tree(ax, tree_path, trait=TRAIT):
    tree = bt.loadNexus(tree_path, tip_regex=r'\|([0-9]+\-[0-9]+\-[0-9]+)')
    if tree.root.absoluteTime is None:
        latest_tip = max(
            (k for k in tree.Objects if k.branchType == "leaf"),
            key=lambda x: x.numdate
        )
        tree.setAbsoluteTime(latest_tip.numdate)

    base_lw, tip_size = 1.0, 15
    for k in tree.Objects:
        if not (hasattr(k, "traits") and trait in k.traits):
            continue
        key      = k.traits[trait]
        pp       = k.traits.get(f"{trait}.prob", 1.0)
        base_col = COLORS.get(key, "#B9B9B9")
        colour   = clr.LinearSegmentedColormap.from_list(
            "fade", ["#B9B9B9", base_col]
        )(pp)

        x  = k.absoluteTime
        xp = k.parent.absoluteTime if k.parent else x
        y  = k.y

        if k.branchType == "leaf":
            ax.scatter(x, y, s=tip_size,  fc=colour, ec="none", zorder=11)
            ax.scatter(x, y, s=tip_size*1.8, fc="k", ec="none", zorder=10)
        elif k.branchType == "node":
            lw = base_lw + len(k.leaves)*0.005
            ax.plot([x, x], [k.children[-1].y, k.children[0].y],
                    lw=lw, color=colour, zorder=9)
        ax.plot([xp, x], [y, y], lw=base_lw, color=colour, zorder=9)

    ax.set_xlim(YEAR_MIN, YEAR_MAX)
    ax.set_xticks(YEAR_TICKS)
    ax.tick_params(axis="x", which="major", labelbottom=True, length=6)
    ax.set_xlabel("", labelpad=8)

    ax.set_ylim(-5, tree.ySpan + 5)
    ax.set_yticks([])
    sns.despine(ax=ax, left=True, bottom=False)
    ax.grid(False)

    handles = [mlines.Line2D([], [], color=COLORS[k], marker="o",
                             markersize=14, markeredgewidth=2.2,
                             markerfacecolor=COLORS[k], label=LABEL_MAP[k])
               for k in HABITATS]
    ax.legend(handles=handles,
              loc="center left", bbox_to_anchor=(0.15, 0.45),
              frameon=True, facecolor="white", edgecolor="white",
              prop={"size": UNIFORM_FONTSIZE})
    set_axis_font(ax)

def build_violin_dicts(trees_file, anchor_year, trait=TRAIT):
    all_migrations, all_persistences = [], []

    with open(trees_file) as fh:
        for idx, raw in enumerate(fh, 1):
            if not raw.lower().startswith("tree state_"):
                continue
            try:
                newick = raw.split("= [&R] ")[1].strip()
            except IndexError:
                continue
            try:
                tree = bt.loadNewick(StringIO(newick), absoluteTime=False)
                tree.setAbsoluteTime(anchor_year)
            except Exception:
                continue

            for node in tree.Objects:
                if not node.traits:
                    continue
                child = node.traits.get(trait)
                par   = (node.parent.traits.get(trait)
                         if (node.parent and node.parent.traits) else "root")
                if child != par:
                    all_migrations.append(
                        dict(tree_number=idx, parent_state=par, child_state=child))

            for tip in tree.Objects:
                if tip.branchType != "leaf" or not tip.traits:
                    continue
                state = tip.traits.get(trait)
                cur   = tip.parent
                while cur and cur.traits:
                    par_state = cur.traits.get(trait)
                    if par_state != state:
                        days = (tip.absoluteTime - cur.absoluteTime) * 365.25
                        all_persistences.append(
                            dict(tree_number=idx, child_state=state,
                                 persistence_years=days/365.25))
                        break
                    cur = cur.parent

    mig_df = pd.DataFrame(all_migrations)
    per_df = pd.DataFrame(all_persistences)

    intro = (mig_df.groupby(["tree_number","child_state"]).size()
                  .reset_index(name="count"))
    expo  = (mig_df.groupby(["tree_number","parent_state"]).size()
                  .reset_index(name="count"))
    pers  = (per_df.groupby(["tree_number","child_state"])["persistence_years"]
                   .mean().reset_index())

    intro_dict  = {h: intro.loc[intro.child_state  == h, "count"].values
                   for h in HABITATS}
    export_dict = {h: expo .loc[expo .parent_state == h, "count"].values
                   for h in HABITATS}
    pers_dict   = {h: pers .loc[pers .child_state  == h, "persistence_years"].values
                   for h in HABITATS}
    return intro_dict, export_dict, pers_dict

def summarize_violin_stats_split(intro_dict, export_dict, pers_dict):
    """
    Build three DataFrames (Intro, Export, Persistence), each with:
      state, Mean, Median, q25, q75
    """
    def stats(vals):
        vals = np.asarray(vals, float)
        if vals.size == 0:
            return np.nan, np.nan, np.nan, np.nan
        mean   = float(np.nanmean(vals))
        median = float(np.nanmedian(vals))
        q25    = float(np.nanpercentile(vals, 25))
        q75    = float(np.nanpercentile(vals, 75))
        return mean, median, q25, q75

    intro_rows, export_rows, pers_rows = [], [], []

    for state in HABITATS:
        mean, median, q25, q75 = stats(intro_dict.get(state, []))
        intro_rows.append({
            "state": state,
            "Mean":  mean,
            "Median": median,
            "q25":  q25,
            "q75":  q75,
        })

        mean, median, q25, q75 = stats(export_dict.get(state, []))
        export_rows.append({
            "state": state,
            "Mean":  mean,
            "Median": median,
            "q25":  q25,
            "q75":  q75,
        })

        mean, median, q25, q75 = stats(pers_dict.get(state, []))
        pers_rows.append({
            "state": state,
            "Mean":  mean,
            "Median": median,
            "q25":  q25,
            "q75":  q75,
        })

    intro_df  = pd.DataFrame(intro_rows)
    export_df = pd.DataFrame(export_rows)
    pers_df   = pd.DataFrame(pers_rows)
    return intro_df, export_df, pers_df

def make_four_panel_sidebar(tree_full         ="bssvs_rates.tree",
                            trees_file_violin ="bssvs500.trees",
                            anchor_year       =2025.2520547945205,
                            out_name          ="Habitat_four_panel_sidebar"):
    for fp in [tree_full, trees_file_violin]:
        if not Path(fp).exists():
            raise FileNotFoundError(fp)

    intro_dict, export_dict, pers_dict = build_violin_dicts(
        trees_file_violin, anchor_year, trait=TRAIT
    )

    intro_df, export_df, pers_df = summarize_violin_stats_split(
        intro_dict, export_dict, pers_dict
    )
    intro_df.to_csv(f"{out_name}_intro_summary.csv",       index=False)
    export_df.to_csv(f"{out_name}_export_summary.csv",      index=False)
    pers_df.to_csv(f"{out_name}_persistence_summary.csv",   index=False)

    fig = plt.figure(figsize=(28, 16), facecolor="white")
    gs  = fig.add_gridspec(1, 2, width_ratios=[3.1, 2.6], wspace=0.06)

    ax_tree = fig.add_subplot(gs[0, 0])
    draw_tree(ax_tree, tree_full, trait=TRAIT)

    right_gs = gs[0, 1].subgridspec(3, 1, hspace=0.04)

    ax_v1 = fig.add_subplot(right_gs[0, 0])
    ax_v2 = fig.add_subplot(right_gs[1, 0], sharex=ax_v1)
    ax_v3 = fig.add_subplot(right_gs[2, 0], sharex=ax_v1)

    panels = [
        ("Exportations",      export_dict, ax_v1),
        ("Introductions",     intro_dict,  ax_v2),
        ("Persistence (yrs)", pers_dict,   ax_v3),
    ]

    positions   = list(range(len(HABITATS)))
    xticklabels = [LABEL_MAP[k] for k in HABITATS]

    for (ylab, data_dict, v_ax) in panels:
        for pos, state in enumerate(HABITATS):
            vals = data_dict.get(state, [])
            if len(vals) == 0:
                continue
            vp = v_ax.violinplot(vals, positions=[pos], widths=0.85,
                                 showmedians=True, bw_method=0.5,
                                 showextrema=False)
            body = vp["bodies"][0]
            body.set_facecolor(COLORS[state]); body.set_edgecolor(COLORS[state])
            body.set_alpha(0.88)
            vp["cmedians"].set_edgecolor("black")

        v_ax.set_ylabel(ylab, fontsize=UNIFORM_FONTSIZE)
        v_ax.set_xticks(positions)
        v_ax.set_xticklabels(xticklabels, rotation=22, ha="right",
                             fontsize=UNIFORM_FONTSIZE-1)

        v_ax.grid(False)
        sns.despine(ax=v_ax, left=False, bottom=False)
        if ylab in ("Exportations", "Introductions"):
            v_ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            v_ax.set_ylim(bottom=0)

        set_axis_font(v_ax)

    plt.setp(ax_v1.get_xticklabels(), visible=False)
    plt.setp(ax_v2.get_xticklabels(), visible=False)
    ax_v1.tick_params(axis="x", which="both", length=0)
    ax_v2.tick_params(axis="x", which="both", length=0)
    ax_v3.set_xlabel("", fontsize=UNIFORM_FONTSIZE)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.985, bottom=0.10,
                        wspace=0.06, hspace=0.00)

    plt.tight_layout(rect=[0.00, 0.00, 1.00, 0.98])
    fig.savefig(f"{out_name}.pdf")
    fig.savefig(f"{out_name}.png", dpi=300, bbox_inches="tight")
    plt.show()

    return intro_df, export_df, pers_df

def run_all_subsamples():
    """
    Run the four–panel sidebar figure for all three subsamples.

    Adjust the file paths below if your directory layout / filenames differ.
    """
    anchor_year = 2025.2520547945205

    subsamples = [
        dict(
            tag="Subsample1",
            tree_full=(
                "/scratch/ss11645/old_combined/subsampled1/DTA/"
                "Habitat_Host_GC/bssvs_rates/bssvs_rates_Subsample1.tree"
            ),
            trees_file_violin=(
                "/scratch/ss11645/old_combined/subsampled1/DTA/"
                "Habitat_Host_GC/bssvs_rates/"
                "bssvs_rates_Subsample1_combined.trees"
            ),
        ),
        dict(
            tag="Subsample2",
            tree_full=(
                "/scratch/ss11645/old_combined/subsampled2/DTA/"
                "Habitat_Host_GC/bssvs_rates/bssvs_rates_Subsample2.tree"
            ),
            trees_file_violin=(
                "/scratch/ss11645/old_combined/subsampled2/DTA/"
                "Habitat_Host_GC/bssvs_rates/"
                "bssvs_rates_Subsample2_combined.trees"
            ),
        ),
        dict(
            tag="Subsample3",
            tree_full=(
                "/scratch/ss11645/old_combined/subsampled3/DTA/"
                "Habitat_Host_GC/bssvs_rates/bssvs_rates_Subsample3.tree"
            ),
            trees_file_violin=(
                "/scratch/ss11645/old_combined/subsampled3/DTA/"
                "Habitat_Host_GC/bssvs_rates/"
                "bssvs_rates_Subsample3_combined.trees"
            ),
        ),
    ]

    for cfg in subsamples:
        tag = cfg["tag"]
        print(f"\n── Processing {tag} ─────────────────────────────────────")
        make_four_panel_sidebar(
            tree_full=cfg["tree_full"],
            trees_file_violin=cfg["trees_file_violin"],
            anchor_year=anchor_year,
            out_name=f"Habitat_four_panel_sidebar_{tag}",
        )

if __name__ == "__main__":
    run_all_subsamples()
