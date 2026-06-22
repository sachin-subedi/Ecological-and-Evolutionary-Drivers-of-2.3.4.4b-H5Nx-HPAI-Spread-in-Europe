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

COLORS_FULL = {
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
HABITATS_FULL  = list(COLORS_FULL.keys())
LABEL_MAP_FULL = {k: k.replace("_", " ") for k in HABITATS_FULL}

# ── Grouped palette — MAIN FIGURE ─────────────────────────────────────────────
HABITAT_GROUP_MAP = {
    "Wetland"        : "WT",
    "Farm"           : "FA",
    "Grassland"      : "GW",
    "Coastal"        : "CM",
    "Marine"         : "CM",
    "Forest"         : "FO",
    "Woodland"       : "FO",
    "Shrubland"      : "FO",
    "Human_Modified" : "FO",
    "Rock"           : "FO",
    "Urban"          : "UB",
}

COLORS_GROUPED = {
    "WT" : "#006D2C",
    "FA" : "#CD5C5C",
    "GW" : "#26A69A",
    "CM" : "#0072B2",
    "FO" : "#BCBD22",
    "UB" : "#E6842A",
}

LABEL_GROUPED = {
    "WT" : "Wetland",
    "FA" : "Farm",
    "GW" : "Grassland",
    "CM" : "Coastal",
    "FO" : "Forest",
    "UB" : "Urban",
}

GROUPS_ORDER = ["CM", "FA", "FO", "GW", "UB", "WT"]

YEAR_MIN, YEAR_MAX = 2015.5, 2026
YEAR_TICKS         = np.arange(2016, 2026, 1)

# ── Datasets ──────────────────────────────────────────────────────────────────
DATASETS = [
    {
        "tag"        : "equal",
        "tree_full"  : "bssvs_rates_equal.tree",
        "trees_violin": "bssvs_rates/bssvs_rates_equal_combined.trees",
    },
    {
        "tag"        : "proportional",
        "tree_full"  : "bssvs_rates_proportional.tree",
        "trees_violin": "bssvs_rates_proportional_combined.trees",
    },
    {
        "tag"        : "stratified",
        "tree_full"  : "bssvs_rates_stratified.tree",
        "trees_violin": "bssvs_rates_stratified_combined.trees",
    },
]

ANCHOR_YEAR = 2025.2520547945205


def set_axis_font(ax):
    ax.tick_params(axis="both", which="major", labelsize=UNIFORM_FONTSIZE)
    ax.xaxis.label.set_size(UNIFORM_FONTSIZE)
    ax.yaxis.label.set_size(UNIFORM_FONTSIZE)

def _draw_tree_core(ax, tree_path, colors_dict, group_map=None,
                    legend_order=None, label_map=None, trait=TRAIT):
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
        raw_key   = k.traits[trait]
        state_key = group_map[raw_key] if (group_map and raw_key in group_map) else raw_key
        pp        = k.traits.get(f"{trait}.prob", 1.0)
        base_col  = colors_dict.get(state_key, "#B9B9B9")
        colour    = clr.LinearSegmentedColormap.from_list(
            "fade", ["#B9B9B9", base_col])(pp)

        x  = k.absoluteTime
        xp = k.parent.absoluteTime if k.parent else x
        y  = k.y

        if k.branchType == "leaf":
            ax.scatter(x, y, s=tip_size,       fc=colour, ec="none", zorder=11)
            ax.scatter(x, y, s=tip_size * 1.8, fc="k",    ec="none", zorder=10)
        elif k.branchType == "node":
            lw = base_lw + len(k.leaves) * 0.005
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

    if legend_order is None:
        legend_order = list(colors_dict.keys())
    if label_map is None:
        label_map = {k: k for k in colors_dict}

    handles = [
        mlines.Line2D([], [], color=colors_dict[k], marker="o",
                      markersize=14, markeredgewidth=2.2,
                      markerfacecolor=colors_dict[k],
                      label=label_map.get(k, k))
        for k in legend_order if k in colors_dict
    ]
    ax.legend(handles=handles,
              loc="center left", bbox_to_anchor=(0.15, 0.45),
              frameon=True, facecolor="white", edgecolor="white",
              prop={"size": UNIFORM_FONTSIZE})
    set_axis_font(ax)


def draw_tree_grouped(ax, tree_path, trait=TRAIT):
    _draw_tree_core(ax, tree_path,
                    colors_dict  = COLORS_GROUPED,
                    group_map    = HABITAT_GROUP_MAP,
                    legend_order = GROUPS_ORDER,
                    label_map    = LABEL_GROUPED,
                    trait        = trait)


def draw_tree_full(ax, tree_path, trait=TRAIT):
    _draw_tree_core(ax, tree_path,
                    colors_dict  = COLORS_FULL,
                    group_map    = None,
                    legend_order = HABITATS_FULL,
                    label_map    = LABEL_MAP_FULL,
                    trait        = trait)

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
                                 persistence_years=days / 365.25))
                        break
                    cur = cur.parent

    mig_df = pd.DataFrame(all_migrations)
    per_df = pd.DataFrame(all_persistences)

    intro = (mig_df.groupby(["tree_number", "child_state"]).size()
                   .reset_index(name="count"))
    expo  = (mig_df.groupby(["tree_number", "parent_state"]).size()
                   .reset_index(name="count"))
    pers  = (per_df.groupby(["tree_number", "child_state"])["persistence_years"]
                   .mean().reset_index())

    intro_dict  = {h: intro.loc[intro.child_state  == h, "count"].values
                   for h in HABITATS_FULL}
    export_dict = {h: expo .loc[expo .parent_state == h, "count"].values
                   for h in HABITATS_FULL}
    pers_dict   = {h: pers .loc[pers .child_state  == h, "persistence_years"].values
                   for h in HABITATS_FULL}
    return intro_dict, export_dict, pers_dict


def summarize_violin_stats_split(intro_dict, export_dict, pers_dict):
    def stats(vals):
        vals = np.asarray(vals, float)
        if vals.size == 0:
            return np.nan, np.nan, np.nan, np.nan
        return (float(np.nanmean(vals)), float(np.nanmedian(vals)),
                float(np.nanpercentile(vals, 25)), float(np.nanpercentile(vals, 75)))

    rows = {"intro": [], "export": [], "pers": []}
    for state in HABITATS_FULL:
        for key, dct in [("intro", intro_dict), ("export", export_dict), ("pers", pers_dict)]:
            mean, median, q25, q75 = stats(dct.get(state, []))
            rows[key].append({"state": state, "Mean": mean, "Median": median,
                               "q25": q25, "q75": q75})
    return (pd.DataFrame(rows["intro"]),
            pd.DataFrame(rows["export"]),
            pd.DataFrame(rows["pers"]))

def _draw_violin_panels_grouped(fig, right_gs, intro_dict, export_dict, pers_dict):
    ax_v1 = fig.add_subplot(right_gs[0, 0])
    ax_v2 = fig.add_subplot(right_gs[1, 0], sharex=ax_v1)
    ax_v3 = fig.add_subplot(right_gs[2, 0], sharex=ax_v1)

    def pool_groups(raw_dict):
        grouped = {g: [] for g in GROUPS_ORDER}
        for raw_hab, grp in HABITAT_GROUP_MAP.items():
            vals = raw_dict.get(raw_hab, [])
            if len(vals):
                grouped[grp].extend(vals.tolist())
        return {g: np.array(v) for g, v in grouped.items()}

    intro_grouped  = pool_groups(intro_dict)
    export_grouped = pool_groups(export_dict)
    pers_grouped   = pool_groups(pers_dict)

    panels = [
        ("Exportations",      export_grouped, ax_v1),
        ("Introductions",     intro_grouped,  ax_v2),
        ("Persistence (yrs)", pers_grouped,   ax_v3),
    ]
    positions   = list(range(len(GROUPS_ORDER)))
    xticklabels = [LABEL_GROUPED[g] for g in GROUPS_ORDER]

    for (ylab, data_dict, v_ax) in panels:
        for pos, grp in enumerate(GROUPS_ORDER):
            vals = data_dict.get(grp, np.array([]))
            if len(vals) == 0:
                continue
            vp   = v_ax.violinplot(vals, positions=[pos], widths=0.85,
                                   showmedians=True, bw_method=0.5,
                                   showextrema=False)
            body = vp["bodies"][0]
            body.set_facecolor(COLORS_GROUPED[grp])
            body.set_edgecolor(COLORS_GROUPED[grp])
            body.set_alpha(0.88)
            vp["cmedians"].set_edgecolor("black")

        v_ax.set_ylabel(ylab, fontsize=UNIFORM_FONTSIZE)
        v_ax.set_xticks(positions)
        v_ax.set_xticklabels(xticklabels, rotation=22, ha="right",
                             fontsize=UNIFORM_FONTSIZE - 1)
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


def _draw_violin_panels_full(fig, right_gs, intro_dict, export_dict, pers_dict):
    ax_v1 = fig.add_subplot(right_gs[0, 0])
    ax_v2 = fig.add_subplot(right_gs[1, 0], sharex=ax_v1)
    ax_v3 = fig.add_subplot(right_gs[2, 0], sharex=ax_v1)

    panels = [
        ("Exportations",      export_dict, ax_v1),
        ("Introductions",     intro_dict,  ax_v2),
        ("Persistence (yrs)", pers_dict,   ax_v3),
    ]
    positions   = list(range(len(HABITATS_FULL)))
    xticklabels = [LABEL_MAP_FULL[k] for k in HABITATS_FULL]

    for (ylab, data_dict, v_ax) in panels:
        for pos, state in enumerate(HABITATS_FULL):
            vals = data_dict.get(state, [])
            if len(vals) == 0:
                continue
            vp   = v_ax.violinplot(vals, positions=[pos], widths=0.85,
                                   showmedians=True, bw_method=0.5,
                                   showextrema=False)
            body = vp["bodies"][0]
            body.set_facecolor(COLORS_FULL[state])
            body.set_edgecolor(COLORS_FULL[state])
            body.set_alpha(0.88)
            vp["cmedians"].set_edgecolor("black")

        v_ax.set_ylabel(ylab, fontsize=UNIFORM_FONTSIZE)
        v_ax.set_xticks(positions)
        v_ax.set_xticklabels(xticklabels, rotation=22, ha="right",
                             fontsize=UNIFORM_FONTSIZE - 1)
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

def _make_figure(tree_full, trees_violin, out_name, grouped=True):
    for fp in [tree_full, trees_violin]:
        if not Path(fp).exists():
            raise FileNotFoundError(fp)

    print(f"  Building violin data from {trees_violin} ...")
    intro_dict, export_dict, pers_dict = build_violin_dicts(
        trees_violin, ANCHOR_YEAR, trait=TRAIT)

    if grouped:
        intro_df, export_df, pers_df = summarize_violin_stats_split(
            intro_dict, export_dict, pers_dict)
        intro_df.to_csv(f"{out_name}_intro_summary.csv",      index=False)
        export_df.to_csv(f"{out_name}_export_summary.csv",    index=False)
        pers_df.to_csv(f"{out_name}_persistence_summary.csv", index=False)

    fig     = plt.figure(figsize=(28, 16), facecolor="white")
    gs      = fig.add_gridspec(1, 2, width_ratios=[3.1, 2.6], wspace=0.06)
    ax_tree = fig.add_subplot(gs[0, 0])

    if grouped:
        draw_tree_grouped(ax_tree, tree_full, trait=TRAIT)
    else:
        draw_tree_full(ax_tree, tree_full, trait=TRAIT)

    right_gs = gs[0, 1].subgridspec(3, 1, hspace=0.04)
    if grouped:
        _draw_violin_panels_grouped(fig, right_gs, intro_dict, export_dict, pers_dict)
    else:
        _draw_violin_panels_full(fig, right_gs, intro_dict, export_dict, pers_dict)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.985, bottom=0.10,
                        wspace=0.06, hspace=0.00)
    plt.tight_layout(rect=[0.00, 0.00, 1.00, 0.98])
    fig.savefig(f"{out_name}.pdf")
    fig.savefig(f"{out_name}.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)
    print(f"  ✓ saved → {out_name}.png / .pdf")

if __name__ == "__main__":
    for ds in DATASETS:
        tag          = ds["tag"]
        tree_full    = ds["tree_full"]
        trees_violin = ds["trees_violin"]

        print(f"\n==== {tag.upper()} — main figure ====")
        _make_figure(tree_full, trees_violin,
                     out_name=f"Habitat_{tag}_main",
                     grouped=True)

        print(f"\n==== {tag.upper()} — supplementary figure ====")
        _make_figure(tree_full, trees_violin,
                     out_name=f"Habitat_{tag}_supplementary",
                     grouped=False)
