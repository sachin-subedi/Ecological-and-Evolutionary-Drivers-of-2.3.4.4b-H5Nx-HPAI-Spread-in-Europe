from io      import StringIO
from pathlib import Path
import numpy  as np
import pandas as pd
import baltic as bt

TRAIT = "Habitat"

HABITATS_FULL = [
    "Coastal", "Farm", "Forest", "Grassland", "Human_Modified",
    "Marine", "Shrubland", "Urban", "Woodland", "Rock", "Wetland",
]

DATASETS = [
    {"tag": "equal_rep1",        "trees_violin": "bssvs_rates_equal_combined.trees"},
    {"tag": "equal_rep2",        "trees_violin": "bssvs_rates_equal_combined.trees"},
    {"tag": "equal_rep3",        "trees_violin": "bssvs_rates_equal_combined.trees"},
    {"tag": "proportional_rep1", "trees_violin": "bssvs_rates_proportional_combined.trees"},
    {"tag": "proportional_rep2", "trees_violin": "bssvs_rates_proportional_combined.trees"},
    {"tag": "proportional_rep3", "trees_violin": "bssvs_rates_proportional_combined.trees"},
    {"tag": "stratified_rep1",   "trees_violin": "bssvs_rates_stratified_combined.trees"},
    {"tag": "stratified_rep2",   "trees_violin": "bssvs_rates_stratified_combined.trees"},
    {"tag": "stratified_rep3",   "trees_violin": "bssvs_rates_stratified_combined.trees"},
]

ANCHOR_YEAR = 2025.2520547945205

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
        return (round(float(np.nanmean(vals)),        3),
                round(float(np.nanmedian(vals)),       3),
                round(float(np.nanpercentile(vals, 25)), 3),
                round(float(np.nanpercentile(vals, 75)), 3))

    rows = {"intro": [], "export": [], "pers": []}
    for state in HABITATS_FULL:                       # all 11 habitats, always
        for key, dct in [("intro", intro_dict), ("export", export_dict), ("pers", pers_dict)]:
            mean, median, q25, q75 = stats(dct.get(state, []))
            rows[key].append({"state": state, "Mean": mean, "Median": median,
                               "q25": q25, "q75": q75})
    return (pd.DataFrame(rows["intro"]),
            pd.DataFrame(rows["export"]),
            pd.DataFrame(rows["pers"]))

def _make_tables(trees_violin, out_name):
    if not Path(trees_violin).exists():
        raise FileNotFoundError(trees_violin)

    print(f"  Building violin data from {trees_violin} ...")
    intro_dict, export_dict, pers_dict = build_violin_dicts(
        trees_violin, ANCHOR_YEAR, trait=TRAIT)

    intro_df, export_df, pers_df = summarize_violin_stats_split(
        intro_dict, export_dict, pers_dict)

    intro_df.to_csv(f"{out_name}_intro_summary.csv",      index=False, float_format="%.3f")
    export_df.to_csv(f"{out_name}_export_summary.csv",    index=False, float_format="%.3f")
    pers_df.to_csv(f"{out_name}_persistence_summary.csv", index=False, float_format="%.3f")

    print(f"  ✓ saved → {out_name}_*.csv")
    return intro_df, export_df, pers_df


if __name__ == "__main__":
    all_intro, all_export, all_pers = [], [], []

    for ds in DATASETS:
        tag          = ds["tag"]
        trees_violin = ds["trees_violin"]

        scheme, rep = tag.rsplit("_", 1)   # e.g. "equal", "rep1"

        print(f"\n==== {tag.upper()} ====")
        intro_df, export_df, pers_df = _make_tables(
            trees_violin, out_name=f"Habitat_{tag}")

        for df, bucket in [(intro_df,  all_intro),
                           (export_df, all_export),
                           (pers_df,   all_pers)]:
            df = df.copy()
            df.insert(0, "tag",    tag)
            df.insert(1, "scheme", scheme)
            df.insert(2, "rep",    rep)
            bucket.append(df)

    if all_intro:
        intro_all  = pd.concat(all_intro,  ignore_index=True)
        export_all = pd.concat(all_export, ignore_index=True)
        pers_all   = pd.concat(all_pers,   ignore_index=True)

        intro_all .to_csv("Habitat_ALL_intro_summary.csv",       index=False, float_format="%.3f")
        export_all.to_csv("Habitat_ALL_export_summary.csv",      index=False, float_format="%.3f")
        pers_all  .to_csv("Habitat_ALL_persistence_summary.csv", index=False, float_format="%.3f")

        intro_all ["metric"] = "introductions"
        export_all["metric"] = "exportations"
        pers_all  ["metric"] = "persistence"
        master = pd.concat([intro_all, export_all, pers_all], ignore_index=True)
        master.to_csv("Habitat_ALL_combined_summary.csv", index=False, float_format="%.3f")
        print("\n  ✓ combined CSVs saved → Habitat_ALL_*.csv")
