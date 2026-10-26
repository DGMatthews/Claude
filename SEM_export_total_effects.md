# Total effects of every SEM variable on performance

`total_effects.csv` gives, for every variable in the suction-feeding path model (SEM), its total effect on each performance outcome:
- `Max_sub_pressure`: peak sub-ambient pressure, Pa (negative = suction);
- `Total_ingested_volume`: water volume ingested, mm³.

It is for combining with the QTL results, to predict how much each QTL affects performance through the traits it changes.

`export_info.md` in the same folder names the model. The QTL traits are in `QTL_traits_per_fish.csv` (see `SEM_export_QTL_traits.md`).

## Where the numbers come from

- **Fitting:** `5_SEM_path_analysis.R` fits the model (see `SEM_export_path.md`).
- **Direct effects:** each arrow's direct effect is the slope of the child's equation with respect to the parent, at a reference point. The reference point is every variable at 90% of its largest value, except the four-bar inputs, which are at their medians. Identities (exact calculations) are linearized the same way.
- **Total effects:** the total effect of x on an outcome is the sum, over every path from x to the outcome, of the products of the direct effects along the path. With all direct effects in a matrix B, the totals are (I − B)⁻¹ − I.
- **Writing the file:** `6_SEM_exports.R` writes it from script 5's saved results.

## Columns

| column | meaning |
|---|---|
| `outcome`, `variable` | the performance outcome, and the variable whose effect is given |
| `role`, `group_guess` | as in the other exports (exogenous, fitted, identity, latent; morphology, kinematics, ...) |
| `has_path` | FALSE if there is no path from the variable to this outcome (all effects 0) |
| `direct`, `indirect`, `total` | effects in raw units: outcome units per unit of the variable. `indirect` = `total` − `direct` |
| `total_std` | `total` × SD(variable) / SD(outcome), with SDs across strikes (as script 5 reports) |
| `total_std_low`, `total_std_high` | bootstrap percentile interval of `total_std` over fish, if bootstrapped |
| `total_std_fish` | `total` × SD(variable) / SD(outcome), with SDs across fish of the per-fish values in `QTL_traits_per_fish.csv` (by default each fish's maximum over its strikes): the scale of a QTL analysis on those values |
| `sd_variable_*`, `sd_outcome_*` | the SDs used (across strikes, and across fish of the per-fish values) |
| `var_share_multiplier`, `var_share_multiplier_fish` | `total_std`², `total_std_fish`² (see below) |

## Combining with QTL results

Use the fish versions (`*_fish`) when the QTL were mapped on `QTL_traits_per_fish.csv`, whose values are by default each fish's maximum over its strikes.
- The direct and total effects themselves come from the strike-level model, and are slopes, so they apply to maxima as well.
- Only the SDs, and so the standardized values, are taken from the per-fish values.

**Signed effect of an allele.** If a QTL shifts trait x by `a` (its additive effect in raw units of x), the predicted shift in the outcome is `a × total` (outcome units). In SD units: `a_std × total_std_fish`, where `a_std = a / sd_variable_fish`.

**Share of the outcome's variance.** If the QTL explains `PVE_x` of the variance of x, the predicted share of the outcome's variance due to that QTL, acting through x, is:

> `PVE_outcome ≈ PVE_x × total_std_fish²` (= `PVE_x × var_share_multiplier_fish`)

This uses the **squared** standardized effect, not `PVE_x × total_std_fish`:
- PVE is a share of variance, and the standardized effect is a ratio of SDs, so it has to be squared to become a ratio of variances.
- Example: x explaining 10% of trait variance, with total_std_fish = 0.5, gives 0.10 × 0.25 = 2.5% of the outcome's variance.

### Cautions

1. **Don't add up traits on the same chain.** A trait's total effect already includes everything downstream of it.
   - If a QTL hits gape and also the flow that gape drives, adding gape's contribution to flow's counts the same effect twice.
   - Use the most upstream trait the QTL affects.
   - Or, for several traits, add the QTL's direct shifts δ to each trait (from a scan conditioned on the upstream traits). Then propagate them with the direct effects in `path_edges.csv`: Δ = (I − B)⁻¹ δ.
2. **A QTL can affect performance through paths the model doesn't contain.** These predictions are what the model's paths imply. They are not a substitute for mapping pressure and ingested volume directly. Comparing the two is a useful test of the model.
3. **The effects are linearized at the reference point.** Strong nonlinearities (products such as flow = volume rate / gape area, and the four-bar geometry) make them local.
   - Standardized total effects above 1, and so multipliers above 1, can occur, especially for `Total_ingested_volume`, which is a product (T × gape area × flow). Use them to rank traits.
   - Don't read multipliers above 1 as a literal share of variance.
4. **Measurement covariates are not biology.** `pressure_node_zR` and `static_node_zR_tmax` (how far the PIV node was from the mouth) have large effects on pressure, but they describe the measurement. Leave them out of QTL predictions.
5. **The bootstrap intervals are strike-level.** They are for `total_std`, not `total_std_fish`. They show which effects are distinguishable from 0, and the fish-level effects have about the same relative uncertainty.
6. **The Effort parts:** `Behavior_fish` (the fish's mean Effort) is the fish-level, heritable-looking part. `Effort_strike` is strike-to-strike variation, which isn't a fish trait; its per-fish maximum mostly reflects how many strikes the fish made.

## Partial effects by anatomical structure (`partial_effects_*.csv`)

These split each phenotype's total effect on each outcome by the parts of the feeding apparatus it acts through, for genome-wide performance maps of each structure separately.

### Definition

For a group G of variables (e.g. every hyoid variable), the partial effect of phenotype x on an outcome through G is the sum over all paths from x to the outcome that pass through at least one variable of G. It is calculated as:
- the total effect, minus
- the total effect left after cutting every arrow into G's variables.

If x itself belongs to G (e.g. `Hyoid_vel_tmax` for the hyoid group), all its effect goes through G.

### The groups (`partial_effect_groups.csv` lists every variable in each)

| group | variables matched by name |
|---|---|
| `neurocranium` | Neurocranium*, Epaxial_CSA, Supraoccipital_crest_depth |
| `hyoid` | Hyoid*, Sternohyoideus_CSA, Hypaxial_CSA, Urohyal_length, Ceratohyal_length |
| `premaxilla` | Protrusion*, Length_premaxilla_asc_arm |
| `maxilla` | Maxilla*, Oral_* (oral four-bar: mandible → maxilla → premaxilla) |
| `mandible` | Mandible*, Protractor_hyoideus_CSA |
| `opercular` | Opercular*, Operculum_length, IOP*, RA_process_length (opercular four-bar: neurocranium → mandible) |
| `gape` | R_*, A_* (gape radius and area; downstream of mandible and premaxilla) |
| `ram` | U_ram*, Cov_R_Uram |
| `timing` | T, ttpg, Time_*, tmax_from_Tstart, Relative_*, Kinetic_Synchronization |
| `effort` | Effort_strike, Behavior_fish |
| `none` | the paths that pass through none of the groups above (e.g. the PIV node distance acting directly on pressure) |

The groups are set by `partialGroups` in `6_SEM_exports.R` (regular expressions on variable names).

### Files

- **`partial_effects_by_group.csv`:** one row per outcome × phenotype × group (groups with no path omitted; `none` always listed). Its columns:
  - `partial`: raw units.
  - `partial_std`: standardized across strikes, with bootstrap interval `partial_std_low` / `partial_std_high`.
  - `partial_std_fish`: standardized with the per-fish SDs, as `total_std_fish`.
  - `share_of_total`: `partial / total`.
  - `var_share_multiplier_fish`: `partial_std_fish²`.
  - `variable_in_group`: TRUE if the phenotype itself belongs to the group.
- **`partial_effects_wide_<outcome>.csv`:** one row per phenotype, one column per group (`partial_std_fish`, 0 where there is no path), plus `total_std_fish`. Use it directly for the per-structure maps.

### Using them with the QTL results

For a QTL on phenotype x with allelic effect `a_std` (in SD units of x) and PVE_x:
- predicted shift in the outcome through structure G = `a_std × partial_std_fish(x, G)`;
- predicted share of the outcome's variance through G ≈ `PVE_x × partial_std_fish(x, G)²`.

A genome-wide map for structure G uses these values in place of the total effect. The same cautions apply as for total effects: use the most upstream trait a locus affects, and leave out the measurement covariates.

### Cautions specific to partial effects

1. **The groups overlap, so the partial effects don't add up to the total.** A path from jaw length through the mandible and then the gape counts in both `mandible` and `gape`. `none` plus the effect through any group equals the total, but the groups among themselves are not a partition.
2. **Paths can have opposite signs.** A group's partial effect can exceed the total effect, or have the opposite sign, when other paths cancel part of it. For example, protrusion velocity's effect through `gape` can be larger than its total.
3. **The bootstrap intervals are for the strike-level standardized values** (`partial_std`). Use them to judge which structure-specific effects are distinguishable from 0.
