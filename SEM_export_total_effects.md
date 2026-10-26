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
