# QTL trait sheets from the SEM

These files hold every variable in the final suction-feeding path model (SEM), including the variables used only to calculate others (for example the four-bar link lengths behind the kinematic transmission coefficients). Each is a trait for QTL mapping in the ACxTRC hybrid cross.

`export_info.md` in the same folder names the model.

## Where the data come from

- **Kinematics** (gape, protrusion, hyoid, mandible and neurocranium movements, their velocities, timing, ram speed): `2_ACxTRC_Kinematic_Measurements_Sept2026.R`, from tracked landmarks in high-speed video, one row per strike.
  - Values called `*_tmax` are taken at the time of peak flow velocity (`tmaxVel`). `*_Max` and `*_t_mean` are the strike's maximum and mean.
- **Morphology** (muscle cross-sections, bone and link lengths, ranges of motion, four-bar rest angles, static KT): from the scans and the kinematics script. One value per fish.
- **Flow and pressure**: `3_PIV_Flow_Measurements.R`, from particle image velocimetry (PIV), corrected to the mouth aperture with a decay exponent k = 0.61.
  - `Max_sub_pressure` = peak sub-ambient pressure (Pa, negative = suction).
  - `U_flow_*` = flow speeds (mm/s).
  - `Total_ingested_volume` is in mm³.
  - Videos flagged as PIV problems in the data sheet have no flow or pressure.
- **Effort**: `Effort_score`, a latent score from `Neurocranium_Adj_vel_mean`, `Hyoid_Adj_vel_mean` and strike duration `T`. It is split into:
  - `Behavior_fish`: the fish's mean score, its consistent striking intensity;
  - `Effort_strike`: each strike's deviation from that mean.
- **Assembly:** `4_SEM_dataset.R` combines all of these. `5_SEM_path_analysis.R` fits the model and adds the Effort scores. `6_SEM_exports.R` writes these files.

**No size correction has been applied.** Every sheet includes the fish's standard length, `SL` in mm, for the size corrections. It comes from script 4 (`slColumn`: the scan file, else the data sheet). `SL` is not part of the SEM.

## Files

### `QTL_traits_per_fish.csv` (use this for mapping)

- **One row per fish** (`Fish_ID`), with `Family`, `Animal.ID`, `SL` and `n_strikes` (strikes recorded).
- **Each strike-level trait is the fish's maximum over its strikes:** the largest value in the trait's usual direction.
  - For traits that are mostly negative, that means the most negative value. For example, `Max_sub_pressure`, where more negative = stronger suction.
  - `max_direction` in the variable list says which applies.
- **Each trait's maximum can come from a different strike,** so a row describes each fish's capacity in every trait, not one particular strike.
- **Fish-level traits** (morphology) are the same in every strike, so the value is simply the value.
- **NA:** a trait is NA when the fish has no strike with that value. Coverage differs a lot between traits; see `n_fish_with_value` below. Pressure and flow exist for fewer fish than kinematics.
- **Setting:** this is `perFishSummary = "max"` in `6_SEM_exports.R`. With `"mean"`, the main file holds means instead.

### `QTL_traits_per_fish_mean.csv`

The same, with each trait's mean over the fish's strikes instead of its maximum, for comparison. If `perFishSummary = "mean"`, this file is `QTL_traits_per_fish_max.csv` instead.

### `QTL_traits_per_strike.csv`

- **One row per strike** (`vidName`), with the same IDs and traits.
- **Use it for** a mapping model with repeated measures (fish as a random effect), or to compute other summaries, such as a fish's best strike.

### `QTL_variable_list.csv`: one row per trait

| column | meaning |
|---|---|
| `variable` | the column name in the trait sheets (`SL` is listed with role `size (not in the SEM)`) |
| `role` | `exogenous` (no arrows in), `fitted` (has an equation), `identity` (calculated from other variables), `outcome` (performance), `Effort indicator`, `latent score` |
| `group_guess` | rough category from the name (morphology, kinematics, timing, buccal volume, flow, ram, effort, measurement covariate, performance) |
| `level` | `fish` = one value per fish (morphology); `strike` = varies between strikes |
| `n_strikes_with_value`, `n_fish_with_value` | coverage |
| `fish_value` | what `QTL_traits_per_fish.csv` holds for the trait: `value` (fish-level), `maximum over strikes` or `mean over strikes` |
| `max_direction` | `largest` or `most negative`: which end the maximum takes (traits whose median is negative take the most negative value) |
| `sd_strikes`, `sd_fish` | spread across strikes, and across fish (of the per-fish values) |
| `equation` | how the variable is predicted or calculated in the model (empty for exogenous variables) |

## How to use them

1. **Map every trait** in `QTL_traits_per_fish.csv`. Identity traits (e.g. `A_tmax`, the four-bar angles) are exact functions of other traits, so expect their QTL to overlap with the QTL of their inputs.
2. **Report each QTL peak's percentage of variance explained (PVE)** on the same per-fish scale (maxima), and its additive allelic effect, both raw and in SD units of the trait. The total-effects export (`total_effects.csv`, see `SEM_export_total_effects.md`) turns these into predicted effects on pressure and ingested volume.
3. **Strike-level traits carry measurement noise** (kinematics, flow, pressure).
   - A maximum over more strikes tends to be larger, and less noisy about the fish's capacity, so `n_strikes` can be used as a covariate or a weight.
   - The mean sheet is the alternative if the maxima look driven by single outlying strikes.
4. **Exclude the measurement covariates from biological interpretation.** `pressure_node_zR` and `static_node_zR_tmax` are how far the PIV measurement point was from the mouth, so they describe the measurement rather than the fish. Mapping them is a useful check for confounds.
