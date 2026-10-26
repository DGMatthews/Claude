# Notes for script 5 (`5_SEM_path_analysis.R`)

Every step that involves an SEM is in script 5. Script 4 (`4_SEM_dataset.R`) only builds the dataset.

## Input

- Read `SEM_dataset.rds` (per strike) or `SEM_dataset_byFish.rds` (per fish), not the CSVs. The `.rds` keeps `Prey_num` as an ordered factor (one < few < several < many).
- Satiation = strike number (`_E##`), a proxy: higher strike numbers mean more food eaten since fasting began.
- No size correction anywhere in the SEM. Size corrections are only for the genetic screening of each parameter.

## Effort: two-step hybrid (in the script)

1. lavaan, measurement model only: `Effort =~ Neurocranium_Adj_vel_mean + Hyoid_Adj_vel_mean + T`, FIML (`missing = "ml"`), indicators standardized.
   - Just-identified with 3 indicators (df = 0): no fit test, so check the loadings.
2. **Bartlett factor scores** (`lavPredict(..., method = "Bartlett")`) = `Effort_score`. They are unbiased for the factor, so slopes on Effort aren't attenuated (regression scores would be).
3. Effort's arrows into `Hyoid_vel_tmax`, `Neurocranium_linear_vel_tmax` and the timing variables are ordinary path equations.
   - Effort is split into `Behavior_fish` (the fish mean) and `Effort_strike`; both are exogenous.
   - `Effort_score ~ Prey_num + Satiation + Fasting_time` is off by default (`effortCovariates`): these explained ~2% of strike effort.
4. The bootstrap redoes step 1 in every resample, so the error in the scores is in the confidence intervals.

## Decisions in the draft to check

- **Mandible check equation.** The docs have `Mandible_ang_predicted(tmax)` on the left. In the data, that column *is* the four-bar prediction (the change in θ4 from rest). So the script uses the measured `Mandible_ang_tmax` as the response, with the prediction as a predictor.
- **`Neurocranium_linear_vel_tmax ~ Epaxial_CSA * Effort * I(1/Neurocranium_MA)`.** The docs write `/Neurocranium_MA`. In an R formula `/` means nesting, so every ratio is written with `I()`.
- **`#` equations use R's `*`.** This gives all main effects and interactions, so the three-way products have 7 terms each. At fish level that's a lot of parameters for the sample size.
- **`*` (physics) equations use one product term (`I(...)`).** In the pressure equation each term is converted to Pa, so the coefficients should be about 1 if the physics holds.
- **Units of `Total_ingested_volume`.** The coefficients are 1 only if it is in mm³.
- **Linearization reference point.** Every variable is set to 90% of its largest-magnitude value. Identities are then recalculated from their parents' values, so the reference point is consistent (setting `refType`, with `"median"` as the alternative).
- **`Prey_num` in the path effects.** Its effect is the mean change per step up the order: (many − one) / 3.
- **T is an Effort indicator and also a predictor (exogenous in the path model).** The d-separation test will flag claims between Effort and the timing variables if Effort really drives T. If it does, a path from Effort to T, or to the timing variables, may be needed.
- **d-separation.**
  - Claims between two exogenous variables aren't tested; the morphology variables are correlated through size anyway.
  - Claims that identities make trivially true are skipped and listed.
  - With this many variables there are about 1000 claims, so Fisher's C is sensitive. The list of claims with the smallest p shows which arrows are missing.
- **Identities are checked against the data**, including the four-bar values, which use script 2's own functions. A row with a large `median_rel_diff` means script 2 or script 3 calculates the variable differently from the docs.

## Comparing models (end of script 5, `model_comparison.csv`)

Each row of `variants` is one model: time basis × width model × flow version × KT version × equation set (`structures`). Script 4 saves every time basis × width model × flow version combination in the `.rds`.

There are three kinds of fit, each comparable across a different set of models:

1. **Fisher's C** (absolute fit): p > 0.05 means the model isn't rejected. This works for any model, but C grows with the number of claims.
2. **AIC / BIC** (likelihood-based, Shipley & Douma 2020): the summed log-likelihoods of the path equations, each refitted by ML on the rows every model in the group shares.
   - Only valid between models whose responses are the same data (`aic_group`): KT version, width model, and structures that only change predictors.
   - The time basis and the flow version change the measured responses, so those models fall in different groups.
   - `model_comparison_equations.csv` shows which equations drive a difference.
3. **Cross-validated whole-model R²** of each outcome (predictive fit), on the rows every comparable model predicts (`cv_group_<outcome>`).
   - This is the comparison across groups. `Max_sub_pressure` is the same data in every model; `Total_ingested_volume` changes with the flow version.

Suggested order:
- Drop models with Fisher's C p < 0.05, if any pass.
- Within an `aic_group`, take the lowest AIC (ΔAIC < 2 = about equal).
- Across groups, take the highest cross-validated R² for both outcomes.

## Variance explained in each outcome

- **Local R²:** the outcome's own equation, using its measured direct predictors (e.g. pressure from the measured flow and gape).
- **Whole-model R²:** the outcome predicted from the exogenous variables only, through every equation in the chain, using fixed effects only. This is the share of variance the model as a whole explains.
- **Cross-validated whole-model R²:** the same, with each fold's fish predicted from models fitted without them. This is the honest number to report.

## Run time

The bootstrap refits everything (lavaan plus 21 mixed models) for every resample. Use `nBoot = 200` to check the script, and 1000+ for the final run, which will be slow (tens of minutes per time basis).

## Buccal width model

Only `noWidth` (no lateral expansion of the buccal cavity) is used: it is the default in scripts 4 and 5, and the only width model script 4 saves.
- Script 2: the bar model's expansion reached ~1.5 × W_head, and its retraction fit was poor (R² ≈ 0).
- Script 3: within videos, noWidth tracked the measured flow best (shape and volume rate). noRetract predicted about half the measured volume rate.
- Scripts 2 and 3 still calculate all three models as diagnostics. To compare them in the SEM, add `"main"` / `"noRetract"` to script 4's `widthModels` and set `compareWidthModels <- TRUE` in script 5.

## Strike maxima (`_Max`)

Each maximum has the same predictors as the variable at tmax, and takes the place of the tmax value as the predictor of the strike mean.

| Maximum | Equation (as at tmax) | Predicts |
|---|---|---|
| `R_Max` | `Protrusion_Max + Mandible_depression_Max` | `R_t_mean`, `R_tmax` |
| `Protrusion_Max` | `Neurocranium_tip_disp_Max + Maxilla_tip_disp_Max` | `Protrusion_t_mean`, `Protrusion_tmax`, `V_buccal_Max` |
| `Hyoid_Max` | `Relative_tmax * Hyoid_ROM` | `Hyoid_depression_t_mean`, `Hyoid_depression_tmax`, `V_buccal_Max` |
| `Neurocranium_Max` | `Relative_tmax * Neurocranium_ROM` | `Neurocranium_rotation_tmax` |
| `Mandible_depression_Max` | identity: `Mandible_length * sin(Mandible_Max)` | `Mandible_depression_t_mean`, `R_Max` |
| `A_Max` | identity: `(pi * W_jaw / 2) * R_Max` | `V_buccal_Max` |
| `Neurocranium_tip_disp_Max` | identity: `Neurocranium_Max * Neurocranium_length` | `Protrusion_Max` |

- **Arrows into the tmax values:** the Max → tmax arrows are `maxArrows`; `character(0)` removes them.
- **Exogenous inputs:** `Mandible_Max` and `Maxilla_tip_disp_Max` (from script 2's `Maxilla_tip_disp_peak`).
- **Built in script 4:** `R_Max` comes from `A_Max` (the spline peak where there is no `W_jaw`). Script 4 also calculates `Neurocranium_tip_disp_Max`.
- **Peak time basis:** the `_tmax` columns are the maxima, so the Max equations and arrows are dropped there. Every `_Max` name is replaced by its `_tmax` column (`maxToTmax`).

## Testing the Muller decay exponent

The flow correction is U_mouth = U_node × (1 + (z/R)²)^k, with Muller's k = 1.5.

1. Run script 3. The check "Muller exponent k estimated within videos" regresses log(node velocity / kinematic prediction) on log(1 + z²) within each video.
   - The prediction is the no-width buccal volume, in the earth frame.
   - Each video's scale is removed, so the estimate doesn't use differences between videos.
2. If `k_within` (ownLag) is clearly different from 1.5, set `mullerExponentAlt` to it in script 3 and rerun script 3. That adds the versions `staticK` and `mouthFrontCorrK`.
3. Rerun script 4. It adds the flow version `staticK+mouthFrontCorrK`, and skips it when script 3 didn't make it.
4. In script 5, set `compareMullerK <- TRUE`, with `compareDistance <- TRUE` for the same model without the distance covariates.
5. Compare the cross-validated pressure R², and check whether `static_node_zR_tmax` still matters.
   - If k was the problem, distance should matter much less with the new exponent.
   - If distance still matters as much, the blur / node-placement explanation is more likely.

The pressure correction g(s) is unchanged: it comes from Muller's potential, so it has no free exponent.

## Speed

Almost all the run time is refitting mixed models: about 2 refits per d-separation claim (~1600 claims per model), one refit of every equation per cross-validation fold, and one refit of everything per bootstrap sample.
- `nCores` runs the claims, folds and bootstrap samples in parallel. The default is all cores but one.
- `dsepVariants` limits the d-separation test to the models that need it; `character(0)` = none. It is the slowest part of a model without a bootstrap.
- `nBoot <- 0` gives no confidence intervals and leaves everything else the same.
- Bootstrapped models run last. `model_comparison_partial.csv` is saved after each model. The run time of each model and of the bootstrap is printed.
