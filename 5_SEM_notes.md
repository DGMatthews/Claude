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
3. `Effort_score ~ Prey_num + Satiation + Fasting_time` and Effort's arrows into `Hyoid_vel_tmax` and `Neurocranium_linear_vel_tmax` are ordinary path equations.
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
