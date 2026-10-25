# Notes for script 5 (SEM)

Every step that involves an SEM goes in script 5. Script 4 (`4_SEM_dataset.R`) only builds the dataset.

## Input

- Read `SEM_dataset.rds` (per strike) or `SEM_dataset_byFish.rds` (per fish), not the CSVs.
  - Each file is a list with `$tmax` and `$tmaxVel` (values at tmaxVel are saved under the `_tmax` names; the `timeBasis` column says which).
  - The `.rds` keeps `Prey_num` as an **ordered factor**. A CSV loses that.
- Prey number is ordered categorical. In `lm`/`lmer` an ordered factor gets polynomial contrasts (`.L`, `.Q`, ...).
  - If you want treatment contrasts instead, set them explicitly. Don't convert it to numeric.
- Satiation = strike number (`_E##`), used as a proxy: higher strike numbers mean more food eaten since fasting began.
- No size correction anywhere in the SEM. Size corrections are only for the genetic screening of each parameter.

## Effort (latent variable): two-step hybrid

piecewiseSEM has no latent variables, so Effort is made first in lavaan and then used as an observed variable.

1. **Measurement model in lavaan.** Only the measurement part goes in lavaan:
   ```r
   mEffort <- 'Effort =~ Neurocranium_Adj_vel_mean + Hyoid_Adj_vel_mean + T'
   fitEffort <- cfa(mEffort, data = d, missing = "ml", std.lv = TRUE)   # FIML: rows missing one indicator are kept
   ```
   - With 3 indicators the model is just-identified (df = 0): the loadings are estimated, but there is no fit test.
   - Check that the loadings have the expected signs and that none is near 0.
   - Standardize the indicators first if their scales differ a lot (T is in s, the velocities in mm/s).
2. **Factor scores: use Bartlett scores.**
   ```r
   d$Effort_score <- NA
   d$Effort_score[lavInspect(fitEffort, "case.idx")] <- lavPredict(fitEffort, method = "Bartlett")[, "Effort"]
   ```
   - Bartlett scores are unbiased for the factor, so their slopes are not biased when Effort is used as a predictor later. Regression scores bias those slopes.
   - Use `case.idx` so the scores line up with the right rows.
   - Set the sign so Effort increases with `Neurocranium_Adj_vel_mean`.
3. **piecewiseSEM with the scores** as an observed variable:
   - `Effort_score ~ Prey_num + Satiation + Fasting_time` (the covariates go here, not in the lavaan step)
   - `Effort_score` as a predictor in the `Hyoid_vel_tmax` and `Neurocranium_linear_vel_tmax` equations
   - Per strike: use `lmer(... + (1 | Fish_ID))` in every equation.
4. **Standard errors.** The scores are estimates, and step 3 treats them as known, so its SEs are too small.
   - Bootstrap the whole two-step procedure: refit steps 1–3 on each resample.
   - Per strike, resample fish (`Fish_ID`), not strikes.

## Datasets

- Fit the model at both time bases (`$tmax` and `$tmaxVel`) and compare.
- Script 4 prints, for each equation, how many rows have every variable. The `Effort_indicators` row is the minimum; FIML can use rows missing some indicators.
- Width model and flow version are set in script 4 (`widthModel`, `flowVersion`).
  - Choosing the width model needs the aperture-corrected measured flow (the corrected mouth-front version in script 3), still to come.
