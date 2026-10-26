# SEM path export: nodes, arrows and coefficients

These files describe the fitted path model (structural equation model, SEM) of suction feeding in ACxTRC hybrid fish: every variable and every arrow, with its strength. They are for redrawing the path, with each line's grey level showing how strong that connection is.

`export_info.md` in the same folder names the model these files come from, with its sample size and fit.

## Where the numbers come from

The pipeline runs six scripts in order:
- **`2_ACxTRC_Kinematic_Measurements_Sept2026.R`:** kinematics from tracked landmarks, one row per strike.
- **`3_PIV_Flow_Measurements.R`:** flow and pressure from particle image velocimetry (PIV). Both are corrected to the mouth aperture with a flow-decay exponent of k = 0.61, estimated from the data.
- **`4_SEM_dataset.R`:** combines kinematics, flow, pressure and morphology into one row per strike.
- **`5_SEM_path_analysis.R`:** fits the model.
- **`6_SEM_exports.R`:** writes these files from script 5's saved results, without refitting.

How script 5 fits the model:
- **Piecewise SEM.** Each equation is fitted separately, as a mixed model with a random intercept for fish (`lme4::lmer`), or as `lm` when the response has one value per fish.
- **Identities are not fitted.** These are variables calculated exactly from others (e.g. gape area = π · jaw width / 2 · gape radius, or the four-bar angles).
- **Effort is a latent variable.** It is measured by three indicators in `lavaan`, then split into a fish-level part (`Behavior_fish`, the fish's mean score) and a strike-level part (`Effort_strike`).
- **Direct effects are slopes.** The direct effect of a parent on a child is the derivative of the child's equation with respect to the parent. It is evaluated at a reference point: every variable at 90% of its largest value, except the four-bar link lengths and rest angles, which are at their medians.
  - Many equations contain products and interactions, so this slope is a local, linearised effect at that point. It is not a single regression coefficient.

## Files

### `path_nodes.csv`: one row per variable

| column | meaning |
|---|---|
| `variable` | variable name, as in the scripts |
| `role` | `exogenous` (no arrows in: morphology, timing inputs, covariates), `fitted` (has a regression equation), `identity` (calculated exactly from its parents), `outcome` (performance: `Max_sub_pressure` = peak sub-ambient pressure, Pa; `Total_ingested_volume` = water volume ingested, mm³), `latent` (the Effort score) |
| `group_guess` | rough category assigned from the name, for layout and colour: morphology, kinematics, timing, buccal volume, flow, ram, effort, measurement covariate, performance. Check it and edit it if needed |
| `equation` | the right-hand side of the variable's equation, in R formula syntax. `*` means main effects plus interaction, `I(...)` is a single calculated term, and `:` marks an interaction term in the coefficient table |
| `n_rows` | strikes used to fit the equation |
| `R2_marginal`, `R2_conditional` | variance explained by the fixed effects, and by fixed + fish random effects (Nakagawa & Schielzeth) |
| `mean_strikes`, `sd_strikes` | mean and SD across strikes |

### `path_edges.csv`: one row per arrow (parent → child)

| column | meaning |
|---|---|
| `from`, `to` | parent and child |
| `edge_type` | `fitted` (regression arrow), `identity` (exact calculation, no error), `measurement (loading)` (Effort score → its indicators), `decomposition` (Effort score → its fish-mean and strike parts) |
| `direct` | direct effect in raw units: change in the child per unit of the parent, at the reference point |
| `direct_std` | standardized direct effect = `direct` × SD(parent) / SD(child), both SDs across strikes. Use this for line strength. For loadings it is the standardized loading |
| `direct_std_low`, `direct_std_high` | bootstrap percentile interval over fish, if the model was bootstrapped; otherwise NA |
| `regression_terms` | the terms of the child's equation that contain the parent, with estimate and p, e.g. `R_tmax = 4.2 (p 0.003)` |
| `min_p` | the smallest p among those terms (for a loading, the loading's p) |
| `abs_std` | absolute value of `direct_std` |
| `plot_weight` | `abs_std` divided by the 95th percentile of `abs_std` over the fitted arrows, capped at 1 |
| `sign` | `+` or `-` |

### `path_equation_terms.csv`

Script 5's full coefficient table: one row per term of every fitted equation, with estimate, SE, p, n, and bootstrap interval if bootstrapped.

### `path_model.dot`

A ready-made Graphviz drawing of the whole path:
- Line darkness and width follow `plot_weight`.
- Identity arrows and identity variables are dashed, and the outcomes are bold.
- Render it with `dot -Tpdf path_model.dot -o path_model.pdf`, or paste it into an online Graphviz viewer.

## How to redraw the path

- **Line darkness:** set it from `plot_weight`, with grey = 1 − `plot_weight` (0 = black). Example in R:
  ```r
  edges$grey <- grey(1 - edges$plot_weight)
  ```
  - `plot_weight` is scaled to the 95th percentile, so a few very strong arrows don't wash out the rest.
  - To scale on the raw strength instead, use `abs_std / max(abs_std)`.
- **Line style:**
  - Draw identity arrows dashed (or leave them out of a simplified figure). Their "strength" is a calculation, not an estimated relationship.
  - The loadings are on a different scale (correlation-like), so draw them in their own style.
- **Sign:** `sign` can set a second channel, for example colour or a solid vs. dotted line.
- **Layout:** layers that follow the biology work well:

  > morphology → timing and Effort → kinematics → buccal volume → flow → performance

  `group_guess` gives a first pass at these layers, and the measurement covariates (`*_node_zR`: how far the PIV node was from the mouth) can sit off to the side.
- **Uncertainty:** with a bootstrap, an arrow whose interval (`direct_std_low` to `direct_std_high`) spans 0 can be drawn lighter or dotted.
- **Simplified figure:** keep the fitted arrows with `abs_std` above a threshold, plus the arrows into the outcomes.

## Cautions

- `direct_std` values above 1 are possible. Some equations are products (e.g. flow = volume rate / gape area), so a local slope × SD ratio can exceed 1. They still rank the arrows correctly.
- The direct effect of a parent that appears in an interaction depends on the other variable's value at the reference point. Its regression terms are listed separately in `regression_terms`.
- The model's fit test (Fisher's C, in `export_info.md`) rejects it, as is usual with ~1700 independence claims. The strengths are still the estimated effects, but there are correlations the model doesn't include.
