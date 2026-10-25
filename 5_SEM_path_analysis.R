
#############################################################################################
# Path analysis (SEM) of the suction-feeding model                                          #
#############################################################################################
  # FIRST DRAFT. Equations follow SEM_Equations_Dynamic_KT / SEM_Equations_Static_KT (ktVersion below)
  # Input: the .rds from 4_SEM_dataset.R (one row per strike, or per fish with the _byFish file)
  #
  # Steps, run separately for each time basis (values at tmax and at tmaxVel):
  #   1. Effort (latent variable): lavaan measurement model only, with FIML (missing = "ml"),
  #      then Bartlett factor scores (Effort_score), used as an observed variable from here on
  #   2. Path equations: one regression per equation (lmer with a random intercept for fish when there are several strikes per fish,
  #      lm when there is one row per fish). Each equation uses the rows that have all of its variables
  #   3. Identities (*** in the docs): not fitted, the equation itself is the path (coefficients preset, no error term)
  #      Each one is checked against the data, so a mismatch shows an upstream problem
  #   4. Test of the whole model: d-separation claims and Fisher's C (Shipley's test, as in piecewiseSEM)
  #   5. Path effects: each equation linearized at reference values (first derivative, see the note at the top of the docs),
  #      then direct effects multiplied along every path to get indirect and total effects on the two outcomes
  #   6. Bootstrap over fish: steps 1, 2 and 5 repeated on resampled fish for confidence intervals
  #      This includes the factor scores, so their estimation error is in the intervals
  #
  # Why piecewiseSEM::psem isn't used for the fit:
  #   - Its basis set treats every term of a formula as a variable, so the nonlinear physics terms (I(...)) and the
  #     identities (*** equations) give claims that are trivially true or can't be tested
  #   - The d-separation test here is the same test (each claim = the response's own model + the other variable,
  #     conditioned on both variables' parents), but it knows which variables are identities
  #   - The mixed models are the same ones psem would fit (lme4)
  #
  # No size correction anywhere: the outcomes are physical. Size corrections are only for the genetic screening
  #
  # Needs: lavaan, lme4, lmerTest (p-values for the mixed models)

rm(list=ls())

library(lavaan)
library(lme4)
library(lmerTest)


##############
#  Settings  #
##############

  pathData <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data'   ### CHANGE

  datasetFile <- file.path(pathData, 'SEM_dataset.rds')   ### CHANGE  'SEM_dataset.rds' (per strike) or 'SEM_dataset_byFish.rds' (per fish)
  outFolder   <- file.path(pathData, 'SEM_results')       ### CHANGE  folder for the results (made if needed)

  # Script 2, for the four-bar functions and sign settings (so the four-bar identities match script 2 exactly)
  script2 <- '2_ACxTRC_Kinematic_Measurements_Sept2026.R'   ### CHANGE  path to script 2

  # Time bases to run (the dataset has both)
  timeBases <- c("tmax", "tmaxVel")

  # Kinematic transmission: "dynamic" (Oral_KT_tmax, Opercular_KT_tmax) or "static" (Oral_4bar_KT, Opercular_4bar_KT at 30 and 10 degrees)
  ktVersion <- "dynamic"   ### CHANGE

  # Random effect for the per-strike models (ignored with one row per fish)
  randomTerm <- "(1 | Fish_ID)"   ### CHECK  e.g. "(1 | Family/Fish_ID)" to add family

  # Density of water (kg/m^3), for the pressure equation
  rho <- 1000

  # Effort indicators are put on the same scale (mean 0, SD 1) before the measurement model. T is in s, the others aren't
  standardizeIndicators <- TRUE

  # Reference values for the linearization (docs: "90% of max to avoid edge fitting problems")
    # "fracMax": refFrac x the value with the largest magnitude (keeps the sign, so pressure uses 90% of the most negative)
    # "median": the median
    # Identities are then recalculated from their parents' reference values, so the reference point is consistent
  refType <- "fracMax"   ### CHECK
  refFrac <- 0.9

  # Bootstrap
  nBoot  <- 200    ### CHANGE  1000+ for the final run (200 is for checking the script)
  ciLevel <- 0.95
  set.seed(1)

  # The two performance outcomes. Effects on these are reported
  outcomes <- c("Max_sub_pressure", "Total_ingested_volume")


######################################
#  Equations (from the SEM documents) #
######################################

  ### CHECK all equations against the current SEM_Equations documents
  # R formula notation:
    # a * b      main effects of a and b plus their interaction (the "#" equations: additive and multiplicative effects)
    # I(a * b)   one term: the product itself (the physics-based "*" equations, where the product is the physical quantity)
    # I(a / b)   a ratio as one term. In a formula a plain "/" means nesting, not division, so ratios always need I()
  # {OralKT} and {OpercularKT} are replaced by the KT columns for ktVersion

  # Path equations (fitted). Name = response
  pathEqs <- c(

    # Performance 1 (*): Muller et al. eq. 39, expanded. Each term in Pa (velocities mm/s -> m/s, R mm -> m)
      # If the physics holds exactly every coefficient is 1
    Max_sub_pressure = paste("I(-0.5 * rho * (U_flow_ff_tmax / 1000)^2) + I(0.5 * rho * (U_ram_tmax / 1000)^2)",
                             "+ I(rho * R_tmax * U_flow_ff_acc_tmax / 1e6) + I(rho * R_tmax * U_ram_acc_tmax / 1e6)",
                             "+ I(rho * U_flow_ff_tmax * R_vel_tmax / 1e6) + I(rho * U_ram_tmax * R_vel_tmax / 1e6)"),

    # Performance 2 (*): T<A><U_flow_ff> + T Cov(A, U_flow_ff). Coefficients are 1 if Total_ingested_volume is in mm^3
    Total_ingested_volume = "I(T * A_mean * U_flow_ff_mean) + I(T * Cov_A_Uflowff)",   ### CHECK units of Total_ingested_volume

    # Gape (#)
    R_tmax     = "Protrusion_tmax + Mandible_depression_tmax",
    R_vel_tmax = "Protrusion_vel_tmax + Mandible_dep_vel_tmax",
    R_t_mean   = "Protrusion_t_mean + Mandible_depression_t_mean",

    # Protrusion (#)
    Protrusion_tmax     = "Neurocranium_tip_disp_tmax + Maxilla_tip_disp_tmax",
    Protrusion_vel_tmax = "Neurocranium_linear_vel_tmax + Maxilla_linear_vel_tmax",
    Protrusion_t_mean   = "Relative_ttpg * Protrusion_ROM",
    Protrusion_ROM      = "Length_premaxilla_asc_arm",

    # Mean mandible depression (#)
    Mandible_depression_t_mean = "Relative_ttpg * Mandible_ROM + ttpg * Mandible_dep_vel_tmax + Hyoid_depression_t_mean * I(Protractor_hyoideus_CSA / Mandible_length)",

    # Measured flow (#)
    U_flow_ef_meas_tmax     = "Relative_tmax * Kinetic_Synchronization + I(V_buccal_vel_tmax / A_tmax)",
    U_flow_ef_meas_acc_tmax = "Relative_tmax * Kinetic_Synchronization + I(V_buccal_acc_tmax / A_tmax) + I(V_buccal_vel_tmax * A_vel_tmax / A_tmax^2)",
    U_flow_ef_meas_mean     = "Kinetic_Synchronization + U_flow_ff_predicted_mean_bigGape",

    # Mandible rotation (#)
    Mandible_angular_vel_tmax = "Neurocranium_angular_vel_tmax * {OpercularKT} * IOP_ligament_CSA + Hyoid_vel_tmax * I(Protractor_hyoideus_CSA / Mandible_length)",

    # Hyoid (#)
    Hyoid_vel_tmax          = "Hypaxial_CSA * Effort_score * Sternohyoideus_CSA + Urohyal_length",
    Hyoid_depression_tmax   = "Relative_tmax * Hyoid_ROM",
    Hyoid_depression_t_mean = "Relative_ttpg * Hyoid_ROM",

    # Neurocranium (#). Docs: Epaxial_CSA * Effort / Neurocranium_MA, written as a product with 1/MA
    Neurocranium_linear_vel_tmax = "Epaxial_CSA * Effort_score * I(1 / Neurocranium_MA)",   ### CHECK
    Neurocranium_rotation_tmax   = "Relative_tmax * Neurocranium_ROM",

    # Maxilla (*): output angular velocity = input angular velocity x KT
    Maxilla_angular_vel_tmax = "I(Mandible_angular_vel_tmax * {OralKT})",

    # Effort (#), on the factor scores from step 1
    Effort_score = "Prey_num + Satiation + Fasting_time")

  # Check equations: how well the four-bar models predict the measured angles
    # Fitted and reported, but not part of the path model (nothing downstream uses them)
    # Docs write the left side of the mandible check as Mandible_ang_predicted(tmax). In the data that column IS the four-bar
    # prediction (change in theta4 from rest), so the measured Mandible_ang_tmax is the response here
  checkEqs <- c(
    Maxilla_ang_tmax  = "Maxilla_ang_predicted_tmax",
    Mandible_ang_tmax = "Mandible_ang_predicted_tmax * IOP_ligament_CSA + Hyoid_depression_tmax * I(Protractor_hyoideus_CSA / Mandible_length)")   ### CHECK

  # Identities (***): response = R expression of other columns. Used for the path effects and checked against the data
    # Only the ones connected to the path equations are used (the static and dynamic KT are both listed)
    # The four-bar functions and the sign settings come from script 2
  oralLinks <- "Mandible_input_length, Maxilla_length, Maxilla_palatine_length, Oral_fixed_length"
  operLinks <- "Operculum_length, IOP_link_length, RA_process_length, Opercular_fixed_length"
  identityEqs <- c(
    U_flow_ff_tmax           = "U_flow_ef_meas_tmax - U_ram_tmax",
    U_flow_ff_acc_tmax       = "U_flow_ef_meas_acc_tmax - U_ram_acc_tmax",
    U_flow_ff_mean           = "U_flow_ef_meas_mean - U_ram_mean",
    Mandible_dep_vel_tmax    = "Mandible_length * Mandible_angular_vel_tmax * cos(Mandible_ang_tmax)",
    Mandible_depression_tmax = "Mandible_length * sin(Mandible_ang_tmax)",
    A_mean                   = "(pi * W_jaw / 2) * R_t_mean",
    Cov_A_Uflowff            = "(pi * W_jaw / 2) * Cov_R_Uflowff",
    Cov_R_Uflowff            = "Cov_R_Uflowefmeas - Cov_R_Uram",
    Neurocranium_angular_vel_tmax = "Neurocranium_linear_vel_tmax / Neurocranium_length",
    Neurocranium_tip_disp_tmax    = "Neurocranium_rotation_tmax * Neurocranium_length",
    Neurocranium_MA          = "Supraoccipital_crest_depth / Neurocranium_length",
    Relative_tmax            = "tmax_from_Tstart / T",
    Relative_ttpg            = "ttpg / T",
    Kinetic_Synchronization  = "apply(cbind(Time_hyoid, ttpg, Time_cranial), 1, sd) / T",
    # Oral four-bar (dynamic)
    Oral_theta2_tmax  = "Oral_theta2_rest + oralOpenSign * Mandible_ang_tmax",
    Oral_KT_tmax      = paste0("maxillaSign * fourBarKTinst(Oral_theta2_tmax, ", oralLinks, ", output = 'theta3')"),
    Oral_theta3_rest  = paste0("fourBarTheta3(Oral_theta2_rest, ", oralLinks, ")"),
    Maxilla_ang_predicted_tmax = paste0("maxillaSign * wrapAngle(fourBarTheta3(Oral_theta2_tmax, ", oralLinks, ") - Oral_theta3_rest)"),
    # Opercular four-bar (dynamic)
    Opercular_theta2_tmax = "Opercular_theta2_rest + opercularOpenSign * Neurocranium_rotation_tmax",
    Opercular_KT_tmax     = paste0("fourBarKTinst(Opercular_theta2_tmax, ", operLinks, ", output = 'theta4')"),
    Opercular_theta4_rest = paste0("fourBarTheta4(Opercular_theta2_rest, ", operLinks, ")"),
    Mandible_ang_predicted_tmax = paste0("wrapAngle(fourBarTheta4(Opercular_theta2_tmax, ", operLinks, ") - Opercular_theta4_rest)"),
    # Static KT (30 degrees oral, 10 degrees opercular)
    Oral_4bar_KT      = paste0("maxillaSign * wrapAngle(fourBarTheta3(Oral_theta2_rest + oralOpenSign * oralStaticAngle, ", oralLinks,
                               ") - Oral_theta3_rest) / (oralOpenSign * oralStaticAngle)"),
    Opercular_4bar_KT = paste0("wrapAngle(fourBarTheta4(Opercular_theta2_rest + opercularOpenSign * opercularStaticAngle, ", operLinks,
                               ") - Opercular_theta4_rest) / (opercularOpenSign * opercularStaticAngle)"))

  # Names in the equations that are constants, not columns
  constantNames <- c("rho", "pi", "oralOpenSign", "opercularOpenSign", "maxillaSign", "oralStaticAngle", "opercularStaticAngle")

  # Effort measurement model (step 1). Indicators only: the covariates are in the Effort_score path equation
  effortIndicators <- c("Neurocranium_Adj_vel_mean", "Hyoid_Adj_vel_mean", "T")

  # Pairs never tested in the d-separation (the factor scores are made from the indicators, so they're dependent by construction)
  dsepSkipPairs <- lapply(effortIndicators, function(v) c("Effort_score", v))


##### Four-bar functions and sign settings from script 2
  src2 <- readLines(script2)
  fbStart <- grep("wrapAngle <- function", src2, fixed = TRUE)
  fbEnd   <- grep("# Correlation for the checks", src2, fixed = TRUE) - 1
  if(length(fbStart) != 1 || length(fbEnd) != 1) {
    stop("Can't find the four-bar functions in script 2 (wrapAngle ... before '# Correlation for the checks')")
  }
  eval(parse(text = src2[fbStart:fbEnd]))
  signLines <- grep("^\\s*(oralOpenSign|opercularOpenSign|maxillaSign|oralStaticAngle|opercularStaticAngle)\\s*<-", src2, value = TRUE)
  if(length(signLines) != 5) {
    stop("Expected 5 four-bar sign / angle settings in script 2, found ", length(signLines))
  }
  eval(parse(text = signLines))


##### KT version
  ktCols <- list(dynamic = c(OralKT = "Oral_KT_tmax", OpercularKT = "Opercular_KT_tmax"),
                 static  = c(OralKT = "Oral_4bar_KT", OpercularKT = "Opercular_4bar_KT"))[[ktVersion]]
  for(nm in names(ktCols)) {
    pathEqs <- gsub(paste0("{", nm, "}"), ktCols[[nm]], pathEqs, fixed = TRUE)
  }


###############
#  Functions  #
###############

  # Variables (columns) in an equation's right side
  varsOf <- function(rhs) {
    setdiff(all.vars(parse(text = rhs)), constantNames)
  }

  # Fit one equation: lmer with randomTerm if there are several rows per fish, otherwise lm
    # Uses the rows that have every variable in the equation (and Fish_ID)
  fitEq <- function(response, rhs, d, extraTerms = NULL) {
    rhsNow  <- paste(c(rhs, extraTerms), collapse = " + ")
    varsNow <- c(response, varsOf(rhsNow), "Fish_ID")
    dNow <- d[complete.cases(d[, varsNow, drop = FALSE]), , drop = FALSE]
    if(multiPerFish) {
      lmer(as.formula(paste(response, "~", rhsNow, "+", randomTerm)), data = dNow)
    } else {
      lm(as.formula(paste(response, "~", rhsNow)), data = dNow)
    }
  }

  # Prediction from the fixed effects only (the "average fish")
  predFixed <- function(m, newdata) {
    if(inherits(m, "merMod")) predict(m, newdata = newdata, re.form = NA) else predict(m, newdata = newdata)
  }

  # Coefficient table: estimate, SE, p for each fixed-effect term
  coefTable <- function(m, eqName) {
    ct <- summary(m)$coefficients
    data.frame(equation = eqName, term = rownames(ct), estimate = ct[, "Estimate"], se = ct[, "Std. Error"],
               p = if(ncol(ct) >= 4) ct[, ncol(ct)] else NA_real_, n = nobs(m), row.names = NULL)
  }

  # R^2. Mixed models: marginal (fixed effects only) and conditional (fixed + random), Nakagawa & Schielzeth
  rSquared <- function(m) {
    if(inherits(m, "merMod")) {
      varF  <- var(predict(m, re.form = NA))
      varRE <- sum(sapply(VarCorr(m), function(v) sum(diag(v))))
      varE  <- sigma(m)^2
      return(c(R2_marginal = varF / (varF + varRE + varE), R2_conditional = (varF + varRE) / (varF + varRE + varE)))
    }
    c(R2_marginal = summary(m)$r.squared, R2_conditional = NA_real_)
  }

  # Order the variables so every variable comes after its parents. Stops if there is a cycle
  topoOrder <- function(parents) {
    done <- character(0)
    while(length(done) < length(parents)) {
      ready <- names(parents)[!(names(parents) %in% done) & sapply(parents, function(p) all(p %in% done))]
      if(length(ready) == 0) {
        stop("The path model has a cycle among: ", paste(setdiff(names(parents), done), collapse = ", "))
      }
      done <- c(done, ready)
    }
    done
  }

  # Reference value for the linearization
  refValue <- function(x) {
    if(is.factor(x)) {
      return(factor(levels(x)[ceiling(nlevels(x) / 2)], levels = levels(x), ordered = is.ordered(x)))
    }
    x <- x[is.finite(x)]
    if(length(x) == 0) return(NA_real_)
    if(refType == "median") median(x) else refFrac * x[which.max(abs(x))]
  }

  # Derivative of fun (a function of a one-row data frame) with respect to column x, at refRow
    # Numeric columns: central difference. Ordered factor (Prey_num): mean change per step, (top level - bottom level) / (levels - 1)
  derivAt <- function(fun, refRow, x, scaleNow) {
    if(is.factor(refRow[[x]])) {
      levs <- levels(refRow[[x]])
      lo <- refRow; lo[[x]] <- factor(levs[1], levels = levs, ordered = TRUE)
      hi <- refRow; hi[[x]] <- factor(levs[length(levs)], levels = levs, ordered = TRUE)
      return((fun(hi) - fun(lo)) / (length(levs) - 1))
    }
    h <- 1e-4 * (if(is.finite(scaleNow) && scaleNow > 0) scaleNow else 1)
    up <- refRow; up[[x]] <- up[[x]] + h
    dn <- refRow; dn[[x]] <- dn[[x]] - h
    (fun(up) - fun(dn)) / (2 * h)
  }

  # SD used to standardize effects (ordered factor: SD of the level number, to match the per-step derivative)
  sdOf <- function(x) {
    if(is.factor(x)) sd(as.numeric(x), na.rm = TRUE) else sd(x, na.rm = TRUE)
  }


  ##### Step 1: Effort factor scores
    # Measurement model only, FIML so a strike missing one indicator is kept
    # Bartlett scores: unbiased for the factor, so slopes on Effort_score later aren't attenuated (regression scores would be)
    # Sign set so Effort increases with Neurocranium_Adj_vel_mean
  fitEffort <- function(d) {
    ind <- d[, effortIndicators]
    if(standardizeIndicators) {
      ind[] <- lapply(ind, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    }
    keep <- which(rowSums(!is.na(ind)) > 0)

    fit <- cfa(paste("Effort =~", paste(effortIndicators, collapse = " + ")), data = ind[keep, ], missing = "ml", std.lv = TRUE)
    scores <- rep(NA_real_, nrow(d))
    scores[keep[lavInspect(fit, "case.idx")]] <- lavPredict(fit, method = "Bartlett")[, "Effort"]

    flipNow <- if(cor(scores, d[[effortIndicators[1]]], use = "pairwise.complete.obs") < 0) -1 else 1
    loadings <- standardizedSolution(fit)
    loadings <- loadings[loadings$op == "=~", c("rhs", "est.std", "se", "pvalue")]
    loadings$est.std <- flipNow * loadings$est.std
    list(score = flipNow * scores, loadings = loadings, fit = fit)
  }


  ##### Steps 2 and 5 (and 4 if doDsep): fit everything on one dataset and calculate the path effects
  runModel <- function(d, doDsep = FALSE) {

    effortNow <- fitEffort(d)
    d$Effort_score <- effortNow$score

    # Step 2: path and check equations
    fits <- list()
    for(y in names(pathEqs)) {
      fits[[y]] <- tryCatch(fitEq(y, pathEqs[[y]], d), error = function(e) stop("Equation ", y, " failed: ", conditionMessage(e)))
    }
    checkFits <- lapply(setNames(names(checkEqs), names(checkEqs)), function(y) fitEq(y, checkEqs[[y]], d))

    # Reference row: every column at its reference value, then identities recalculated from their parents in path order
    refRow <- as.data.frame(lapply(d[, setdiff(names(d), c("vidName", "Fish_ID", "Family", "Animal.ID", "timeBasis", "Hyoid_depth_source")), drop = FALSE],
                                   function(x) if(is.numeric(x) || is.factor(x)) refValue(x) else NA))
    for(v in intersect(nodeOrder, names(activeIdentities))) {
      refRow[[v]] <- eval(parse(text = activeIdentities[[v]])[[1]], refRow, globalenv())
    }

    # Step 5: direct effects B[y, x] = dy/dx at the reference row, for every arrow x -> y
    B <- matrix(0, length(nodeOrder), length(nodeOrder), dimnames = list(nodeOrder, nodeOrder))
    for(y in names(nodeParents)) {
      if(length(nodeParents[[y]]) == 0) next
      funNow <- if(y %in% names(pathEqs)) {
        local({ m <- fits[[y]]; function(nd) predFixed(m, nd) })
      } else {
        local({ ex <- parse(text = activeIdentities[[y]])[[1]]; function(nd) eval(ex, nd, globalenv()) })
      }
      for(x in nodeParents[[y]]) {
        B[y, x] <- derivAt(funNow, refRow, x, sdOf(d[[x]]))
      }
    }
    badB <- which(!is.finite(B), arr.ind = TRUE)
    if(nrow(badB) > 0) {
      warning("Direct effects that couldn't be calculated at the reference values (set to 0): ",
              paste(paste0(colnames(B)[badB[, 2]], " -> ", rownames(B)[badB[, 1]]), collapse = ", "))
      B[!is.finite(B)] <- 0
    }

    # Total effects: (I - B)^-1 - I, element [y, x] = sum over all paths from x to y of the product of the direct effects
    totalE <- solve(diag(nrow(B)) - B) - diag(nrow(B))
    dimnames(totalE) <- dimnames(B)

    # Standardized (x 1 SD of the cause / 1 SD of the outcome)
    sds <- sapply(nodeOrder, function(v) sdOf(d[[v]]))
    totalStd  <- totalE * outer(1 / sds, sds)
    directStd <- B * outer(1 / sds, sds)

    out <- list(effort = effortNow, fits = fits, checkFits = checkFits, refRow = refRow,
                B = B, totalE = totalE, directStd = directStd, totalStd = totalStd,
                coefs = do.call(rbind, lapply(names(fits), function(y) coefTable(fits[[y]], y))))
    if(doDsep) {
      out$dsep <- dSepTest(d, fits)
    }
    out
  }


  ##### Step 4: d-separation test
    # Every pair of variables with no arrow between them is a claim: independent, given both variables' parents
    # The response of a claim is the one later in the path order. It must have a fitted equation:
      # claims where both variables are exogenous (no equation) aren't tested
      # identities are never the response (given their parents they have no variation left)
    # Claim test: the response's own equation + the parents of the other variable + the other variable,
      # compared with the same model without the other variable (likelihood ratio for mixed models, F test for lm)
    # Skipped as untestable: the other variable is fully determined by the variables already in the model
      # (an identity of them, or linearly predicted with R^2 > 0.9999)
  dSepTest <- function(d, fits) {

    # Is v fixed by the variables in condSet (through identities)?
    determinedBy <- function(v, condSet) {
      if(v %in% condSet) return(TRUE)
      if(v %in% names(activeIdentities)) return(all(sapply(nodeParents[[v]], determinedBy, condSet = condSet)))
      FALSE
    }

    claims <- data.frame()
    for(a in seq_along(nodeOrder)) {
      for(b in seq_along(nodeOrder)) {
        if(b <= a) next
        x <- nodeOrder[a]; y <- nodeOrder[b]   # y is later in the path order
        if(x %in% nodeParents[[y]] || y %in% nodeParents[[x]]) next
        if(any(sapply(dsepSkipPairs, function(p) setequal(p, c(x, y))))) next
        # Response must be a fitted equation (prefer the later one)
        if(!(y %in% names(pathEqs))) {
          if(x %in% names(pathEqs)) { tmp <- x; x <- y; y <- tmp } else next
        }

        condVars <- union(varsOf(pathEqs[[y]]), nodeParents[[x]])
        rowNow <- data.frame(response = y, other = x, conditioning = paste(setdiff(condVars, x), collapse = ", "),
                             n = NA_integer_, p = NA_real_, note = "", stringsAsFactors = FALSE)

        if(determinedBy(x, setdiff(condVars, x))) {
          rowNow$note <- "skipped: other variable is an identity of the conditioning variables"
          claims <- rbind(claims, rowNow)
          next
        }

        extraNow <- setdiff(nodeParents[[x]], c(varsOf(pathEqs[[y]]), x))
        varsAll  <- unique(c(y, varsOf(pathEqs[[y]]), extraNow, x, "Fish_ID"))
        dClaim <- d[complete.cases(d[, varsAll, drop = FALSE]), , drop = FALSE]
        rowNow$n <- nrow(dClaim)
        if(nrow(dClaim) < 10) {
          rowNow$note <- "skipped: fewer than 10 rows"
          claims <- rbind(claims, rowNow)
          next
        }

        # Linear determination check (e.g. U_flow_ef_meas = U_flow_ff + U_ram when both are conditioned on)
        if(is.numeric(dClaim[[x]])) {
          mmNow <- model.matrix(as.formula(paste("~", paste(c(pathEqs[[y]], extraNow), collapse = " + "))), dClaim)
          r2Now <- summary(lm(dClaim[[x]] ~ mmNow))$r.squared
          if(is.finite(r2Now) && r2Now > 0.9999) {
            rowNow$note <- "skipped: other variable is linearly determined by the conditioning variables"
            claims <- rbind(claims, rowNow)
            next
          }
        }

        res <- tryCatch({
          m0 <- fitEq(y, pathEqs[[y]], dClaim, extraTerms = extraNow)
          m1 <- fitEq(y, pathEqs[[y]], dClaim, extraTerms = c(extraNow, x))
          a01 <- suppressMessages(anova(m0, m1))
          a01[2, ncol(a01)]
        }, error = function(e) { rowNow$note <<- paste("failed:", conditionMessage(e)); NA_real_ })
        rowNow$p <- res
        claims <- rbind(claims, rowNow)
      }
    }

    tested <- claims[!is.na(claims$p), ]
    fisherC <- -2 * sum(log(pmax(tested$p, 1e-300)))
    dfC <- 2 * nrow(tested)
    list(claims = claims[order(claims$p), ], fisherC = fisherC, df = dfC,
         p = if(dfC > 0) pchisq(fisherC, dfC, lower.tail = FALSE) else NA_real_)
  }


######################
#  Read the data     #
######################

  datasets <- readRDS(datasetFile)
  dir.create(outFolder, showWarnings = FALSE)

  # Several rows per fish (per-strike data): mixed models. One row per fish: lm
  multiPerFish <- anyDuplicated(datasets[[1]]$Fish_ID) > 0
  print(paste0("Unit of analysis: ", if(multiPerFish) paste0("strike (lmer with ", randomTerm, ")") else "fish (lm)"))
  print(paste0("KT: ", ktVersion, " (", paste(ktCols, collapse = ", "), ")"))

  if(!is.ordered(datasets[[1]]$Prey_num)) {
    stop("Prey_num isn't an ordered factor: read the .rds from script 4, not a .csv")
  }


#######################################
#  Path structure (same for all runs) #
#######################################

  # Parents of every variable in the path model
    # Identities are kept if their response is used by a path equation or by another kept identity
  dataNames <- c(names(datasets[[1]]), "Effort_score")
  usedVars <- unique(c(names(pathEqs), unlist(lapply(pathEqs, varsOf))))
  activeIdentities <- character(0)
  repeat {
    newIds <- setdiff(intersect(names(identityEqs), usedVars), names(activeIdentities))
    if(length(newIds) == 0) break
    activeIdentities <- c(activeIdentities, identityEqs[newIds])
    usedVars <- unique(c(usedVars, unlist(lapply(identityEqs[newIds], varsOf))))
  }
  if(any(names(activeIdentities) %in% names(pathEqs))) {
    stop("Variables with both a path equation and an identity: ", paste(intersect(names(activeIdentities), names(pathEqs)), collapse = ", "))
  }

  nodeParents <- setNames(vector("list", length(usedVars)), usedVars)
  for(v in usedVars) {
    nodeParents[[v]] <- if(v %in% names(pathEqs)) varsOf(pathEqs[[v]]) else if(v %in% names(activeIdentities)) varsOf(activeIdentities[[v]]) else character(0)
  }
  nodeOrder <- topoOrder(nodeParents)

  missingVars <- setdiff(c(usedVars, unlist(lapply(checkEqs, varsOf)), names(checkEqs), effortIndicators), dataNames)
  if(length(missingVars) > 0) {
    stop("Columns not in the dataset (rerun script 4?): ", paste(missingVars, collapse = ", "))
  }

  exogenous <- names(nodeParents)[sapply(nodeParents, length) == 0]
  print(paste0(length(pathEqs), " path equations, ", length(activeIdentities), " identities, ", length(exogenous), " exogenous variables"))


#######################################
#  Run                                #
#######################################

  results <- list()

  for(tb in timeBases) {

    print("")
    print(paste0("################ Time basis: ", tb, " ################"))
    d <- datasets[[tb]]


    ##### Identity check: each identity recalculated from the data, compared with the column
      # Large differences mean the upstream scripts don't calculate the variable the way the equation says
      # Every identity whose columns are in the data, including ones not used in the path model (e.g. the other KT version)
    idsToCheck <- identityEqs[sapply(names(identityEqs), function(v) all(c(v, varsOf(identityEqs[[v]])) %in% names(d)))]
    idCheck <- do.call(rbind, lapply(names(idsToCheck), function(v) {
      calc <- eval(parse(text = idsToCheck[[v]])[[1]], d, globalenv())
      ok <- is.finite(calc) & is.finite(d[[v]])
      relDiff <- abs(calc[ok] - d[[v]][ok]) / pmax(abs(d[[v]][ok]), 1e-12)
      data.frame(identity = v, in_path_model = v %in% names(activeIdentities), n = sum(ok), median_rel_diff = if(any(ok)) median(relDiff) else NA,
                 frac_over_1pct = if(any(ok)) mean(relDiff > 0.01) else NA)
    }))
    print("Identity check (relative difference between the equation and the column; should be ~0):")
    print(idCheck, row.names = FALSE, digits = 3)


    ##### Main fit (steps 1, 2, 4, 5)
    main <- runModel(d, doDsep = TRUE)

    print("Effort measurement model (standardized loadings):")
    print(main$effort$loadings, row.names = FALSE, digits = 3)
    if(!isTRUE(lavInspect(main$effort$fit, "post.check"))) {
      print("WARNING: Effort model has a negative variance or a non-positive-definite matrix (Heywood case). Check the loadings")
    }

    eqSummary <- do.call(rbind, lapply(names(main$fits), function(y) {
      data.frame(equation = y, n = nobs(main$fits[[y]]), t(rSquared(main$fits[[y]])))
    }))
    print("Path equations (rows used and R^2):")
    print(eqSummary, row.names = FALSE, digits = 3)

    print("Check equations (four-bar predicted vs. measured angles):")
    for(y in names(main$checkFits)) {
      print(paste0(y, ": n = ", nobs(main$checkFits[[y]]), ", R2 = ", round(rSquared(main$checkFits[[y]])[1], 3)))
      print(coefTable(main$checkFits[[y]], y)[, c("term", "estimate", "se", "p")], row.names = FALSE, digits = 3)
    }

    claimsNow <- main$dsep$claims
    print(paste0("d-separation: ", sum(!is.na(claimsNow$p)), " claims tested, ", sum(is.na(claimsNow$p)), " skipped. ",
                 "Fisher's C = ", round(main$dsep$fisherC, 2), ", df = ", main$dsep$df, ", p = ", signif(main$dsep$p, 3),
                 " (p > 0.05: the data are consistent with the model)"))
    print("Smallest p-values (candidate missing paths):")
    print(head(claimsNow[!is.na(claimsNow$p), c("response", "other", "n", "p")], 10), row.names = FALSE, digits = 3)


    ##### Step 6: bootstrap over fish
      # Fish are resampled with all their strikes. A fish drawn twice gets two IDs, so the random effect treats them as two fish
    fishIDs <- unique(d$Fish_ID)
    bootTotal <- array(NA_real_, dim = c(nBoot, length(outcomes), length(nodeOrder)), dimnames = list(NULL, outcomes, nodeOrder))
    bootTotalStd <- bootTotal
    bootCoefs <- list()
    nFail <- 0
    for(k in seq_len(nBoot)) {
      pick <- sample(fishIDs, length(fishIDs), replace = TRUE)
      dBoot <- do.call(rbind, lapply(seq_along(pick), function(j) {
        rowsNow <- d[d$Fish_ID == pick[j], , drop = FALSE]
        rowsNow$Fish_ID <- paste0(pick[j], "_b", j)
        rowsNow
      }))
      runNow <- tryCatch(suppressWarnings(suppressMessages(runModel(dBoot))), error = function(e) NULL)
      if(is.null(runNow)) {
        nFail <- nFail + 1
        next
      }
      bootTotal[k, , ]    <- runNow$totalE[outcomes, ]
      bootTotalStd[k, , ] <- runNow$totalStd[outcomes, ]
      bootCoefs[[k]] <- setNames(runNow$coefs$estimate, paste(runNow$coefs$equation, runNow$coefs$term, sep = " : "))
      if(k %% 50 == 0) print(paste0("Bootstrap ", k, " / ", nBoot))
    }
    if(nFail > 0) {
      print(paste0(nFail, " of ", nBoot, " bootstrap samples failed to fit and were left out"))
    }
    lowQ <- (1 - ciLevel) / 2
    upQ  <- 1 - lowQ


    ##### Coefficients with bootstrap intervals
    coefsNow <- main$coefs
    bootCoefMat <- do.call(rbind, lapply(bootCoefs, function(b) b[paste(coefsNow$equation, coefsNow$term, sep = " : ")]))
    if(!is.null(bootCoefMat)) {
      coefsNow$boot_low  <- apply(bootCoefMat, 2, quantile, lowQ, na.rm = TRUE)
      coefsNow$boot_high <- apply(bootCoefMat, 2, quantile, upQ, na.rm = TRUE)
    }


    ##### Effects on the outcomes: direct, indirect, total (raw units and standardized), with bootstrap intervals
    effects <- do.call(rbind, lapply(outcomes, function(o) {
      causes <- setdiff(nodeOrder[main$totalE[o, ] != 0], o)
      if(length(causes) == 0) return(NULL)
      data.frame(outcome = o, cause = causes, exogenous = causes %in% exogenous,
                 direct = main$B[o, causes], indirect = main$totalE[o, causes] - main$B[o, causes], total = main$totalE[o, causes],
                 total_low = apply(bootTotal[, o, causes, drop = FALSE], 3, quantile, lowQ, na.rm = TRUE),
                 total_high = apply(bootTotal[, o, causes, drop = FALSE], 3, quantile, upQ, na.rm = TRUE),
                 total_std = main$totalStd[o, causes],
                 total_std_low = apply(bootTotalStd[, o, causes, drop = FALSE], 3, quantile, lowQ, na.rm = TRUE),
                 total_std_high = apply(bootTotalStd[, o, causes, drop = FALSE], 3, quantile, upQ, na.rm = TRUE),
                 row.names = NULL)
    }))
    effects <- effects[order(effects$outcome, -abs(effects$total_std)), ]
    for(o in outcomes) {
      print(paste0("Largest standardized total effects on ", o, ":"))
      print(head(effects[effects$outcome == o, c("cause", "direct", "total", "total_std", "total_std_low", "total_std_high")], 15),
            row.names = FALSE, digits = 3)
    }


    ##### Save
    write.csv(coefsNow,  file.path(outFolder, paste0("coefficients_", tb, ".csv")), row.names = FALSE)
    write.csv(effects,   file.path(outFolder, paste0("effects_", tb, ".csv")), row.names = FALSE)
    write.csv(claimsNow, file.path(outFolder, paste0("dsep_claims_", tb, ".csv")), row.names = FALSE)
    write.csv(eqSummary, file.path(outFolder, paste0("equations_", tb, ".csv")), row.names = FALSE)
    write.csv(idCheck,   file.path(outFolder, paste0("identity_check_", tb, ".csv")), row.names = FALSE)

    results[[tb]] <- list(main = main, coefs = coefsNow, effects = effects, idCheck = idCheck,
                          bootTotal = bootTotal, bootTotalStd = bootTotalStd, nBootFailed = nFail)
  }

  results$settings <- list(datasetFile = datasetFile, ktVersion = ktVersion, randomTerm = randomTerm, refType = refType,
                           refFrac = refFrac, nBoot = nBoot, standardizeIndicators = standardizeIndicators,
                           pathEqs = pathEqs, checkEqs = checkEqs, identities = activeIdentities)
  saveRDS(results, file.path(outFolder, paste0("SEM_results_", ktVersion, if(multiPerFish) "" else "_byFish", ".rds")))
  print(paste0("Results saved in: ", normalizePath(outFolder)))
