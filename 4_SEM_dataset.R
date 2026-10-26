
#############################################################################################
# Dataset for the path analysis (piecewiseSEM)                                              #
#############################################################################################
  # Joins the outputs of scripts 1-3 and the feeding data sheet into one table per time basis:
    # SEM_dataset_tmax.csv:    every instantaneous variable at tmax (time of peak pressure)
    # SEM_dataset_tmaxVel.csv: the same variables at tmaxVel (time of peak flow velocity)
  # Both files use the SAME column names (the "_tmax" names from SEM_Equations), so one model formula runs on either
    # The timeBasis column says which time a file is
  # Missing values stay NA. piecewiseSEM fits each equation on the rows that have that equation's variables,
    # so a video without pressure still contributes to the kinematic and morphological paths

  # Choices that aren't settled yet are settings (### CHANGE): one row per strike or per fish, the epaxial section,
    # the buccal width model, the measured-flow version, and whether videos with no real strike end are kept

  # No size correction, on purpose: the SEM models physical outcomes, so raw sizes are used. Size correction is only for the genetic mapping
  # Not done here: the nonlinear-path linearization (done in the SEM)

rm(list=ls())


##############
#  Settings  #
##############

  pathData <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data'

  kinFile   <- file.path(pathData, 'kinematicsAll_sleapOutput.csv')         # 2_ACxTRC_Kinematic_Measurements_Sept2026.R
  flowFile  <- file.path(pathData, 'flowMeasurements.csv')                  # 3_PIV_Flow_Measurements.R
  scanFile  <- file.path(pathData, 'scan_measurements.csv')                 # 1_Scan_measurements_from_landmarks.R
  sheetFile <- file.path(pathData, 'ACxTRCfeedingdata_Oct2026_FIXED.csv')   # feeding data sheet (same file as script 2)

  outPrefix <- file.path(pathData, 'SEM_dataset')   # saves <outPrefix>_tmax.csv, <outPrefix>_tmaxVel.csv and <outPrefix>.rds
    # Use the .rds in script 5: it keeps Prey_num as an ordered factor, which a .csv can't store
    # The .rds has one dataset per time basis x width model x flow version (names like "tmax|main|mouthFrontCorr"),
      # so script 5 can compare them. The CSVs are only the widthModel / flowVersion chosen below

  # One row per strike (keep Fish_ID for a random effect in each equation) or per fish
    # Per fish: every value comes from the fish's strike with the most negative Max_sub_pressure (its strongest suction)
  unitOfAnalysis <- "strike"   ### CHANGE  "strike" or "fish"

  # Per fish, for fish with no strike that has pressure
  fishNoPressure <- "largestGape"   ### CHANGE  "largestGape" (use the strike with the largest A_Max) or "drop" (leave the fish out)

  # Prey_num is an ordered categorical variable. Its levels, lowest to highest, as written in the data sheet
    # Matching ignores capitals and spaces at either end ("Few " = "few")
    # The script stops if the sheet has a value that is in neither list
  preyLevels <- c("none", "one", "few", "several", "many")   ### CHECK  spelling in the data sheet

  # Values in the prey column that flag a problem with the PIV (not a prey count)
    # These videos keep their kinematics, but their pressure and flow are set to NA, and Prey_num is NA (not recorded)
    # PIV_flagged = TRUE marks them
  pivFlagValues <- c("bad", "dontusethis", "bubble", "tube")   ### CHANGE

  # Which epaxial cross-section is Epaxial_CSA
  epaxialSection <- "crest"    ### CHANGE  "crest" (highest point of the supraoccipital crest) or "flat" (where the neurocranium flattens out)

  # Buccal width models (V_buccal and U_flow_ff_predicted) and measured-flow versions (script 3) to save in the .rds
    # A flow version is one of script 3's versions ("static", "staticRaw", "mouthFrontCorr", "mouthFront"; "" = script 3's mainVersion),
    # or two joined with "+": "instant+means"
      # instant: values at tmax / tmaxVel (flow and acceleration, earth and fish frame)
      # means:   strike means, covariances and total ingested volume
      # e.g. "static+mouthFrontCorr" (default): the aperture-corrected static node at tmax (where the correction matters for pressure),
        # and the corrected mouth-front flow for the strike means (it covers ~95% of each strike, and its size matches the flow
        # predicted from buccal volume)
        # The corrected static node only covers part of each strike (frames within maxZoverR of the mouth), so its own strike means
        # exist for few videos. The uncorrected node (staticRaw) is ~1.4R from the mouth, so its fish-frame means are mostly negative
  # Buccal width models. Only noWidth (no lateral expansion) is kept: in script 2 the bar model's expansion reached ~1.5 x W_head and
    # its retraction fit was poor (R^2 ~ 0), and in script 3 noWidth tracked the measured flow best within videos (shape and volume rate)
    # Add "main" (bar model) or "noRetract" back to compare them in script 5 (compareWidthModels)
  widthModels  <- c("noWidth")   ### CHANGE  any of "noWidth", "main", "noRetract"
  flowVersions <- c("static+mouthFrontCorr", "static+staticRaw", "static", "staticRaw", "mouthFrontCorr",
                    "staticK+mouthFrontCorrK", "mouthFrontCorrK")   ### CHANGE
    # "mouthFrontCorrK": the mouth-front flow (k = mullerExponentAlt) on its own; script 5's flow measurement check uses it
    # "staticK+mouthFrontCorrK": the same as the default, corrected with script 3's second decay exponent (mullerExponentAlt)
      # Only made when script 3 was run with mullerExponentAlt set; otherwise skipped with a message

  # Time bases. "tmax" (peak pressure) and "tmaxVel" (peak velocity): every variable at the same instant
    # "peak": every variable at its own maximum over the strike (needs script 2's _peak columns and script 3's _Max columns)
      # Variables script 5 calculates from others (identities such as U_flow_ff, Mandible_depression, the four-bar angles) are left
      # NA here and recalculated from the peaks in script 5
  timeBases <- c("tmax", "tmaxVel", "peak")   ### CHANGE

  # The combination saved to the CSVs and printed in the checks (also the default in script 5)
  widthModel  <- "noWidth"         ### CHANGE
  flowVersion <- "staticK+mouthFrontCorrK"   ### CHANGE  (k = 0.61; needs mullerExponentAlt set in script 3)

  # Videos where the gape never dropped back below the threshold (Tend_censored). Their duration-based values are already NA
  dropCensored <- FALSE        ### CHANGE  TRUE drops these videos entirely

  # Buccal volume shrinking at peak velocity (V_buccal_vel_tmaxVel < 0) isn't trusted: that video's buccal-volume values
    # (V_buccal_vel, V_buccal_acc, V_buccal_Max, U_flow_ff_predicted_mean_bigGape) are NA in every dataset, tmax and tmaxVel
    # Checked separately for each width model. buccal_excluded = TRUE marks them
  dropNegBuccalVel <- TRUE     ### CHANGE

  # Max_sub_pressure: "corrected" = the PIV peak corrected to the aperture (script 3's Max_sub_pressure_corr; NA when the pressure node
    # was behind the mouth or beyond maxZoverR), or "raw" = the PIV peak at its own node. The raw value is always kept as Max_sub_pressure_raw
  pressureVersion <- "corrected"   ### CHANGE  "corrected" or "raw"

  # Effort inputs from the feeding data sheet: name in this dataset = column in the data sheet
    # Satiation isn't in the data sheet: it's the strike number for the fish (XX in ..._EXX of the video name), set below
  effortCols <- c(Prey_num     = "Prey.number",         ### CHECK
                  Fasting_time = "Days.without.food")   ### CHECK


#####################
#  Read the data    #
#####################

  kin   <- read.csv(kinFile, stringsAsFactors = FALSE)
  flow  <- read.csv(flowFile, stringsAsFactors = FALSE)
  scan  <- read.csv(scanFile, stringsAsFactors = FALSE)
  sheet <- read.csv(sheetFile, stringsAsFactors = FALSE)

  # Same data-sheet rows as script 2 (equivVidRow refers to this filtered sheet)
  vidData <- sheet[which(sheet$Generation == 'F2' & sheet$Event != '-'), ]

  # Only videos that were processed
  kin <- kin[kin$processed %in% TRUE, ]
  if(dropCensored) {
    kin <- kin[!(kin$Tend_censored %in% TRUE), ]
  }
  print(paste0(nrow(kin), " processed videos"))


##### Scan values not already in the kinematics output (CSAs, crest depth, urohyal, premaxilla)
  # Match on family and animal ID the same way script 2 does (A prefix, periods instead of + and -)
  scan$Animal.ID <- as.character(scan$Animal.ID)
  noA <- which(substring(scan$Animal.ID, 1, 1) != 'A')
  scan$Animal.ID[noA] <- paste0('A', scan$Animal.ID[noA])
  scan$Animal.ID <- gsub('\\+', '.', scan$Animal.ID)
  scan$Animal.ID <- gsub('\\-', '.', scan$Animal.ID)

  scanRow <- match(paste(kin$family, kin$Animal.ID), paste(scan$Family, scan$Animal.ID))

  epaxialCol <- c(crest = "Epaxial_CSA_crest", flat = "Epaxial_CSA_flat")[[epaxialSection]]
  scanExtra  <- c(Epaxial_CSA = epaxialCol, Hypaxial_CSA = "Hypaxial_CSA", Sternohyoideus_CSA = "Sternohyoideus_CSA",
                  IOP_ligament_CSA = "IOP_ligament_CSA", Protractor_hyoideus_CSA = "Protractor_hyoideus_CSA",
                  Supraoccipital_crest_depth = "Supraoccipital_crest_depth", Urohyal_length = "Urohyal_length",
                  Length_premaxilla_asc_arm = "Length_premaxilla_asc_arm")
  for(nm in names(scanExtra)) {
    kin[[nm]] <- scan[[scanExtra[[nm]]]][scanRow]
  }
  print(paste0(sum(!is.na(scanRow)), " videos matched to a scan"))


##### Flow (one row per kinematic video, matched on vidName)
  flowRow <- match(kin$vidName, flow$vidName)
  # Column name in the flow output for one script 3 version
  flowName <- function(colNow, versionNow) paste0(colNow, if(versionNow == "") "" else paste0("_", versionNow))
  # Script 3 versions for the instantaneous values and for the strike means of a flow version ("a+b", or one version for both)
  flowParts <- function(flowNow) {
    partsNow <- strsplit(flowNow, "+", fixed = TRUE)[[1]]
    if(length(partsNow) == 1) partsNow <- c(partsNow, partsNow)
    if(flowNow == "") partsNow <- c("", "")
    if(length(partsNow) != 2) stop("Flow version must be one version or 'instant+means': ", flowNow)
    setNames(partsNow, c("instant", "means"))
  }


##### Effort inputs from the data sheet
  for(nm in names(effortCols)) {
    if(effortCols[[nm]] %in% names(vidData)) {
      kin[[nm]] <- vidData[[effortCols[[nm]]]][kin$equivVidRow]
    } else {
      kin[[nm]] <- NA
      print(paste0("Data sheet has no column '", effortCols[[nm]], "': ", nm, " is NA"))
    }
  }

  # Fasting_time must be a number. Cells that aren't (text, typos) become NA, and are listed so they can be fixed in the data sheet
  fastRaw <- trimws(as.character(kin$Fasting_time))
  kin$Fasting_time <- suppressWarnings(as.numeric(fastRaw))
  badFast <- !is.na(fastRaw) & fastRaw != "" & is.na(kin$Fasting_time)
  if(any(badFast)) {
    badCounts <- table(fastRaw[badFast])
    print(paste0(sum(badFast), " videos have a ", effortCols[["Fasting_time"]], " value that isn't a number (set to NA): ",
                 paste(paste0("'", names(badCounts), "' = ", badCounts), collapse = ", ")))
    print(head(kin$vidName[badFast], 10))
  }

  # Videos whose prey column flags a PIV problem (pivFlagValues): pressure and flow are removed in buildDataset
  flagCount <- table(tolower(trimws(as.character(kin$Prey_num)))[tolower(trimws(as.character(kin$Prey_num))) %in% tolower(trimws(pivFlagValues))])
  if(length(flagCount) > 0) {
    print(paste0(sum(flagCount), " videos with a PIV problem in the prey column (pressure and flow set to NA): ",
                 paste(names(flagCount), flagCount, sep = " = ", collapse = ", ")))
  }

  # Satiation = strike number for the fish (E01 = 1, E02 = 2, ...), from the end of the video name
    # The name can have one more part after it (take or segment number: _000000, _000001, _S0001, ...), which is skipped
    # Names that don't match: the Event column from the data sheet (E01, ...) instead
  kin$Satiation <- suppressWarnings(as.numeric(sub("^.*_E([0-9]+)(_[A-Za-z0-9]+)?$", "\\1", kin$vidName)))
  fromEvent <- is.na(kin$Satiation)
  if("Event" %in% names(kin)) {
    kin$Satiation[fromEvent] <- suppressWarnings(as.numeric(sub("^E", "", kin$Event[fromEvent])))
  }
  if(any(fromEvent & !is.na(kin$Satiation))) {
    print(paste0(sum(fromEvent & !is.na(kin$Satiation)), " videos: Satiation from the data sheet's Event column (no _E## in the name)"))
  }
  if(any(is.na(kin$Satiation))) {
    print(paste0(sum(is.na(kin$Satiation)), " videos have no strike number in the name or the Event column, so their Satiation is NA:"))
    print(head(kin$vidName[is.na(kin$Satiation)], 10))
  }


#######################################
#  Variables computed here            #
#######################################

  # Neurocranium mechanical advantage
  kin$Neurocranium_MA <- kin$Supraoccipital_crest_depth / kin$Neurocranium_length

  # <Neurocranium_Adj_vel(t)> = mean of Neurocranium_linear_vel(t) * MA / Epaxial_CSA, Tstart to peak cranial elevation
    # MA and CSA are constants for each fish, so this is exactly the plain mean times MA / CSA
  kin$Neurocranium_Adj_vel_mean <- kin$Neurocranium_linear_vel_mean * kin$Neurocranium_MA / kin$Epaxial_CSA

  # <Hyoid_Adj_vel(t)> = mean of Hyoid_vel(t) / (Hypaxial_CSA * Sternohyoideus_CSA), Tstart to peak depression
  kin$Hyoid_Adj_vel_mean <- kin$Hyoid_vel_mean / (kin$Hypaxial_CSA * kin$Sternohyoideus_CSA)


#######################################
#  Build one dataset per time basis   #
#######################################

  # Width-model suffix for the buccal-volume columns
  widthSuffixOf <- function(widthNow) c(main = "", noRetract = "_noRetract", noWidth = "_noWidth")[[widthNow]]

  # Variables that are calculated at a time: name in the SEM (with _tmax) = column stem in the outputs (time suffix added below)
  timeVarsKin <- c("R", "R_vel", "A", "A_vel", "Protrusion", "Protrusion_vel",
                   "Mandible_ang", "Mandible_angular_vel", "Mandible_depression", "Mandible_dep_vel", "Mandible_tip_disp",
                   "Neurocranium_rotation", "Neurocranium_angular_vel", "Neurocranium_tip_disp", "Neurocranium_linear_vel",
                   "Maxilla_ang", "Maxilla_angular_vel", "Maxilla_tip_disp", "Maxilla_linear_vel",
                   "Hyoid_depression", "Hyoid_vel",
                   "U_ram", "U_ram_acc",
                   "Oral_theta2", "Oral_theta3", "Oral_theta4", "Oral_KT", "Maxilla_ang_predicted",
                   "Opercular_theta2", "Opercular_theta3", "Opercular_theta4", "Opercular_KT", "Mandible_ang_predicted",
                   "Relative")
  timeVarsBuccal <- c("V_buccal_vel", "V_buccal_acc")                        # + width-model suffix
  timeVarsFlow   <- c("U_flow_ef_meas", "U_flow_ef_meas_acc", "U_flow_ff", "U_flow_ff_acc")   # + flow-version suffix

  # One value per strike (no time)
  strikeVars <- c("Ttotal", "Ttpg", "Relative_ttpg", "Kinetic_Synchronization", "Time_hyoid", "Time_cranial",
                  "R_t_mean", "A_mean", "A_Max", "Protrusion_t_mean", "Mandible_depression_t_mean", "Hyoid_depression_t_mean",
                  "Hyoid_Max", "Protrusion_Max",   # strike maxima, for the Max equations in script 5 (+ R_Max, Maxilla_tip_disp_Max below)
                  "Mandible_Max", "Mandible_depression_Max", "Neurocranium_Max",
                  "U_ram_mean", "Cov_R_Uram",
                  "Hyoid_vel_mean", "Neurocranium_linear_vel_mean", "Neurocranium_Adj_vel_mean", "Hyoid_Adj_vel_mean",
                  "Prey_num", "Fasting_time", "Satiation")
  strikeVarsBuccal <- c("V_buccal_Max", "U_flow_ff_predicted_mean_bigGape")   # + width-model suffix
  strikeVarsFlow   <- c("U_flow_ef_meas_mean", "U_flow_ff_mean", "U_flow_ff_mean_bigGape", "U_flow_ef_meas_mean_bigGape",
                        "Cov_R_Uflowefmeas", "Cov_R_Uflowff", "Cov_A_Uflowff", "Total_ingested_volume")   # + flow-version suffix

  # One value per fish (the same in every strike)
  fishVars <- c("Protrusion_ROM", "Mandible_ROM", "Hyoid_ROM", "Neurocranium_ROM",
                "Oral_4bar_KT", "Opercular_4bar_KT", "Oral_theta2_rest", "Oral_theta3_rest", "Opercular_theta2_rest", "Opercular_theta4_rest",
                "Oral_input_ratio", "Oral_coupler_ratio", "Oral_output_ratio",
                "Opercular_input_ratio", "Opercular_coupler_ratio", "Opercular_output_ratio", "IOP_link_length",
                "Mandible_length", "Mandible_input_length", "Maxilla_length", "Maxilla_palatine_length", "Oral_fixed_length",
                "Operculum_length", "RA_process_length", "Opercular_fixed_length", "Neurocranium_length",
                "W_jaw", "W_head", "D_head", "Ceratohyal_length",
                "Epaxial_CSA", "Hypaxial_CSA", "Sternohyoideus_CSA", "IOP_ligament_CSA", "Protractor_hyoideus_CSA",
                "Supraoccipital_crest_depth", "Neurocranium_MA", "Urohyal_length", "Length_premaxilla_asc_arm")

  # Flags kept for subsetting
  flagVars <- c("Tend_censored", "Mouth_open_corrected", "Hyoid_depth_source")

  # Stop if any expected column is missing from the inputs (a name changed in an earlier script)
  if(!(pressureVersion %in% c("corrected", "raw"))) {
    stop("pressureVersion must be 'corrected' or 'raw'")
  }
  if(!(widthModel %in% widthModels) || !(flowVersion %in% flowVersions)) {
    stop("widthModel and flowVersion must be in widthModels and flowVersions")
  }
  # Flow versions made only with some script 3 settings (the K versions): skipped when script 3 didn't make them
  optionalFlow <- grepl("K(\\+|$)", flowVersions)
  notMade <- optionalFlow & !sapply(flowVersions, function(f) all(paste0("U_flow_ef_meas_tmax_", flowParts(f)) %in% names(flow)))
  if(any(notMade)) {
    print(paste0("Flow versions skipped (not in script 3's output; set mullerExponentAlt there): ", paste(flowVersions[notMade], collapse = ", ")))
    flowVersions <- flowVersions[!notMade]
  }
  missingKin  <- setdiff(c(as.vector(outer(timeVarsKin, c("_tmax", "_tmaxVel"), paste0)),
                           unlist(lapply(widthModels, function(w) c(as.vector(outer(timeVarsBuccal, c("_tmax", "_tmaxVel"), function(a, b) paste0(a, b, widthSuffixOf(w)))),
                                                                     paste0(strikeVarsBuccal, widthSuffixOf(w))))),
                           "Hyoid_tmax_filled", "Hyoid_tmaxVel_filled", "tmax", "tmaxVel", "Tstart",
                           setdiff(strikeVars, names(effortCols)), fishVars, flagVars),
                         names(kin))
  missingFlow <- setdiff(c(unlist(lapply(flowVersions, function(f) c(flowName(as.vector(outer(timeVarsFlow, c("_tmax", "_tmaxVel"), paste0)), flowParts(f)[["instant"]]),
                                                                     flowName(strikeVarsFlow, flowParts(f)[["means"]])))), "Max_sub_pressure", "Max_sub_pressure_corr",
                           "static_node_zR_tmax", "static_node_zR_tmaxVel", "pressure_node_zR"),
                         names(flow))
  if(length(missingKin) > 0 || length(missingFlow) > 0) {
    stop("Columns not found. Kinematics: ", paste(missingKin, collapse = ", "), ". Flow: ", paste(missingFlow, collapse = ", "))
  }


  buildDataset <- function(timeNow, widthNow, flowNow) {

    widthSuffix <- widthSuffixOf(widthNow)
    out <- data.frame(vidName = kin$vidName, Fish_ID = kin$Unique.ID, Family = kin$family, Animal.ID = kin$Animal.ID,
                      timeBasis = timeNow, widthModel = widthNow, flowVersion = flowNow, stringsAsFactors = FALSE)

    # Instantaneous values, saved under the _tmax name whichever time this is
      # "peak": each variable's maximum over the strike (NA where the column doesn't exist; script 5 recalculates identities)
    colOr <- function(df, colNow, rows = NULL) {
      if(!(colNow %in% names(df))) return(rep(NA_real_, nrow(out)))
      if(is.null(rows)) df[[colNow]] else df[[colNow]][rows]
    }
    kinPeakName  <- c(A = "A_Max", Hyoid_depression = "Hyoid_Max", Relative = "Relative_tmax")   # peaks saved under other names
    flowPeakName <- c(U_flow_ef_meas = "U_flow_ef_meas_Max", U_flow_ef_meas_acc = "U_flow_ef_meas_acc_Max")
    for(v in timeVarsKin) {
      colNow <- if(timeNow != "peak") paste0(v, "_", timeNow) else if(v %in% names(kinPeakName)) kinPeakName[[v]] else paste0(v, "_peak")
      out[[paste0(v, "_tmax")]] <- colOr(kin, colNow)
    }
    for(v in timeVarsBuccal) out[[paste0(v, "_tmax")]] <- colOr(kin, paste0(v, "_", if(timeNow == "peak") "peak" else timeNow, widthSuffix))
    for(v in timeVarsFlow) {
      stemNow <- if(timeNow != "peak") paste0(v, "_", timeNow) else if(v %in% names(flowPeakName)) flowPeakName[[v]] else paste0(v, "_peak")
      out[[paste0(v, "_tmax")]] <- colOr(flow, flowName(stemNow, flowParts(flowNow)[["instant"]]), flowRow)
    }
    out$Hyoid_tmax_filled <- if(timeNow == "peak") NA else kin[[paste0("Hyoid_", timeNow, "_filled")]]
    # Time of tmax (or tmaxVel) from Tstart (s), for Relative_tmax = tmax / T in script 5 ("peak": tmax, the pressure peak)
    out$tmax_from_Tstart <- kin[[if(timeNow == "peak") "tmax" else timeNow]] - kin$Tstart

    # Per strike and per fish
    for(v in strikeVars)       out[[v]] <- kin[[v]]
    # Strike maxima script 2 doesn't save under a _Max name
      # R_Max: from A_Max (= pi W_jaw / 2 x the largest gape), so the A_Max identity in script 5 holds exactly; the spline peak without W_jaw
    out$R_Max <- ifelse(!is.na(kin$A_Max) & !is.na(kin$W_jaw), kin$A_Max / (pi * kin$W_jaw / 2), colOr(kin, "R_peak"))
    out$Maxilla_tip_disp_Max <- colOr(kin, "Maxilla_tip_disp_peak")
    for(v in strikeVarsBuccal) out[[v]] <- kin[[paste0(v, widthSuffix)]]

    # Buccal volume shrinking at tmaxVel: no buccal-volume values for this video (dropNegBuccalVel)
    out$buccal_excluded <- dropNegBuccalVel & (kin[[paste0("V_buccal_vel_tmaxVel", widthSuffix)]] < 0) %in% TRUE
    out[out$buccal_excluded, c(paste0(timeVarsBuccal, "_tmax"), strikeVarsBuccal)] <- NA
    for(v in strikeVarsFlow)   out[[v]] <- flow[[flowName(v, flowParts(flowNow)[["means"]])]][flowRow]
    for(v in fishVars)         out[[v]] <- kin[[v]]
    out$Neurocranium_tip_disp_Max <- out$Neurocranium_Max * out$Neurocranium_length   # as Neurocranium_tip_disp_tmax (identity in script 5)
    for(v in flagVars)         out[[v]] <- kin[[v]]

    # Outcome 1 (Pa, negative = sub-ambient). Doesn't depend on the time basis
    out$Max_sub_pressure_raw <- flow$Max_sub_pressure[flowRow]
    out$Max_sub_pressure     <- if(pressureVersion == "corrected") flow$Max_sub_pressure_corr[flowRow] else out$Max_sub_pressure_raw
    # The peak corrected with script 3's second decay (pressureDecayAlt), for script 5's comparePressureK (NA if script 3 didn't make it)
    out$Max_sub_pressure_K   <- colOr(flow, "Max_sub_pressure_corrK", flowRow)

    # Distance of the static velocity node from the mouth at tmax (or tmaxVel), in multiples of the gape radius R
      # (+ = in front of the mouth). The same in every flow version; a covariate for the uncorrected static flow in script 5
    out$static_node_zR_tmax <- flow[[paste0("static_node_zR_", if(timeNow == "peak") "tmaxVel" else timeNow)]][flowRow]   # peak: at peak velocity

    # Distance of the pressure peak node from the mouth (multiples of R), at the pressure peak. A covariate for pressure in script 5
    out$pressure_node_zR <- flow$pressure_node_zR[flowRow]

    # SEM names for the per-strike variables
    names(out)[names(out) == "Ttotal"] <- "T"
    names(out)[names(out) == "Ttpg"]   <- "ttpg"

    # Prey_num as an ordered factor (levels in preyLevels order, not alphabetical)
      # Blank cells are NA
    preyNow <- tolower(trimws(as.character(out$Prey_num)))
    preyNow[preyNow == ""] <- NA
    levelsNow <- tolower(trimws(preyLevels))
    flagsNow  <- tolower(trimws(pivFlagValues))
    if(any(!(na.omit(preyNow) %in% c(levelsNow, flagsNow)))) {
      stop("Prey_num values not in preyLevels or pivFlagValues: ", paste(setdiff(unique(na.omit(preyNow)), c(levelsNow, flagsNow)), collapse = ", "))
    }

    # PIV problem flagged in the prey column: no pressure or flow for this video, and no prey count
    out$PIV_flagged <- preyNow %in% flagsNow
    pivCols <- c(paste0(timeVarsFlow, "_tmax"), strikeVarsFlow, "Max_sub_pressure", "Max_sub_pressure_raw", "Max_sub_pressure_K", "static_node_zR_tmax", "pressure_node_zR")
    out[out$PIV_flagged, pivCols] <- NA
    preyNow[out$PIV_flagged] <- NA
    out$Prey_num <- factor(preyNow, levels = levelsNow, ordered = TRUE)

    # One row per fish: the whole row of the fish's strike with the most negative Max_sub_pressure
      # Fish with no pressure in any strike: the strike with the largest A_Max, or left out (fishNoPressure)
      # bestStrikeBasis records which. n_strikes = how many strikes the fish has
    if(unitOfAnalysis == "fish") {
      pickRow <- sapply(split(seq_len(nrow(out)), out$Fish_ID), function(rowsNow) {
        if(any(!is.na(out$Max_sub_pressure[rowsNow]))) {
          return(rowsNow[which.min(out$Max_sub_pressure[rowsNow])])
        }
        if(fishNoPressure == "largestGape" && any(!is.na(out$A_Max[rowsNow]))) {
          return(rowsNow[which.max(out$A_Max[rowsNow])])
        }
        NA_integer_
      })
      fishOut <- out[na.omit(pickRow), ]
      fishOut$bestStrikeBasis <- ifelse(!is.na(fishOut$Max_sub_pressure), "pressure", "largestGape")
      fishOut$n_strikes <- as.vector(table(out$Fish_ID)[fishOut$Fish_ID])
      if(any(is.na(pickRow))) {
        print(paste0(sum(is.na(pickRow)), " fish left out (no strike with pressure", if(fishNoPressure == "largestGape") " or A_Max", ")"))
      }
      return(fishOut)
    }
    out
  }


  # Effort (latent variable) is not made here: it's fitted in script 5 (see 5_SEM_notes.md)
  # Every time basis x width model x flow version, named "timeBasis|widthModel|flowVersion"
  if("peak" %in% timeBases && !("R_peak" %in% names(kin))) {
    print("The kinematics output has no _peak columns (rerun script 2 in runMode 'all'): the peak time basis is skipped")
    timeBases <- setdiff(timeBases, "peak")
  }
  combos <- expand.grid(timeBasis = timeBases, widthModel = widthModels, flowVersion = flowVersions, stringsAsFactors = FALSE)
  allDatasets <- list()
  for(k in seq_len(nrow(combos))) {
    allDatasets[[paste(combos$timeBasis[k], combos$widthModel[k], combos$flowVersion[k], sep = "|")]] <-
      buildDataset(combos$timeBasis[k], combos$widthModel[k], combos$flowVersion[k])
  }

  # The chosen combination, for the checks and the CSVs
  datasets <- lapply(setNames(timeBases, timeBases), function(tNow) allDatasets[[paste(tNow, widthModel, flowVersion, sep = "|")]])


############
#  Checks  #
############

  # Rows that have every variable of each path equation (SEM_Equations, path equations only)
    # Low numbers show which equation is limited by missing data
  semEquations <- list(
    Pressure            = c("Max_sub_pressure", "U_flow_ff_tmax", "U_ram_tmax", "R_tmax", "U_flow_ff_acc_tmax", "U_ram_acc_tmax", "R_vel_tmax"),
    Ingested_volume     = c("Total_ingested_volume", "T", "A_mean", "U_flow_ff_mean", "Cov_A_Uflowff"),
    Gape_tmax           = c("R_tmax", "Protrusion_tmax", "Mandible_depression_tmax"),
    Gape_vel_tmax       = c("R_vel_tmax", "Protrusion_vel_tmax", "Mandible_dep_vel_tmax"),
    Gape_mean           = c("R_t_mean", "Protrusion_t_mean", "Mandible_depression_t_mean"),
    Protrusion_tmax     = c("Protrusion_tmax", "Neurocranium_tip_disp_tmax", "Maxilla_tip_disp_tmax"),
    Protrusion_vel_tmax = c("Protrusion_vel_tmax", "Neurocranium_linear_vel_tmax", "Maxilla_linear_vel_tmax"),
    Protrusion_mean     = c("Protrusion_t_mean", "Relative_ttpg", "Protrusion_ROM"),
    Mandible_dep_mean   = c("Mandible_depression_t_mean", "Relative_ttpg", "Mandible_ROM", "ttpg", "Mandible_dep_vel_tmax", "Hyoid_depression_t_mean", "Protractor_hyoideus_CSA", "Mandible_length"),
    Flow_tmax           = c("U_flow_ef_meas_tmax", "Relative_tmax", "Kinetic_Synchronization", "V_buccal_vel_tmax", "A_tmax"),
    Flow_acc_tmax       = c("U_flow_ef_meas_acc_tmax", "Relative_tmax", "Kinetic_Synchronization", "V_buccal_acc_tmax", "V_buccal_vel_tmax", "A_vel_tmax", "A_tmax"),
    Flow_mean           = c("U_flow_ef_meas_mean", "Kinetic_Synchronization", "U_flow_ff_predicted_mean_bigGape"),
    Mandible_angvel     = c("Mandible_angular_vel_tmax", "Neurocranium_angular_vel_tmax", "Opercular_KT_tmax", "IOP_ligament_CSA", "Hyoid_vel_tmax", "Protractor_hyoideus_CSA", "Mandible_length"),
    Hyoid_vel_tmax      = c("Hyoid_vel_tmax", "Hypaxial_CSA", "Sternohyoideus_CSA", "Urohyal_length"),            # + Effort (script 5)
    Hyoid_dep_tmax      = c("Hyoid_depression_tmax", "Relative_tmax", "Hyoid_ROM"),
    Neuro_linvel_tmax   = c("Neurocranium_linear_vel_tmax", "Epaxial_CSA", "Neurocranium_MA"),                     # + Effort (script 5)
    Neuro_rotation_tmax = c("Neurocranium_rotation_tmax", "Relative_tmax", "Neurocranium_ROM"),
    Maxilla_angvel      = c("Maxilla_angular_vel_tmax", "Mandible_angular_vel_tmax", "Oral_KT_tmax"),
    Maxilla_ang_check   = c("Maxilla_ang_tmax", "Maxilla_ang_predicted_tmax"),
    Mandible_ang_check  = c("Mandible_ang_predicted_tmax", "Opercular_theta4_tmax", "Opercular_theta4_rest", "IOP_ligament_CSA", "Hyoid_depression_tmax", "Protractor_hyoideus_CSA", "Mandible_length"),
    Protrusion_ROM      = c("Protrusion_ROM", "Length_premaxilla_asc_arm"),
    Effort_indicators   = c("Neurocranium_Adj_vel_mean", "Hyoid_Adj_vel_mean", "T"),                               # all three (FIML in script 5 can use fewer)
    Effort_covariates   = c("Prey_num", "Fasting_time", "Satiation"))

  for(tNow in names(datasets)) {
    d <- datasets[[tNow]]
    print(paste0("===== ", tNow, ": ", nrow(d), " rows, ", length(unique(d$Fish_ID)), " fish"))
    if(dropNegBuccalVel) {
      print(paste0(sum(d$buccal_excluded), " videos without buccal-volume values (volume shrinking at tmaxVel, width model ", widthModel, ")"))
    }
    coverage <- sapply(semEquations, function(v) sum(complete.cases(d[, v])))
    print(data.frame(rows_with_all_variables = coverage))
    ### CHECK: a 0 above is expected when an equation uses a column that is all NA here. Those are identities
    # (or, for the peak basis, columns with no peak in script 2) that script 5 recalculates from their parents.
    allNA <- unique(unlist(semEquations))
    allNA <- allNA[allNA %in% names(d) & sapply(allNA, function(v) v %in% names(d) && all(is.na(d[[v]])))]
    if(length(allNA) > 0) {
      print(paste0("All NA here (recalculated in script 5 if an identity): ", paste(allNA, collapse = ", ")))
    }

    # Ranges, to catch unit problems (e.g. flow still in m/s, pressure positive)
    print(summary(d[, c("Max_sub_pressure", "U_flow_ef_meas_tmax", "U_flow_ff_tmax", "U_ram_tmax", "R_tmax", "A_tmax",
                        "V_buccal_vel_tmax", "Total_ingested_volume", "T")]))
    if(any(d$Max_sub_pressure > 0, na.rm = TRUE)) {
      print(paste0(sum(d$Max_sub_pressure > 0, na.rm = TRUE), " rows with positive Max_sub_pressure (expected negative)"))
    }
  }


##########
#  Save  #
##########

  fishTag <- if(unitOfAnalysis == "fish") "_byFish" else ""
  for(tNow in names(datasets)) {
    fileNow <- paste0(outPrefix, "_", tNow, fishTag, ".csv")
    write.csv(datasets[[tNow]], fileNow, row.names = FALSE)
    print(paste0("Saved: ", fileNow))
  }
  # All combinations with their column types (Prey_num ordered), for script 5
  saveRDS(allDatasets, paste0(outPrefix, fishTag, ".rds"))
  print(paste0("Saved: ", outPrefix, fishTag, ".rds (", length(allDatasets), " datasets: ", paste(names(allDatasets), collapse = ", "), ")"))
  print(paste0("Settings: unit = ", unitOfAnalysis, ", pressure = ", pressureVersion, ", epaxial = ", epaxialSection, ", width model = ", widthModel,
               ", flow version = ", if(flowVersion == "") "script 3 mainVersion" else flowVersion, " (CSVs)", ", drop censored = ", dropCensored))
