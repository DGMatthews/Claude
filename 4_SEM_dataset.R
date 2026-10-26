
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
      # e.g. "static+staticRaw": the aperture-corrected static node at tmax (where the correction matters for pressure),
        # and the uncorrected static node for the strike means. The corrected node only covers part of each strike
        # (frames within maxZoverR of the mouth), so its strike means exist for few videos and its ingested volume is extrapolated
  widthModels  <- c("main", "noRetract", "noWidth")   ### CHANGE  any of "main" (bar model), "noRetract", "noWidth"
  flowVersions <- c("static+staticRaw", "static", "staticRaw", "static+mouthFrontCorr", "mouthFrontCorr")   ### CHANGE

  # The combination saved to the CSVs and printed in the checks (also the default in script 5)
  widthModel  <- "main"            ### CHANGE
  flowVersion <- "static+staticRaw"   ### CHANGE

  # Videos where the gape never dropped back below the threshold (Tend_censored). Their duration-based values are already NA
  dropCensored <- FALSE        ### CHANGE  TRUE drops these videos entirely

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
  if(!(widthModel %in% widthModels) || !(flowVersion %in% flowVersions)) {
    stop("widthModel and flowVersion must be in widthModels and flowVersions")
  }
  missingKin  <- setdiff(c(as.vector(outer(timeVarsKin, c("_tmax", "_tmaxVel"), paste0)),
                           unlist(lapply(widthModels, function(w) c(as.vector(outer(timeVarsBuccal, c("_tmax", "_tmaxVel"), function(a, b) paste0(a, b, widthSuffixOf(w)))),
                                                                     paste0(strikeVarsBuccal, widthSuffixOf(w))))),
                           "Hyoid_tmax_filled", "Hyoid_tmaxVel_filled", "tmax", "tmaxVel", "Tstart",
                           setdiff(strikeVars, names(effortCols)), fishVars, flagVars),
                         names(kin))
  missingFlow <- setdiff(c(unlist(lapply(flowVersions, function(f) c(flowName(as.vector(outer(timeVarsFlow, c("_tmax", "_tmaxVel"), paste0)), flowParts(f)[["instant"]]),
                                                                     flowName(strikeVarsFlow, flowParts(f)[["means"]])))), "Max_sub_pressure",
                           "static_node_zR_tmax", "static_node_zR_tmaxVel"),
                         names(flow))
  if(length(missingKin) > 0 || length(missingFlow) > 0) {
    stop("Columns not found. Kinematics: ", paste(missingKin, collapse = ", "), ". Flow: ", paste(missingFlow, collapse = ", "))
  }


  buildDataset <- function(timeNow, widthNow, flowNow) {

    widthSuffix <- widthSuffixOf(widthNow)
    out <- data.frame(vidName = kin$vidName, Fish_ID = kin$Unique.ID, Family = kin$family, Animal.ID = kin$Animal.ID,
                      timeBasis = timeNow, widthModel = widthNow, flowVersion = flowNow, stringsAsFactors = FALSE)

    # Instantaneous values, saved under the _tmax name whichever time this is
    for(v in timeVarsKin)    out[[paste0(v, "_tmax")]] <- kin[[paste0(v, "_", timeNow)]]
    for(v in timeVarsBuccal) out[[paste0(v, "_tmax")]] <- kin[[paste0(v, "_", timeNow, widthSuffix)]]
    for(v in timeVarsFlow)   out[[paste0(v, "_tmax")]] <- flow[[flowName(paste0(v, "_", timeNow), flowParts(flowNow)[["instant"]])]][flowRow]
    out$Hyoid_tmax_filled <- kin[[paste0("Hyoid_", timeNow, "_filled")]]
    # Time of tmax (or tmaxVel) from Tstart (s), for Relative_tmax = tmax / T in script 5
    out$tmax_from_Tstart <- kin[[timeNow]] - kin$Tstart

    # Per strike and per fish
    for(v in strikeVars)       out[[v]] <- kin[[v]]
    for(v in strikeVarsBuccal) out[[v]] <- kin[[paste0(v, widthSuffix)]]
    for(v in strikeVarsFlow)   out[[v]] <- flow[[flowName(v, flowParts(flowNow)[["means"]])]][flowRow]
    for(v in fishVars)         out[[v]] <- kin[[v]]
    for(v in flagVars)         out[[v]] <- kin[[v]]

    # Outcome 1 (Pa, negative = sub-ambient). Doesn't depend on the time basis
    out$Max_sub_pressure <- flow$Max_sub_pressure[flowRow]

    # Distance of the static velocity node from the mouth at tmax (or tmaxVel), in multiples of the gape radius R
      # (+ = in front of the mouth). The same in every flow version; a covariate for the uncorrected static flow in script 5
    out$static_node_zR_tmax <- flow[[paste0("static_node_zR_", timeNow)]][flowRow]

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
    pivCols <- c(paste0(timeVarsFlow, "_tmax"), strikeVarsFlow, "Max_sub_pressure", "static_node_zR_tmax")
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
  combos <- expand.grid(timeBasis = c("tmax", "tmaxVel"), widthModel = widthModels, flowVersion = flowVersions, stringsAsFactors = FALSE)
  allDatasets <- list()
  for(k in seq_len(nrow(combos))) {
    allDatasets[[paste(combos$timeBasis[k], combos$widthModel[k], combos$flowVersion[k], sep = "|")]] <-
      buildDataset(combos$timeBasis[k], combos$widthModel[k], combos$flowVersion[k])
  }

  # The chosen combination, for the checks and the CSVs
  datasets <- list(tmax    = allDatasets[[paste("tmax", widthModel, flowVersion, sep = "|")]],
                   tmaxVel = allDatasets[[paste("tmaxVel", widthModel, flowVersion, sep = "|")]])


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
    coverage <- sapply(semEquations, function(v) sum(complete.cases(d[, v])))
    print(data.frame(rows_with_all_variables = coverage))

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
  print(paste0("Settings: unit = ", unitOfAnalysis, ", epaxial = ", epaxialSection, ", width model = ", widthModel,
               ", flow version = ", if(flowVersion == "") "script 3 mainVersion" else flowVersion, " (CSVs)", ", drop censored = ", dropCensored))
