
#############################################################################################
# Flow and pressure from the PIV output                                                     #
#############################################################################################
  # Reads the PIV pressure and velocity traces for each kinematic video, cleans and smooths them,
    # and calculates the flow and pressure variables for the SEM
  # Run after 2_ACxTRC_Kinematic_Measurements_Sept2026.R: it uses that script's output (strike timing, gape, ram speed)
  # PIV frame 0 = kinematic frame 0 (time 0), both 500 fps. In R, kinematic frame index k (1, 2, ...) = PIV frame k - 1

  # Three versions of the measured flow speed (U_flow_ef_meas), each with its own set of outputs (column suffix):
    # _mouthFront: velocity sampled immediately in front of the mouth in every frame (mouth_front_velocity_trace files). Follows the mouth, but noisier
                     # It is the maximum within the jaw triangle plus a box 1.5 mm in front, so the point can be up to ~1.5 mm from the aperture
    # _mouthFrontCorr: _mouthFront corrected to the aperture with Muller's flow model (below), frame by frame,
                     # using where the maximum was (position columns in the mouth-front files), smoothed over time first
                     # Points at or behind the aperture plane (z <= 0) aren't corrected (factor 1)
    # _static:     velocity at the fixed node where peak velocity happened (peak_node_traces files), corrected to the mouth aperture with
                     # Muller et al.'s flow model: U_mouth = U_node * (z^2 + R^2)^1.5 / R^3
                     # z = distance from the mouth centre to the node along the strike axis, R = gape radius, both from the PIV mouth landmarks
                     # Only frames where the node is close to the mouth are used (maxZoverR, maxOffAxisOverR), since the correction grows quickly with distance
                     # (the same limits apply to _mouthFrontCorr)
    # _staticRaw:  the fixed-node velocity without the correction (reference)
  # The version chosen in mainVersion is also saved without a suffix, under the names the SEM uses (U_flow_ef_meas_tmax, U_flow_ff_mean, ...)

  # Videos included: those in final_results_consistency_filtered.xlsx (the kinematics script reads tmax, tmaxVel and the peak values from it)
    # final_results_all_passed.xlsx only supplies where the peak velocity and pressure nodes were (px)
  # Pressure: Max_sub_pressure = the PIV pipeline's peak pressure (Pa, negative = sub-ambient), only for videos whose pressure analysis passed
    # Max_sub_pressure_corr: the same peak corrected to the aperture with Muller's flow model
      # Muller's velocity profile U(z) = U R^3 / (z^2 + R^2)^1.5 comes from the potential U R (1 - z / sqrt(z^2 + R^2)) on the strike axis,
      # so the unsteady pressure term d/dt(U R) falls with distance as g(s) = 1 - s / sqrt(1 + s^2), s = z / R
      # The 1/2 rho u^2 term falls as (1 + s^2)^-3 (3% of its aperture value at 1.5R), so at the node the pressure is almost all the
      # unsteady term, and p_aperture = p_node / g(s). Same distance limits as the velocity (maxZoverR, maxOffAxisOverR, node in front)
  # Flow is a speed (magnitude, m/s from the PIV, converted to mm/s), assumed to be along the strike axis like U_ram
  # U_flow_ff = U_flow_ef_meas + U_ram (fish frame: speed of water entering the mouth)
    # The measured flow is water moving toward the fish and ram is the fish moving toward the water, so the two add

  # Output: flowMeasurements.csv (one row per kinematic video, matched on vidName) and flowTimeSeries.rds

require(readxl)
require(imputeTS)

rm(list=ls())


##############
#  Settings  #
##############

  # PIV output
  pathPIVfinal     <- 'N:/TRCxAC/Hybrid feeding project/HSV/PIV/PIV Output/Pressure-velocity peak pipeline v2 (velocity-anchored)(0.4 pressure shape)/Final export'
  pivFilteredFile  <- file.path(pathPIVfinal, 'final_results_consistency_filtered.xlsx')   # decides which videos are included (same file the kinematics uses for tmax)
  pivResultsFile   <- file.path(pathPIVfinal, 'final_results_all_passed.xlsx')             # only used for the peak node positions (px)
  pathRawTraces    <- file.path(pathPIVfinal, 'Raw traces')
  pathMouthFront   <- 'N:/TRCxAC/Hybrid feeding project/HSV/PIV/PIV Output/Pressure-velocity peak pipeline v2 (velocity-anchored)(0.4 pressure shape)/Mouth front velocity traces/Traces'
  pathMouthMarks   <- 'N:/TRCxAC/Hybrid feeding project/HSV/PIV/Mouth Landmarks'

  # File name endings (after the PIV video name, e.g. "IL5-PIV_ACxTRC_2024.01.21_F2_A14.46_E05_000000")   ### CHECK once the file lists are in
  suffixPeakTraces <- '_peak_node_traces.csv'
  suffixMouthFront <- '_mouth_front_velocity_trace.csv'
  suffixMouthMarks <- '.mouthlandmarks.csv'

  # Columns in the mouth-front files with the position (px, PIV image) of the maximum in each frame, for _mouthFrontCorr
  frontXcol <- 'mouth_front_velocity_x_px'
  frontYcol <- 'mouth_front_velocity_y_px'

  # Kinematics output (from 2_ACxTRC_Kinematic_Measurements_Sept2026.R)
  pathData      <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data'
  kinFinalData  <- file.path(pathData, 'kinematicsAll_sleapOutput.csv')
  kinTimeSeries <- file.path(pathData, 'kinematicsTimeSeries.rds')
  # Only needed if the kinematics output has no PIV.filename column (runs from before it was added)
  dataSheet     <- file.path(pathData, 'ACxTRCfeedingdata_Oct2026_FIXED.csv')

  # Where to save
  outFlow       <- file.path(pathData, 'flowMeasurements.csv')
  outFlowTS     <- file.path(pathData, 'flowTimeSeries.rds')
  outFlowSkip   <- file.path(pathData, 'flowSkipped.csv')

  # Which velocity version is saved under the SEM names (no suffix)
  mainVersion <- "static"   ### CHANGE  "static", "staticRaw", "mouthFrontCorr" or "mouthFront"
    # static: the hand-checked fixed-node trace, corrected to the aperture. The uncorrected version (staticRaw) is also saved

  # Muller correction of the static node: only frames where the node is within these distances of the mouth (multiples of the gape radius R)
  maxZoverR       <- 2     ### CHANGE  along the strike axis. The correction factor is 2.8 at z = R, 5.9 at 1.5R, 11 at 2R (was 1.5)
  maxOffAxisOverR <- 1     ### CHANGE  sideways from the strike axis (Muller's model is for points on the axis)

  # Cleaning and smoothing, as in the kinematics script
    # Hampel filter: points more than hampelThreshold MADs from the median of a (2 * hampelHalfWindow + 1)-frame window are set to NA
    # Gaps up to maxGap frames are filled linearly, then a smoothing spline with lambda = lambdaRef * (spanRef / span)^3
  hampelHalfWindow  <- 3
  hampelThreshold   <- 3
  madFloorVel       <- 2      ### CHANGE  mm/s. Smallest MAD the Hampel filter uses for velocity (the PIV noise is about 5 mm/s)
  madFloorPressure  <- 0.2    ### CHANGE  Pa
  madFloorLandmark  <- 1      # px, same as the kinematics
  maxGap            <- 3
  # Velocity and pressure have their own lambda: 3b_Flow_smoothing_check.R showed 5e-5 was a little weak for velocity
    # (residuals still following frame-to-frame noise) and too strong for pressure (peak flattened). Recheck with 3b after changing
  lambdaRefVelocity <- 1.5e-4   ### CHANGE  velocity traces (3 x the old 5e-5)
  lambdaRefPressure <- 5e-6     ### CHANGE  pressure trace (0.1 x the old 5e-5). Only used for Pressure_tmax / tmaxVel; Max_sub_pressure is the PIV peak
  lambdaRefMarks    <- 5e-5   ### CHANGE  PIV mouth landmarks
  lambdaRefFrontPos <- 5e-5   ### CHANGE  position of the mouth-front maximum (px), smoothed before the Muller factor so the factor
                              # doesn't jump when the maximum moves to another PIV node. Stronger = smoother correction through the strike
  spanRef           <- 73

  # Which videos get strike means, covariances and total ingested volume (each uses the strike frames that have flow)
    # "beforePeak": the minFramesBeforePeak frames up to and including tmaxVel all have flow, however much of the rest of the strike is missing
      # (the static node loses frames once the fish passes it, and the static velocities were checked by hand)
    # "fraction": at least minStrikeFrac of the Tstart-Tend frames have flow (also used for videos with no tmaxVel)
  coverageRule        <- "beforePeak"   ### CHANGE  "beforePeak" or "fraction"
  minFramesBeforePeak <- 15             ### CHANGE  frames (30 ms at 500 fps)
  minStrikeFrac       <- 0.8            ### CHANGE  only used with coverageRule = "fraction"

  # Big-gape mean: frames with A >= uFlowAreaFrac * A_Max (same as the kinematics script)
  uFlowAreaFrac <- 0.5


#####################
#  Define functions #
#####################

    # Hampel filter: returns x with outlying points set to NA
    hampelNA <- function(x, halfWindow, threshold, madFloor) {
      flagged <- c()
      for(j in seq_along(x)) {
        if(is.na(x[j])) {
          next
        }
        pointsNow <- x[max(1, j - halfWindow):min(length(x), j + halfWindow)]
        madNow    <- max(mad(pointsNow, na.rm = TRUE), madFloor)
        if(abs(x[j] - median(pointsNow, na.rm = TRUE)) > threshold * madNow) {
          flagged <- c(flagged, j)
        }
      }
      x[flagged] <- NA
      x
    }


    # Clean and smooth one trace (one value per kinematic frame, NA where missing)
      # Returns the smoothed value and first derivative in every frame, NA outside the tracked frames and inside gaps longer than maxGap,
        # plus the spline (for values at tmax / tmaxVel)
    smoothTrace <- function(y, timeNow, madFloor, lambdaRefNow) {

      out <- list(y = rep(NA_real_, length(y)), dy = rep(NA_real_, length(y)), fit = NULL, usable = rep(FALSE, length(y)))

      y <- hampelNA(y, hampelHalfWindow, hampelThreshold, madFloor)
      if(sum(!is.na(y)) < 10) {
        return(out)
      }

      # Fill short gaps between the first and last values only
      firstNow <- min(which(!is.na(y)))
      lastNow  <- max(which(!is.na(y)))
      inside   <- firstNow:lastNow
      y[inside] <- na_interpolation(y[inside], option = "linear", maxgap = maxGap)

      goodNow <- which(!is.na(y))
      spanNow <- max(goodNow) - min(goodNow)
      fitNow  <- smooth.spline(x = timeNow[goodNow], y = y[goodNow], lambda = lambdaRefNow * (spanRef / spanNow)^3)

      out$usable[goodNow] <- TRUE
      out$y[goodNow]  <- predict(fitNow, x = timeNow[goodNow], deriv = 0)$y
      out$dy[goodNow] <- predict(fitNow, x = timeNow[goodNow], deriv = 1)$y
      out$fit <- fitNow
      out
    }


    # Smoothed value (deriv 0) or rate (deriv 1) at time tNow. NA if tNow is missing or falls in a frame without a usable value
    traceAt <- function(sm, tNow, timeNow, deriv) {
      if(is.null(sm$fit) || is.na(tNow)) {
        return(NA_real_)
      }
      kNow <- which.min(abs(timeNow - tNow))
      if(!sm$usable[kNow] || tNow < min(timeNow[sm$usable]) || tNow > max(timeNow[sm$usable])) {
        return(NA_real_)
      }
      predict(sm$fit, x = tNow, deriv = deriv)$y
    }


    # Put a PIV trace onto the kinematic frames: value for PIV time t (ms) goes to kinematic frame round(t / 1000 * frameRate) + 1
    toKinFrames <- function(time_ms, value, nFrames, frameRate) {
      out <- rep(NA_real_, nFrames)
      kNow <- round(time_ms / 1000 * frameRate) + 1
      keep <- !is.na(kNow) & !is.na(value) & kNow >= 1 & kNow <= nFrames
      out[kNow[keep]] <- value[keep]
      out
    }


    # Matching key for PIV video names, so small naming differences between the PIV files and the data sheet still match
      # Removes the take number at the end (_000000, _000001, ...), collapses doubled underscores ("__A9.10.02" = "_A9.10.02"),
      # then applies the fixes in pivNameFixes (old text = new text) for typos that can't be caught by a rule
    pivNameFixes <- c("_2.8_11.11_" = "_8.11.11_")   ### CHANGE  keep the same as in the kinematics script
    pivKey <- function(x) {
      x <- sub("_[0-9]{6}$", "", as.character(x))
      x <- gsub("_+", "_", x)
      for(k in seq_along(pivNameFixes)) {
        x <- gsub(names(pivNameFixes)[k], pivNameFixes[[k]], x, fixed = TRUE)
      }
      x
    }


    # Add a row to the skipped log
    addSkip <- function(skippedLog, vidName, section, reason) {
      rbind(skippedLog, data.frame(vidName = vidName, section = section, reason = reason))
    }


##################
#  Read the data #
##################

  kin   <- read.csv(kinFinalData, stringsAsFactors = FALSE)
  kinTS <- readRDS(kinTimeSeries)

  pivFiltered <- read_excel(pivFilteredFile)
  pivFiltered$videoKey <- pivKey(pivFiltered$video)
  pivResults  <- read_excel(pivResultsFile)
  pivResults$videoKey  <- pivKey(pivResults$video)

  # PIV name for each kinematic video
  if(!("PIV.filename" %in% names(kin))) {
    vidDataTemp <- read.csv(dataSheet, stringsAsFactors = FALSE)
    vidData     <- vidDataTemp[which(vidDataTemp$Generation == 'F2' & vidDataTemp$Event != '-'), ]   # same filter as the kinematics script
    kin$PIV.filename <- as.character(vidData$PIV.filename[kin$equivVidRow])
  }

  nVids <- nrow(kin)
  skippedFlow <- data.frame(vidName = character(0), section = character(0), reason = character(0))

  versions <- c("mouthFrontCorr", "mouthFront", "static", "staticRaw")

  # Output columns (one set per velocity version)
  flowCols <- c("U_flow_ef_meas_tmax", "U_flow_ef_meas_tmaxVel", "U_flow_ef_meas_acc_tmax", "U_flow_ef_meas_acc_tmaxVel",
                "U_flow_ff_tmax", "U_flow_ff_tmaxVel", "U_flow_ff_acc_tmax", "U_flow_ff_acc_tmaxVel",
                "U_flow_ef_meas_Max", "U_flow_ef_meas_acc_Max", "U_flow_ef_meas_mean", "U_flow_ff_mean", "U_flow_ef_meas_mean_bigGape", "U_flow_ff_mean_bigGape",
                "Cov_R_Uflowefmeas", "Cov_R_Uflowff", "Cov_A_Uflowff",
                "Total_ingested_volume", "Total_ingested_volume_direct", "U_flow_strike_frac")

  flowData <- data.frame(vidName = kin$vidName, PIV.filename = kin$PIV.filename, pivVideo = NA_character_,
                         Max_sub_pressure = NA_real_, pressure_peak_time = NA_real_,
                         Pressure_tmax = NA_real_, Pressure_tmaxVel = NA_real_,
                         pressure_node_zR = NA_real_, velocity_node_zR = NA_real_,
                         pressure_node_offR = NA_real_, pressure_decay_g = NA_real_, Max_sub_pressure_corr = NA_real_,
                         static_node_zR_tmax = NA_real_, static_node_zR_tmaxVel = NA_real_,
                         static_node_offR_tmax = NA_real_, static_node_offR_tmaxVel = NA_real_,
                         static_Muller_factor_tmax = NA_real_, static_Muller_factor_tmaxVel = NA_real_,
                         static_corrected_frac = NA_real_,
                         front_zR_tmax = NA_real_, front_zR_tmaxVel = NA_real_,
                         front_Muller_factor_tmax = NA_real_, front_Muller_factor_tmaxVel = NA_real_,
                         front_corrected_frac = NA_real_)
  for(colNow in c(flowCols, as.vector(outer(flowCols, paste0("_", versions), paste0)))) {
    flowData[[colNow]] <- NA_real_
  }

  # Within-video comparison of predicted and measured flow (fish frame), for each velocity version and buccal width model
    # predicted = V_buccal'(t) / A(t) from the kinematics (U_flow_ff_predicted*_t), measured = U_flow_ef_meas(t) + U_ram(t)
    # shapeCor: correlation over the strike frames (does the predicted flow rise and fall with the measured flow?)
    # shapeLag / shapeLagCor: the shift (frames, + = measured later) between -maxShapeLag and +maxShapeLag with the highest correlation
    # shapeRatio: median predicted / median measured (size; 1 = right on average)
    # vol*: the same four, for volume rates instead of speeds: predicted V_buccal'(t) (= predicted flow x A(t)) against
      # measured A(t) x U_flow_ff(t). No division by the gape area, so the small-gape frames early in the strike can't dominate
  maxShapeLag   <- 30   # frames (60 ms at 500 fps)
  widthSuffixes <- c(main = "", noRetract = "_noRetract", noWidth = "_noWidth")
  shapeMetrics  <- c("shapeCor", "shapeLag", "shapeLagCor", "shapeRatio", "volCor", "volLag", "volLagCor", "volRatio")
  for(colNow in as.vector(outer(outer(shapeMetrics, names(widthSuffixes), paste, sep = "_"), versions, paste, sep = "_"))) {
    flowData[[colNow]] <- NA_real_
  }

  # Time series: one row per video, one column per kinematic frame
  longestVid <- max(kin$nFrames, na.rm = TRUE)
  tsNames <- c("Pressure_t", "Muller_factor_t", "Muller_factor_front_t", paste0("U_flow_ef_meas_", versions, "_t"))
  flowTS <- lapply(setNames(tsNames, tsNames), function(x) matrix(NA_real_, nrow = nVids, ncol = longestVid, dimnames = list(kin$vidName, NULL)))


######################
#  Loop over videos  #
######################

  for(i in 1:nVids) {

    vidNow <- kin$vidName[i]

    ##### Which PIV video
      if(is.na(kin$PIV.filename[i]) || kin$PIV.filename[i] == "") {
        next   # no PIV for this video
      }
      # Included only if it's in the consistency-filtered results
      filtRowNow <- which(pivFiltered$videoKey == pivKey(kin$PIV.filename[i]))
      if(length(filtRowNow) == 0) {
        next   # not in the filtered PIV results, or no PIV for this video
      }
      if(length(filtRowNow) > 1) {
        skippedFlow <- addSkip(skippedFlow, vidNow, "PIV match", paste0("Matched ", length(filtRowNow), " rows in the filtered PIV results (needs 1)"))
        next
      }
      pivNow <- pivFiltered[filtRowNow, ]
      flowData$pivVideo[i] <- pivNow$video

      # Peak node positions (px) from the all-passed results. NA if the video isn't there (then no node distances or static correction)
      posRowNow <- which(pivResults$videoKey == pivKey(kin$PIV.filename[i]))
      if(length(posRowNow) == 1) {
        for(colNow in c("max_velocity_x_px", "max_velocity_y_px", "pressure_peak_x_px", "pressure_peak_y_px")) {
          pivNow[[colNow]] <- pivResults[[colNow]][posRowNow]
        }
      } else {
        for(colNow in c("max_velocity_x_px", "max_velocity_y_px", "pressure_peak_x_px", "pressure_peak_y_px")) {
          pivNow[[colNow]] <- NA_real_
        }
        skippedFlow <- addSkip(skippedFlow, vidNow, "PIV match", paste0("Matched ", length(posRowNow), " rows in the all-passed PIV results for the node positions (needs 1)"))
      }

      if(is.na(kin$nFrames[i]) || is.na(kin$frameRate[i])) {
        skippedFlow <- addSkip(skippedFlow, vidNow, "Kinematics", "No kinematic frames for this video")
        next
      }

      nNow       <- kin$nFrames[i]
      frameRate  <- kin$frameRate[i]
      timeNow    <- (0:(nNow - 1)) / frameRate   # s, same as the kinematics
      tmaxNow    <- kin$tmax[i]
      tmaxVelNow <- kin$tmaxVel[i]


    ##### Pressure (Pa, negative = sub-ambient)
      # The pipeline's peak value is the outcome. The trace is smoothed only for the values at tmax / tmaxVel
      flowData$Max_sub_pressure[i]   <- pivNow$pressure_peak_Pa
      flowData$pressure_peak_time[i] <- pivNow$pressure_peak_time_ms / 1000

      fileTraces <- file.path(pathRawTraces, paste0(pivNow$video, suffixPeakTraces))
      tracesNow  <- if(file.exists(fileTraces)) read.csv(fileTraces) else NULL
      if(is.null(tracesNow)) {
        skippedFlow <- addSkip(skippedFlow, vidNow, "Files", paste0("No peak node traces file: ", basename(fileTraces)))
      }

      if(!is.null(tracesNow)) {
        pressureNow <- toKinFrames(tracesNow$pressure_time_ms, tracesNow$pressure_Pa, nNow, frameRate)
        smP <- smoothTrace(pressureNow, timeNow, madFloorPressure, lambdaRefPressure)
        flowData$Pressure_tmax[i]    <- traceAt(smP, tmaxNow, timeNow, 0)
        flowData$Pressure_tmaxVel[i] <- traceAt(smP, tmaxVelNow, timeNow, 0)
        flowTS$Pressure_t[i, 1:nNow] <- smP$y
      }


    ##### PIV mouth landmarks: mouth centre, gape radius and strike axis in the PIV image (px)
      fileMarks <- file.path(pathMouthMarks, paste0(pivNow$video, suffixMouthMarks))
      marksOK   <- file.exists(fileMarks)
      if(!marksOK) {
        skippedFlow <- addSkip(skippedFlow, vidNow, "Files", paste0("No mouth landmarks file: ", basename(fileMarks)))
      }

      if(marksOK) {
        marksNow <- read.csv(fileMarks)
        kMarks   <- marksNow$frame_idx + 1   # PIV frame 0 = kinematic frame index 1
        keepM    <- kMarks >= 1 & kMarks <= nNow

        # Each coordinate: onto kinematic frames, Hampel, short gaps filled, smoothed
        markSm <- list()
        for(colNow in c("uj_x", "uj_y", "lj_x", "lj_y", "hinge_x", "hinge_y")) {
          coordNow <- rep(NA_real_, nNow)
          coordNow[kMarks[keepM]] <- marksNow[[colNow]][keepM]
          markSm[[colNow]] <- smoothTrace(coordNow, timeNow, madFloorLandmark, lambdaRefMarks)$y
        }

        mouthX <- (markSm$uj_x + markSm$lj_x) / 2
        mouthY <- (markSm$uj_y + markSm$lj_y) / 2
        gapeR  <- sqrt((markSm$lj_x - markSm$uj_x)^2 + (markSm$lj_y - markSm$uj_y)^2) / 2

        # Strike axis: perpendicular to the jaw-tip line, pointing away from the jaw joint (forward)
        jawUx <- (markSm$lj_x - markSm$uj_x) / (2 * gapeR)
        jawUy <- (markSm$lj_y - markSm$uj_y) / (2 * gapeR)
        axisX <- -jawUy
        axisY <-  jawUx
        flipNow <- ((mouthX - markSm$hinge_x) * axisX + (mouthY - markSm$hinge_y) * axisY) < 0
        axisX[flipNow %in% TRUE] <- -axisX[flipNow %in% TRUE]
        axisY[flipNow %in% TRUE] <- -axisY[flipNow %in% TRUE]

        # Distance of a fixed point from the mouth: along the axis (z, + = in front) and sideways (off), as multiples of R
        nodeDist <- function(px, py) {
          dx <- px - mouthX
          dy <- py - mouthY
          list(zR = (dx * axisX + dy * axisY) / gapeR, offR = abs(dx * jawUx + dy * jawUy) / gapeR)
        }

        # Checks: how far the pressure and velocity peak nodes were from the mouth at their peaks
        kPress <- round(pivNow$pressure_peak_time_ms / 1000 * frameRate) + 1
        kVel   <- round(pivNow$max_velocity_time_ms / 1000 * frameRate) + 1
        if(!is.na(kPress) && kPress >= 1 && kPress <= nNow) {
          dPress <- nodeDist(pivNow$pressure_peak_x_px, pivNow$pressure_peak_y_px)
          flowData$pressure_node_zR[i]   <- dPress$zR[kPress]
          flowData$pressure_node_offR[i] <- dPress$offR[kPress]

          # Peak pressure corrected to the aperture: p_node / g(s), only within the distance limits
          sNow <- dPress$zR[kPress]
          if(!is.na(sNow) && sNow >= 0 && sNow <= maxZoverR && dPress$offR[kPress] <= maxOffAxisOverR) {
            flowData$pressure_decay_g[i]      <- 1 - sNow / sqrt(1 + sNow^2)
            flowData$Max_sub_pressure_corr[i] <- flowData$Max_sub_pressure[i] / flowData$pressure_decay_g[i]
          }
        }
        if(!is.na(kVel) && kVel >= 1 && kVel <= nNow) {
          flowData$velocity_node_zR[i] <- nodeDist(pivNow$max_velocity_x_px, pivNow$max_velocity_y_px)$zR[kVel]
        }

        # Where the static velocity node was relative to the mouth at tmax and tmaxVel (multiples of the gape radius R),
          # and the Muller factor (1 + zR^2)^1.5 there, without the distance limits. For using the uncorrected flow with distance as a covariate
        dStat <- nodeDist(pivNow$max_velocity_x_px, pivNow$max_velocity_y_px)
        for(tNm in c("tmax", "tmaxVel")) {
          kNow <- round(kin[[tNm]][i] * frameRate) + 1
          if(!is.na(kNow) && kNow >= 1 && kNow <= nNow) {
            flowData[[paste0("static_node_zR_", tNm)]][i]        <- dStat$zR[kNow]
            flowData[[paste0("static_node_offR_", tNm)]][i]      <- dStat$offR[kNow]
            flowData[[paste0("static_Muller_factor_", tNm)]][i] <- (dStat$zR[kNow]^2 + 1)^1.5
          }
        }
      }


    ##### Velocity versions (mm/s, one value per kinematic frame)
      velRaw <- list()

      # Fixed node
      if(!is.null(tracesNow)) {
        velRaw$staticRaw <- toKinFrames(tracesNow$velocity_time_ms, tracesNow$velocity_m_per_s * 1000, nNow, frameRate)

        # Corrected to the aperture (needs the mouth landmarks)
        if(marksOK) {
          dNow <- nodeDist(pivNow$max_velocity_x_px, pivNow$max_velocity_y_px)
          mullerFactor <- (dNow$zR^2 + 1)^1.5          # (z^2 + R^2)^1.5 / R^3, with z in units of R
          useNow <- dNow$zR >= 0 & dNow$zR <= maxZoverR & dNow$offR <= maxOffAxisOverR
          mullerFactor[!(useNow %in% TRUE)] <- NA
          velRaw$static <- velRaw$staticRaw * mullerFactor
          flowTS$Muller_factor_t[i, 1:nNow] <- mullerFactor

          # Videos with no strike (no Tstart / Tend) have no strike frames
          if(!is.na(kin$startFrame[i]) && !is.na(kin$endFrame[i])) {
            flowData$static_corrected_frac[i] <- mean(!is.na(mullerFactor[kin$startFrame[i]:kin$endFrame[i]]))
          }
        }
      }

      # In front of the mouth (raw values; this script does its own smoothing)
      fileFront <- file.path(pathMouthFront, paste0(pivNow$video, suffixMouthFront))
      if(file.exists(fileFront)) {
        frontNow <- read.csv(fileFront)
        velRaw$mouthFront <- toKinFrames(frontNow$time_ms, frontNow$mouth_front_velocity_raw_m_per_s * 1000, nNow, frameRate)

        # Corrected to the aperture, frame by frame, from where the maximum was in that frame (needs the mouth landmarks)
        # Files without the position columns (not rerun): no corrected version for this video, logged
        frontPosOK <- all(c(frontXcol, frontYcol) %in% names(frontNow))
        if(!frontPosOK) {
          skippedFlow <- addSkip(skippedFlow, vidNow, "Files", paste0("Mouth-front file has no '", frontXcol, "' / '", frontYcol,
                                                                     "' columns: no mouthFrontCorr. File: ", basename(fileFront)))
        }
        if(marksOK && frontPosOK) {
          # Position of the maximum (px), cleaned and smoothed over time like the mouth landmarks
            # The maximum can jump between PIV nodes from frame to frame; jumps are removed by the Hampel filter and the spline
          frontX <- smoothTrace(toKinFrames(frontNow$time_ms, frontNow[[frontXcol]], nNow, frameRate), timeNow, madFloorLandmark, lambdaRefFrontPos)$y
          frontY <- smoothTrace(toKinFrames(frontNow$time_ms, frontNow[[frontYcol]], nNow, frameRate), timeNow, madFloorLandmark, lambdaRefFrontPos)$y
          dFront <- nodeDist(frontX, frontY)
          frontFactor <- (pmax(dFront$zR, 0)^2 + 1)^1.5   # 1 at or behind the aperture plane
          useFront <- dFront$zR <= maxZoverR & dFront$offR <= maxOffAxisOverR
          frontFactor[!(useFront %in% TRUE)] <- NA
          velRaw$mouthFrontCorr <- velRaw$mouthFront * frontFactor
          flowTS$Muller_factor_front_t[i, 1:nNow] <- frontFactor

          if(!is.na(kin$startFrame[i]) && !is.na(kin$endFrame[i])) {
            flowData$front_corrected_frac[i] <- mean(!is.na(frontFactor[kin$startFrame[i]:kin$endFrame[i]]))
          }
          # Where the maximum was at tmax and tmaxVel (multiples of R), and the factor there
          for(tNm in c("tmax", "tmaxVel")) {
            kNow <- round(kin[[tNm]][i] * frameRate) + 1
            if(!is.na(kNow) && kNow >= 1 && kNow <= nNow) {
              flowData[[paste0("front_zR_", tNm)]][i]            <- dFront$zR[kNow]
              flowData[[paste0("front_Muller_factor_", tNm)]][i] <- frontFactor[kNow]
            }
          }
        }
      } else {
        skippedFlow <- addSkip(skippedFlow, vidNow, "Files", paste0("No mouth front velocity file: ", basename(fileFront)))
      }


    ##### Kinematic series used with the flow (from the kinematics output)
      R_t     <- kinTS$R_t[1, 1:nNow, vidNow]       # gape radius (mm)
      A_t     <- kinTS$A_t[1, 1:nNow, vidNow]       # gape area (mm^2), NA without W_jaw
      U_ram_t <- kinTS$U_ram_t[1, 1:nNow, vidNow]   # ram speed along the strike axis (mm/s)

      strikeK  <- if(is.na(kin$startFrame[i]) || is.na(kin$endFrame[i])) integer(0) else kin$startFrame[i]:kin$endFrame[i]
      censored <- kin$Tend_censored[i] %in% TRUE   # no real strike end: no strike means


    ##### Outputs for each velocity version

      for(verNow in names(velRaw)) {

        smU <- smoothTrace(velRaw[[verNow]], timeNow, madFloorVel, lambdaRefVelocity)
        flowTS[[paste0("U_flow_ef_meas_", verNow, "_t")]][i, 1:nNow] <- smU$y

        resNow <- setNames(rep(NA_real_, length(flowCols)), flowCols)

        # Values at tmax and tmaxVel (mm/s, mm/s^2)
        resNow["U_flow_ef_meas_tmax"]        <- traceAt(smU, tmaxNow, timeNow, 0)
        resNow["U_flow_ef_meas_tmaxVel"]     <- traceAt(smU, tmaxVelNow, timeNow, 0)
        resNow["U_flow_ef_meas_acc_tmax"]    <- traceAt(smU, tmaxNow, timeNow, 1)
        resNow["U_flow_ef_meas_acc_tmaxVel"] <- traceAt(smU, tmaxVelNow, timeNow, 1)

        # Fish frame: add ram, since water moves toward the fish and the fish moves toward the water (identities with the SEM equations)
        resNow["U_flow_ff_tmax"]        <- resNow[["U_flow_ef_meas_tmax"]]        + kin$U_ram_tmax[i]
        resNow["U_flow_ff_tmaxVel"]     <- resNow[["U_flow_ef_meas_tmaxVel"]]     + kin$U_ram_tmaxVel[i]
        resNow["U_flow_ff_acc_tmax"]    <- resNow[["U_flow_ef_meas_acc_tmax"]]    + kin$U_ram_acc_tmax[i]
        resNow["U_flow_ff_acc_tmaxVel"] <- resNow[["U_flow_ef_meas_acc_tmaxVel"]] + kin$U_ram_acc_tmaxVel[i]

        # Over the strike (Tstart-Tend)
        if(length(strikeK) > 0) {
          uS <- smU$y[strikeK]
          resNow["U_flow_strike_frac"] <- mean(!is.na(uS))

          # Within-video shape comparison with the predicted flow (fish frame), for each width model
          measFF <- uS + U_ram_t[strikeK]
          for(wNm in names(widthSuffixes)) {
            predName <- paste0("U_flow_ff_predicted", widthSuffixes[[wNm]], "_t")
            if(!(predName %in% names(kinTS))) next
            predNow <- kinTS[[predName]][1, strikeK, vidNow]
            okNow <- is.finite(predNow) & is.finite(measFF)
            if(sum(okNow) < 10) next
            colEnd <- paste0("_", wNm, "_", verNow)
            flowData[[paste0("shapeCor", colEnd)]][i]   <- cor(predNow[okNow], measFF[okNow])
            flowData[[paste0("shapeRatio", colEnd)]][i] <- median(predNow[okNow]) / median(measFF[okNow])

            # Volume rates (mm^3/s): predicted V_buccal' = predicted flow x A, measured = measured fish-frame flow x A
            aStrike <- A_t[strikeK]
            predQ <- predNow * aStrike
            measQ <- measFF * aStrike
            okQ <- is.finite(predQ) & is.finite(measQ)
            if(sum(okQ) >= 10) {
              flowData[[paste0("volCor", colEnd)]][i]   <- cor(predQ[okQ], measQ[okQ])
              flowData[[paste0("volRatio", colEnd)]][i] <- median(predQ[okQ]) / median(measQ[okQ])
            }

            # Shifted correlations: measured frame k + lag against predicted frame k (speeds, then volume rates)
            lagNows <- -maxShapeLag:maxShapeLag
            lagCors <- matrix(NA_real_, length(lagNows), 2, dimnames = list(NULL, c("shape", "vol")))
            for(j in seq_along(lagNows)) {
              kMeas <- strikeK + lagNows[j]
              inNow <- kMeas >= 1 & kMeas <= nNow
              mNow  <- rep(NA_real_, length(strikeK))
              mNow[inNow] <- smU$y[kMeas[inNow]] + U_ram_t[kMeas[inNow]]
              qNow  <- rep(NA_real_, length(strikeK))
              qNow[inNow] <- mNow[inNow] * A_t[kMeas[inNow]]
              okLag <- is.finite(predNow) & is.finite(mNow)
              if(sum(okLag) >= 10) lagCors[j, "shape"] <- cor(predNow[okLag], mNow[okLag])
              okLagQ <- is.finite(predQ) & is.finite(qNow)
              if(sum(okLagQ) >= 10) lagCors[j, "vol"] <- cor(predQ[okLagQ], qNow[okLagQ])
            }
            for(kind in c("shape", "vol")) {
              if(any(is.finite(lagCors[, kind]))) {
                flowData[[paste0(kind, "Lag", colEnd)]][i]    <- lagNows[which.max(lagCors[, kind])]
                flowData[[paste0(kind, "LagCor", colEnd)]][i] <- max(lagCors[, kind], na.rm = TRUE)
              }
            }
          }
          if(any(!is.na(uS))) {
            resNow["U_flow_ef_meas_Max"] <- max(uS, na.rm = TRUE)
            # Largest acceleration over the strike (for the "peak" time basis in scripts 4 and 5)
            if(any(!is.na(smU$dy[strikeK]))) resNow["U_flow_ef_meas_acc_Max"] <- max(smU$dy[strikeK], na.rm = TRUE)
          }

          # Means, covariances and volume: only with enough frames, and a real strike end
            # Each uses the frames where all its inputs exist. Covariances are (1/n) * sum of products of deviations, as in the kinematics script
            # coverageRule: the minFramesBeforePeak frames up to tmaxVel, or a fraction of the strike (see the settings)
          if(!coverageRule %in% c("beforePeak", "fraction")) {
            stop("coverageRule must be 'beforePeak' or 'fraction'")
          }
          kPeak <- round(tmaxVelNow * frameRate) + 1
          coverageOK <- if(coverageRule == "fraction" || is.na(kPeak)) {   # no tmaxVel: fraction rule
            mean(!is.na(uS)) >= minStrikeFrac
          } else {
            kPeak - minFramesBeforePeak + 1 >= 1 && kPeak <= nNow &&
              all(!is.na(smU$y[(kPeak - minFramesBeforePeak + 1):kPeak]))
          }

          if(!censored && coverageOK) {

            uFF <- uS + U_ram_t[strikeK]   # fish frame (water toward the fish + fish toward the water)
            rS  <- R_t[strikeK]
            aS  <- A_t[strikeK]

            resNow["U_flow_ef_meas_mean"] <- mean(uS, na.rm = TRUE)

            okFF <- !is.na(uFF)
            resNow["U_flow_ff_mean"] <- mean(uFF[okFF])

            okR <- !is.na(uS) & !is.na(rS)
            resNow["Cov_R_Uflowefmeas"] <- mean((rS[okR] - mean(rS[okR])) * (uS[okR] - mean(uS[okR])))
            okRF <- okFF & !is.na(rS)
            resNow["Cov_R_Uflowff"] <- mean((rS[okRF] - mean(rS[okRF])) * (uFF[okRF] - mean(uFF[okRF])))

            okA <- okFF & !is.na(aS)
            if(sum(okA) >= 2) {
              covA  <- mean((aS[okA] - mean(aS[okA])) * (uFF[okA] - mean(uFF[okA])))
              resNow["Cov_A_Uflowff"] <- covA

              # Total ingested volume (mm^3) = T * <A> * <U_flow_ff> + T * Cov(A, U_flow_ff), T = Ttotal
              resNow["Total_ingested_volume"] <- kin$Ttotal[i] * mean(aS[okA]) * mean(uFF[okA]) + kin$Ttotal[i] * covA

              # Check: the integral of A * U_flow_ff over the strike (trapezoid rule). Close to the line above when few frames are missing
              tS <- timeNow[strikeK][okA]
              qS <- aS[okA] * uFF[okA]
              resNow["Total_ingested_volume_direct"] <- sum(diff(tS) * (head(qS, -1) + tail(qS, -1)) / 2)

              # Big-gape mean of the fish-frame flow, to compare with the predicted versions from the kinematics
              bigNow <- okA & aS >= uFlowAreaFrac * kin$A_Max[i]
              if(any(bigNow)) {
                resNow["U_flow_ff_mean_bigGape"]      <- mean(uFF[bigNow])
                resNow["U_flow_ef_meas_mean_bigGape"] <- mean(uS[bigNow])   # same frames, earth frame
              }
            }
          }
        }

        for(colNow in flowCols) {
          flowData[[paste0(colNow, "_", verNow)]][i] <- resNow[[colNow]]
          if(verNow == mainVersion) {
            flowData[[colNow]][i] <- resNow[[colNow]]
          }
        }
      }
  }


############
#  Checks  #
############

  print(paste0(sum(!is.na(flowData$pivVideo)), " of ", nVids, " kinematic videos matched to the consistency-filtered PIV results"))

  # Files missing, and other problems (counted)
  if(nrow(skippedFlow) > 0) {
    print(table(skippedFlow$section))
    print(sort(table(sub(":.*$", "", skippedFlow$reason)), decreasing = TRUE))
  }

  # How far the peak nodes were from the mouth at their peaks (multiples of the gape radius R)
    # Pressure is fine anywhere (the outcome is its own peak). For velocity, values above maxZoverR can't be corrected to the aperture at the peak
  print(summary(flowData$pressure_node_zR))
  print(paste0("Peak pressure corrected to the aperture for ", sum(!is.na(flowData$Max_sub_pressure_corr)), " of ",
               sum(!is.na(flowData$Max_sub_pressure)), " videos with pressure (node in front of the mouth, within maxZoverR and maxOffAxisOverR)"))
  print(summary(1 / flowData$pressure_decay_g))    # correction factor
  print(summary(flowData$Max_sub_pressure_corr))
  print(summary(flowData$velocity_node_zR))
  print(summary(flowData$static_corrected_frac))   # share of strike frames where the static correction could be applied
  print(summary(flowData$front_corrected_frac))    # same for the mouth-front correction (frames beyond maxZoverR / maxOffAxisOverR are NA)
  print(summary(flowData$front_zR_tmaxVel))        # where the mouth-front maximum was at tmaxVel (multiples of R; <= 0 = at or inside the aperture)
  print(summary(flowData$front_Muller_factor_tmaxVel))

  # The three velocity versions: size and agreement (mm/s)
  for(colNow in c("U_flow_ef_meas_tmaxVel", "U_flow_ef_meas_Max", "U_flow_ff_mean")) {
    print(paste0("--- ", colNow))
    vNow <- sapply(versions, function(v) flowData[[paste0(colNow, "_", v)]])
    print(summary(vNow))
    print(cor(vNow, use = "pairwise.complete.obs"))
  }

  # Which buccal width model predicts the measured flow best (big-gape means, fish frame)
    # Correlation (rank of fish) and median ratio (size, 1 = right on average), for each velocity version
  predNames <- c(main = "U_flow_ff_predicted_mean_bigGape", noRetract = "U_flow_ff_predicted_mean_bigGape_noRetract",
                 noWidth = "U_flow_ff_predicted_mean_bigGape_noWidth")
  for(v in versions) {
    measNow <- flowData[[paste0("U_flow_ff_mean_bigGape_", v)]]
    print(paste0("--- measured flow: ", v))
    print(rbind(cor   = sapply(predNames, function(p) cor(kin[[p]], measNow, use = "pairwise.complete.obs")),
                ratio = sapply(predNames, function(p) median(kin[[p]] / measNow, na.rm = TRUE))))
  }

  # Within-video shape: does the predicted flow (buccal volume) follow the measured flow through each strike?
    # Medians across videos, for each velocity version and width model. cor > 0.5 = the share of videos with a clear match
    # A best lag far from 0 means the peaks are offset in time (PIV vs. kinematics alignment, or the volume model's timing)
  for(v in versions) {
    for(kind in c("shape", "vol")) {
      print(paste0("--- within-video predicted vs. measured ", if(kind == "shape") "flow speed" else "volume rate (no division by gape area)",
                   " (fish frame): ", v))
      print(sapply(names(widthSuffixes), function(wNm) {
        colEnd <- paste0("_", wNm, "_", v)
        getNow <- function(m) flowData[[paste0(kind, m, colEnd)]]
        c(n_videos   = sum(!is.na(getNow("Cor"))),
          median_cor = median(getNow("Cor"), na.rm = TRUE),
          share_cor_over_0.5 = mean(getNow("Cor") > 0.5, na.rm = TRUE),
          median_best_lag = median(getNow("Lag"), na.rm = TRUE),
          lag_q25 = unname(quantile(getNow("Lag"), 0.25, na.rm = TRUE)),
          lag_q75 = unname(quantile(getNow("Lag"), 0.75, na.rm = TRUE)),
          share_lag_at_limit = mean(abs(getNow("Lag")) == maxShapeLag, na.rm = TRUE),
          median_best_lag_cor = median(getNow("LagCor"), na.rm = TRUE),
          median_ratio = median(getNow("Ratio"), na.rm = TRUE))
      }), digits = 3)
    }
  }

  # Total ingested volume: decomposition vs. direct integral, for each version (close to 1 when the trace covers the whole strike;
    # above 1 when frames are missing, since the decomposition uses the full duration T)
  for(v in versions) {
    print(paste0("--- ingested volume, decomposition / direct: ", v))
    print(summary(flowData[[paste0("Total_ingested_volume_", v)]] / flowData[[paste0("Total_ingested_volume_direct_", v)]]))
  }


##########
#  Save  #
##########

  write.csv(flowData, outFlow, row.names = FALSE)
  saveRDS(flowTS, outFlowTS)
  write.csv(skippedFlow, outFlowSkip, row.names = FALSE)
  print(paste0("Saved: ", outFlow))
