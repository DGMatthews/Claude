
#This script reads in SLEAP tracked landmarks from kinematic videos, processes them, and returns kinematic values for each video
  # Kinematic values to measure include:

      # Protrusion_ROM
      # Mandible_ROM
      # Hyoid_ROM
      # Neurocranium_ROM
      # Ttotal
      # Ttpg
      # tmax
      # Time_hyoid
      # Time_cranial
      # ⟨A⟩ - Mean gape area
      # ⟨U_flow_ff⟩
      # ⟨U_flow_ff_predicted(t)⟩
      # ⟨U_flow_ef_meas⟩
      # ⟨U_ram⟩
      # ⟨R(t)⟩
      # ⟨Protrusion(t)⟩
      # ⟨Mandible_depression(t)⟩
      # ⟨Hyoid_depression(t)⟩
      # ⟨Neurocranium_Adj_vel(t)⟩
      # ⟨Hyoid_Adj_vel(t)⟩
      # R(tmax)
      # R_vel (tmax)
      # U_flow_ef_meas(tmax)
      # U_flow_ef_meas_acc (tmax)
      # Protrusion(tmax)
      # Protrusion_vel(tmax)
      # Mandible_depression(tmax)
      # Mandible_dep_vel(tmax)
      # Mandible_angular_vel(tmax)
      # Mandible_ang(tmax)
      # Neurocranium_tip_disp(tmax)
      # Maxilla_tip_disp(tmax)
      # Maxilla_linear_vel(tmax)
      # Maxilla_ang(tmax)
      # Maxilla_angular_vel(tmax)
      # Neurocranium_linear_vel(tmax)
      # Neurocranium_rotation(tmax)
      # Neurocranium_angular_vel(tmax)
      # V_buccal_vel(tmax)
      # V_buccal_acc(tmax)
      # A(tmax)
      # A_vel(tmax)
      # Hyoid_depression(tmax)
      # Hyoid_vel(tmax)
      # R(t)
      # U_flow_ef_meas(t)
      # Protrusion(t)
      # hyoid_depression(t)
      # Neurocranium_linear_vel(t)
      # x_mouth(t)
      # Ceratohyal_rotation_angle(t)
      # Ceratohyal_lateral_angle_hypaxial(t)
      # Buccal_width_expansion(t)
      # Neurocranium_tip_disp(t)
      # Neurocranium_rotation(t)


  # Variables that are calculated or output based on the above

      # U_flow_ff(tmax)
      # U_flow_ff_acc (tmax)
      # U_ram(tmax)
      # U_ram(t)
      # U_ram_acc (tmax)
      # Relative_ttpg
      # Relative_tmax
      # Kinetic_Synchronization
      # Hyoid_vel(t)
      # Cov(R, U_ram)
      # Cov(R, U_flow_ff)
      # Cov(R, U_flow_ef_meas)




    
  #Note that all times are relative to the beginning of the video
    #But I can get a time that the event begins from the pressure measurements

#This is how you install the rhdf5 package
#install.packages("BiocManager")
#BiocManager::install("rhdf5")


require(rhdf5)
require(abind)
require(imputeTS)
require(readxl)

rm(list=ls())


# Read in all the data


#choose the folder with the digitized points that you want to analyze
path = 'N:/TRCxAC/Hybrid feeding project/HSV/KINE/SLEAP landmarks/F2'
#Remember to change the path to the correct output location.
setwd(path)
vidDataTemp <- read.csv('C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data/ACxTRCfeedingdata_Oct2026_FIXED.csv')
calTableTemp <- read.csv('N:/TRCxAC/Hybrid feeding project/HSV/Calibrations/CalValues.csv')
calValueTable <- calTableTemp[which(calTableTemp$type=="Kine"),]

Pressure_vel_Data <- read_excel("N:/TRCxAC/Hybrid feeding project/HSV/PIV/PIV Output/Pressure-velocity peak pipeline v2 (velocity-anchored)(0.4 pressure shape)/Final export/final_results_consistency_filtered.xlsx")

#only keep video data from F2 fish, and where the fish actually fed
vidData <- vidDataTemp[which(vidDataTemp$Generation=='F2'&vidDataTemp$Event!='-'),]

#Data from Sleap is stored in a .h5 file, in the HDF5 format
#This is a more complex data format that allows for relatively low memory load while working with large data sets
#The coordinates are stored in a 4th dimensional matrix called "tracks", which is of dimensions [a,b,x,n]
#a is the frame number
#b is the node number (names of the nodes are stored in the "node_names" data set)
#x holds the x,y coordinates of all the points
#n is the number of unique individuals that were tracked through the video (always 1 in this data set)

#Create a list of all the files in the folder
fileNames <- list.files(pattern="\\.h5$", all.files=TRUE, no..=TRUE)
fileNames2 <- gsub("\\.analysis\\.h5$", "",fileNames) #Keep the name of the file, but remove ".analysis.h5" from the end of it


#################
#   Test mode   #
#################
  # Runs the script on only a few videos, for tuning smoothing and debugging
  # NULL = run every video
  # Otherwise a vector of video numbers (positions in fileNames) or video names (fileNames2, without ".analysis.h5")
    # e.g. testVids <- c(1, 25, 300)   or   testVids <- c("2026-05-01_ACxTRC_A12.06_E03")
  testVids <- NULL   ### CHANGE

  # The pooled retraction ratio (beta) needs at least betaMinFrames eye-based frames, which a few test videos may not have
    # In test mode only, this value is used instead when there are too few frames
  testBeta <- 0.3    ### CHANGE  rough guess. Replace with the fitted value once the full dataset has been run

  if(!is.null(testVids)) {

    if(is.character(testVids)) {
      testKeep <- match(testVids, fileNames2)
      if(any(is.na(testKeep))) {
        stop("These test videos are not in the folder: ", paste(testVids[is.na(testKeep)], collapse = ", "))
      }
    } else {
      testKeep <- testVids
      if(any(testKeep < 1 | testKeep > length(fileNames))) {
        stop("Test video numbers must be between 1 and ", length(fileNames))
      }
    }

    fileNames  <- fileNames[testKeep]
    fileNames2 <- fileNames2[testKeep]

    print(paste0("TEST MODE: running ", length(fileNames), " videos"))
  }

#First view the structure of the file
print(h5ls(file=fileNames[1]))


#####################
#  Define functions #
#####################


  #This takes a linear model polynomial fit (saved as a list of coefficients) and turns the coefficients into a fully expressed polynomial function
              #This can then be used to find derivatives of the fit function
          makeLMintoFunction <- function (pc, expr = TRUE) {
            stringexpr <- paste("x", seq_along(pc) - 1, sep = " ^ ")
            stringexpr <- paste(stringexpr, pc, sep = " * ")
            stringexpr <- paste(stringexpr, collapse = " + ")
            if (expr) return(parse(text = stringexpr))
            else return(stringexpr)
          }



    #create a function that identifies outliers in tracking data
          applyHampel <- function(x, halfWindow=3, threshold=3, madFloor=1) { #
            
            flaggedFrames <- c()
            
            for(j in seq_along(x)){ #Loop through each frame
              
              if(is.na(x[j])==TRUE) {
                next
              }
              
              #Find the starting and ending frames of the current window
              minFrameNow <- max(1, j - halfWindow)
              maxFrameNow <- min(length(x), j + halfWindow)
              
              #calculate median and MAD
              pointsNow <- x[minFrameNow:maxFrameNow]
              medianNow <- median(pointsNow, na.rm=TRUE)
              
              MADnow <- mad(pointsNow, na.rm = TRUE)
              if(MADnow<madFloor) {MADnow=madFloor}
              
              
              # check if the frame is more that threshold*MAD from the median
              if(abs(x[j] - medianNow) > threshold * MADnow) {
                flaggedFrames <- c(flaggedFrames,j)
              }
              
            }
            
            return(flaggedFrames)
          }
        
        
        
    # Wrap an angle (radians) into -pi to pi, so differences between angles don't jump by 2*pi
      wrapAngle <- function(angle) {
          ((angle + pi) %% (2*pi)) - pi
        }
        
          # Four-bar output angle theta4 for a given input angle theta2 (radians)
            # Freudenstein equation, using the + solution, same as the SEM document
            # theta2 and theta4 are measured from the line running from the input pivot to the output pivot
        fourBarTheta4 <- function(theta2, inputLength, couplerLength, outputLength, fixedLength) {
          K1 <- fixedLength / inputLength
          K2 <- fixedLength / outputLength
          K3 <- (inputLength^2 - couplerLength^2 + outputLength^2 + fixedLength^2) / (2 * inputLength * outputLength)
          
          A <- K1 - cos(theta2)
          B <- -sin(theta2)
          C <- K2 * cos(theta2) - K3
          
          # Returns NaN if the linkage can't close at this input angle
          2 * atan((B + sqrt(A^2 + B^2 - C^2)) / (A + C))
        }
        
          # Four-bar coupler angle theta3 for a given input angle theta2 (radians)
        fourBarTheta3 <- function(theta2, inputLength, couplerLength, outputLength, fixedLength) {
          theta4 <- fourBarTheta4(theta2, inputLength, couplerLength, outputLength, fixedLength)
          atan2(outputLength * sin(theta4) - inputLength * sin(theta2),
                fixedLength + outputLength * cos(theta4) - inputLength * cos(theta2))
        }
        
          # Instantaneous KT: derivative of the output angle with respect to the input angle, at input angle theta2
            # output = "theta3" for the coupler (oral 4-bar: maxilla), "theta4" for the output link (opercular 4-bar: mandible)
        fourBarKTinst <- function(theta2, inputLength, couplerLength, outputLength, fixedLength, output) {
          theta4 <- fourBarTheta4(theta2, inputLength, couplerLength, outputLength, fixedLength)
          theta3 <- fourBarTheta3(theta2, inputLength, couplerLength, outputLength, fixedLength)
          
          if(output == "theta3") {
            return((inputLength / couplerLength) * sin(theta4 - theta2) / sin(theta3 - theta4))
          }
          if(output == "theta4") {
            return((inputLength / outputLength) * sin(theta2 - theta3) / sin(theta4 - theta3))
          }
          stop("output must be 'theta3' or 'theta4'")
        }
        
        
    # Correlation for the checks. Returns NA (instead of an error that would stop the script) if fewer than 3 videos have both values
    corCheck <- function(x, y) {
      bothNow <- !is.na(x) & !is.na(y)
      if(sum(bothNow) < 3) {
        return(NA_real_)
      }
      cor(x[bothNow], y[bothNow])
    }


    # Add a row to the skipped-analysis log
      # skippedLog = the current log (skippedVids), returns the log with the new row added
      # section = which part of the analysis was skipped, reason = why
    addSkip <- function(skippedLog, vidName, section, reason) {
      rbind(skippedLog, data.frame(vidName = vidName, section = section, reason = reason))
    }
    
    
    
    
    
    
    # Hyoid outputs from a raw depth (or distance) series and its resting zero
      # Used three times: skull-frame depth and nasal-distance excursion inside the loop, and triangle depth after the loop
      # raw_t = depth or distance in each frame (mm), rest = resting value from the scans (mm, NA if no scan)
      # The spline is fit to the raw series, so velocities and Time don't depend on the zero
      # Values that depend on the zero are left NA if depression goes below -negTolerance in any visible frame (the zero must be wrong)
      # Returns a list:
        # values: Time, visible_frac, vel_mean_visible, vel_tmax, Max, tmax, t_mean, vel_mean (NA if not calculated)
        # vel_t, depression_t: time series (NULL if not calculated)
        # log: reasons to add to skippedVids
    hyoidOutputs <- function(raw_t, rest, timeSeriesNow, startPos, endPos, tmax, Tstart,
                             lambdaRefDerived, spanRef, negTolerance, missingMsg) {

      out <- list(values = c(Time = NA_real_, visible_frac = NA_real_, vel_mean_visible = NA_real_, vel_tmax = NA_real_,
                             Max = NA_real_, tmax = NA_real_, t_mean = NA_real_, vel_mean = NA_real_),
                  vel_t = NULL, depression_t = NULL, log = character())

      goodFrames   <- which(!is.na(raw_t))
      strikeFrames <- intersect(startPos:endPos, goodFrames)

      if(length(goodFrames) == 0) {
        out$log <- missingMsg
        return(out)
      }
      if(length(strikeFrames) == 0) {
        out$log <- "Hyoid not visible between Tstart and Tend"
        return(out)
      }

      # Time of peak (s from Tstart), from the raw series so it doesn't depend on the zero
      peakFrame <- strikeFrames[which.max(raw_t[strikeFrames])]
      out$values["Time"] <- timeSeriesNow[peakFrame] - Tstart

      # Fraction of Tstart-Tend frames where the hyoid was visible
      out$values["visible_frac"] <- length(strikeFrames) / length(startPos:endPos)

      rawTmax <- NA_real_

      if(length(goodFrames) >= 10) {

        # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
        spanNow   <- max(goodFrames) - min(goodFrames)
        lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

        rawSmooth <- smooth.spline(x = timeSeriesNow[goodFrames], y = raw_t[goodFrames], lambda = lambdaNow)

        # Velocity (mm/s). Same as d(depression)/dt, because the zero is a constant
        vel_t <- rep(NA_real_, length(raw_t))
        vel_t[goodFrames] <- predict(rawSmooth, x = timeSeriesNow[goodFrames], deriv = 1)$y
        out$vel_t <- vel_t

        # Plain mean velocity over visible frames from Tstart to the peak
        velFrames <- intersect(startPos:peakFrame, goodFrames)
        out$values["vel_mean_visible"] <- mean(vel_t[velFrames])

        if(!is.na(tmax)) {
          if(tmax >= timeSeriesNow[min(goodFrames)] && tmax <= timeSeriesNow[max(goodFrames)]) {
            rawTmax <- predict(rawSmooth, x = tmax, deriv = 0)$y
            out$values["vel_tmax"] <- predict(rawSmooth, x = tmax, deriv = 1)$y
          } else {
            out$log <- c(out$log, "Values at tmax not calculated: tmax outside the frames where the hyoid is visible")
          }
        }
      } else {
        out$log <- c(out$log, "Fewer than 10 visible frames: no velocities or values at tmax")
      }

      # Everything below needs the scan zero. Missing scan values aren't logged, since they aren't a landmarking problem
      if(is.na(rest)) {
        return(out)
      }

      depression_t <- raw_t - rest
      out$depression_t <- depression_t   # saved even if negative, so the problem can be inspected

      if(any(depression_t[goodFrames] < -negTolerance)) {
        out$log <- c(out$log, paste0("FLAG: negative depression (min ", round(min(depression_t[goodFrames]), 3), " mm). Depression values not saved"))
        return(out)
      }

      # Largest depression during the strike (mm)
      out$values["Max"] <- depression_t[peakFrame]

      # Depression at tmax (mm). Raw spline value minus the zero, identical to fitting the spline to depression. NA if tmax wasn't calculated
      out$values["tmax"] <- rawTmax - rest

      # Mean depression over Tstart-Tend (mm). 0 at Tstart/Tend where not visible, gaps filled linearly
      depressionStrike <- depression_t[startPos:endPos]
      if(is.na(depressionStrike[1])) {
        depressionStrike[1] <- 0
      }
      if(is.na(depressionStrike[length(depressionStrike)])) {
        depressionStrike[length(depressionStrike)] <- 0
      }
      depressionStrike <- na_interpolation(depressionStrike, option = "linear")
      out$values["t_mean"] <- mean(depressionStrike)

      # Mean velocity from Tstart to the peak (mm/s) = Max / Time, with depression 0 at Tstart
      if(out$values["Time"] > 0) {
        out$values["vel_mean"] <- out$values["Max"] / out$values["Time"]
      } else {
        out$log <- c(out$log, "FLAG: peak depression is at Tstart, so Hyoid_vel_mean was not calculated")
      }

      out
    }


    # One QC plot panel: a time series with Tstart/Tend (grey dotted) and tmax (red)
      # ySecond (optional) is drawn dashed, for a second version of the same variable
    qcPanel <- function(timeNow, yMain, ySecond, main, ylab, Tstart, Tend, tmax) {

      if(all(is.na(yMain)) && (is.null(ySecond) || all(is.na(ySecond)))) {
        plot.new()
        title(main = paste0(main, " (no data)"))
        return(invisible(NULL))
      }

      yRange <- range(c(yMain, ySecond), na.rm = TRUE)
      plot(timeNow, yMain, type = "l", lwd = 2, ylim = yRange, xlab = "Time (s)", ylab = ylab, main = main)

      if(!is.null(ySecond)) {
        lines(timeNow, ySecond, lty = 2, col = "grey40")
      }

      abline(h = 0, col = "grey85")
      strikeLines <- c(Tstart, Tend)
      abline(v = strikeLines[!is.na(strikeLines)], lty = 3, col = "grey50")
      if(!is.na(tmax)) {
        abline(v = tmax, col = "red")
      }
    }



##################
#  Final prep

nVids <- length(fileNames)

#Predefine dataset to hold all my data
finalData <- data.frame(vidName = fileNames2,
                        Animal.ID = rep(NA_character_, nVids),
                        Event = rep(NA_character_, nVids),
                        eventNum = rep(NA_real_, nVids),
                        family = rep(NA_character_, nVids),
                        Generation = rep("F2", nVids),
                        date = rep(NA_character_, nVids),
                        equivVidRow = rep(NA_real_, nVids),
                        daysNoFood = rep(NA_real_, nVids),
                        cal_lookup = rep(NA_character_, nVids),
                        preyNumber = rep(NA_character_, nVids),
                        frameRate = rep(NA_real_, nVids),
                        cal_pix.mm = rep(NA_real_, nVids),
                        hampelFlags = rep(NA_real_, nVids),
                        Protrusion_ROM = rep(NA_real_, nVids),
                        Mandible_ROM = rep(NA_real_, nVids),
                        Hyoid_ROM = rep(NA_real_, nVids),
                        Neurocranium_ROM = rep(NA_real_, nVids),
                        Protrusion_Max = rep(NA_real_, nVids),
                        Mandible_Max = rep(NA_real_, nVids),
                        Hyoid_Max = rep(NA_real_, nVids),
                        Neurocranium_Max = rep(NA_real_, nVids),
                        Ttotal = rep(NA_real_, nVids),
                        Tstart = rep(NA_real_, nVids),
                        Tend = rep(NA_real_, nVids),
                        Ttpg = rep(NA_real_, nVids),
                        tmax = rep(NA_real_, nVids),
                        tmaxVel = rep(NA_real_, nVids),
                        Time_hyoid = rep(NA_real_, nVids),
                        Time_cranial = rep(NA_real_, nVids),
                        U_flow_ff_mean = rep(NA_real_, nVids),
                        U_ram_mean = rep(NA_real_, nVids),
                        R_t_mean = rep(NA_real_, nVids),
                        Protrusion_t_mean = rep(NA_real_, nVids),
                        Mandible_depression_t_mean = rep(NA_real_, nVids),
                        Hyoid_depression_t_mean = rep(NA_real_, nVids),
                        R_tmax = rep(NA_real_, nVids),
                        R_vel_tmax = rep(NA_real_, nVids),
                        U_flow_ef_meas_tmax = rep(NA_real_, nVids),
                        U_flow_ef_meas_acc_tmax = rep(NA_real_, nVids),
                        Protrusion_tmax = rep(NA_real_, nVids),
                        Protrusion_vel_tmax = rep(NA_real_, nVids),
                        Mandible_depression_tmax = rep(NA_real_, nVids),
                        Mandible_dep_vel_tmax = rep(NA_real_, nVids),
                        Mandible_angular_vel_tmax = rep(NA_real_, nVids),
                        Mandible_ang_tmax = rep(NA_real_, nVids),
                        Mandible_tip_disp_tmax = rep(NA_real_, nVids),
                        Mandible_tip_disp_vel_tmax = rep(NA_real_, nVids),
                        Mandible_tip_disp_Max = rep(NA_real_, nVids),
                        Mandible_depression_Max = rep(NA_real_, nVids),
                        Neurocranium_tip_disp_tmax = rep(NA_real_, nVids),
                        Neurocranium_angular_vel_tmax = rep(NA_real_, nVids),
                        Maxilla_tip_disp_tmax = rep(NA_real_, nVids),
                        Maxilla_linear_vel_tmax = rep(NA_real_, nVids),
                        Maxilla_ang_tmax = rep(NA_real_, nVids),
                        Maxilla_angular_vel_tmax = rep(NA_real_, nVids),
                        Neurocranium_linear_vel_tmax = rep(NA_real_, nVids),
                        Neurocranium_rotation_tmax = rep(NA_real_, nVids),
                        Hyoid_depression_tmax = rep(NA_real_, nVids),
                        Hyoid_vel_tmax = rep(NA_real_, nVids),
                        U_flow_ff_tmax = rep(NA_real_, nVids),
                        U_flow_ff_acc_tmax = rep(NA_real_, nVids),
                        U_ram_tmax = rep(NA_real_, nVids),
                        U_ram_acc_tmax = rep(NA_real_, nVids),
                        Relative_ttpg = rep(NA_real_, nVids),
                        Relative_tmax = rep(NA_real_, nVids),
                        Kinetic_Synchronization = rep(NA_real_, nVids),
                        Mandible_tip_dist_rest = rep(NA_real_, nVids),
                        Oral_theta2_tmax = rep(NA_real_, nVids),
                        Oral_theta3_tmax = rep(NA_real_, nVids),
                        Oral_theta4_tmax = rep(NA_real_, nVids),
                        Oral_KT_tmax = rep(NA_real_, nVids),
                        Opercular_theta2_tmax = rep(NA_real_, nVids),
                        Opercular_theta3_tmax = rep(NA_real_, nVids),
                        Opercular_theta4_tmax = rep(NA_real_, nVids),
                        Opercular_KT_tmax = rep(NA_real_, nVids),
                        Maxilla_ang_predicted_tmax = rep(NA_real_, nVids),
                        Mandible_ang_predicted_tmax = rep(NA_real_, nVids),
                        strikeDirection = rep(NA_real_, nVids),
                        Neurocranium_tip_disp_eye_tmax = rep(NA_real_, nVids),
                        Neurocranium_angular_vel_mean = rep(NA_real_, nVids),
                        Neurocranium_linear_vel_mean = rep(NA_real_, nVids),
                        Neurocranium_eyeNasal_stretch_max = rep(NA_real_, nVids),
                        Neurocranium_eyeQC_flag = rep(NA, nVids),
                        Mandible_tip_disp_ROM = rep(NA_real_, nVids),
                        Mandible_depression_ROM = rep(NA_real_, nVids),
                        Maxilla_Max = rep(NA_real_, nVids),
                        Maxilla_cam_Max = rep(NA_real_, nVids),
                        Maxilla_ang_cam_tmax = rep(NA_real_, nVids),
                        Maxilla_angular_vel_cam_tmax = rep(NA_real_, nVids),
                        Maxilla_tip_disp_cam_tmax = rep(NA_real_, nVids),
                        Maxilla_linear_vel_cam_tmax = rep(NA_real_, nVids),
                        Hyoid_vel_mean = rep(NA_real_, nVids),
                        Hyoid_vel_mean_visible = rep(NA_real_, nVids),
                        Hyoid_visible_frac = rep(NA_real_, nVids),
                        Hyoid_excursion_tmax_backup = rep(NA_real_, nVids),
                        Hyoid_excursion_vel_tmax_backup = rep(NA_real_, nVids),
                        Hyoid_excursion_t_mean_backup = rep(NA_real_, nVids),
                        Hyoid_excursion_Max_backup = rep(NA_real_, nVids),
                        Hyoid_excursion_ROM_backup = rep(NA_real_, nVids),
                        Time_hyoid_backup = rep(NA_real_, nVids),
                        Hyoid_excursion_vel_mean_backup = rep(NA_real_, nVids),
                        Hyoid_excursion_vel_mean_visible_backup = rep(NA_real_, nVids),
                        Hyoid_visible_frac_backup = rep(NA_real_, nVids),
                        nFrames = rep(NA_real_, nVids),
                        startFrame = rep(NA_real_, nVids),
                        endFrame = rep(NA_real_, nVids),
                        Hyoid_depth_source = rep(NA_character_, nVids),
                        V_buccal_vel_tmax = rep(NA_real_, nVids),
                        V_buccal_acc_tmax = rep(NA_real_, nVids),
                        V_buccal_Max = rep(NA_real_, nVids),
                        Buccal_width_expansion_Max = rep(NA_real_, nVids),
                        U_flow_ff_predicted_mean = rep(NA_real_, nVids),
                        U_flow_ff_predicted_mean_bigGape = rep(NA_real_, nVids),
                        Width_hypaxial_floor_frames = rep(NA_real_, nVids),
                        Width_total_floor_frames = rep(NA_real_, nVids),
                        V_buccal_neuroIncluded = rep(NA, nVids),
                        U_ram_mean_frames = rep(NA_real_, nVids),
                        U_ram_body_tmax = rep(NA_real_, nVids),
                        U_ram_body_acc_tmax = rep(NA_real_, nVids),
                        U_ram_body_mean = rep(NA_real_, nVids),
                        U_ram_body_mean_frames = rep(NA_real_, nVids),
                        strikeAxisTilt = rep(NA_real_, nVids),
                        strikeAxisFromNasal = rep(NA, nVids),
                        Kinetic_Synchronization_backup = rep(NA_real_, nVids),
                        Lag_hyoid_gape = rep(NA_real_, nVids),
                        Lag_cranial_gape = rep(NA_real_, nVids),
                        A_tmax = rep(NA_real_, nVids),
                        A_vel_tmax = rep(NA_real_, nVids),
                        A_mean = rep(NA_real_, nVids),
                        A_Max = rep(NA_real_, nVids),
                        Cov_R_Uram = rep(NA_real_, nVids),
                        Cov_R_Uflowff = rep(NA_real_, nVids),
                        Cov_R_Uflowefmeas = rep(NA_real_, nVids)
                        )





#Find the longest video and record how many frames it has
longestVid <- 0
for (i in 1:nVids) {
  info <- h5ls(fileNames[i])
  dims <- as.numeric(strsplit(info$dim[info$name == "tracks"], " x ")[[1]])
  checkLength <- dims[1]
  if (checkLength>longestVid) {longestVid=checkLength}
}


#Create a list of arrays to hold time series data
        # Specifically these variables
                # "R_t",
                # "U_flow_ef_meas_t",
                # "Protrusion_t",
                # "hyoid_depression_t",
                # "Neurocranium_linear_vel_t",
                # "x_mouth_t",
                # "Neurocranium_tip_disp_t",
                # U_ram_t
                # Mandible_depression_t
                # Hyoid_vel_t
                # Mandible_ang_t
                # Maxilla_ang_t
                # Maxilla_tip_disp_t

#Make a list of arrays to hold all the time series. Each array is a time series that has one value x the number of frames x the number of events
finalTimeSeriesData <- list(R_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("R_t", NULL, fileNames2)),
                            U_flow_ef_meas_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("U_flow_ef_meas_t", NULL, fileNames2)),
                            Protrusion_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Protrusion_t", NULL, fileNames2)),
                            Hyoid_depression_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Hyoid_depression_t", NULL, fileNames2)),
                            Neurocranium_linear_vel_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Neurocranium_linear_vel_t", NULL, fileNames2)),
                            x_mouth_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("x_mouth_t", NULL, fileNames2)),
                            Neurocranium_tip_disp_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Neurocranium_tip_disp_t", NULL, fileNames2)),
                            U_ram_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("U_ram_t", NULL, fileNames2)),
                            Mandible_depression_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Mandible_depression_t", NULL, fileNames2)),
                            Mandible_tip_disp_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Mandible_tip_disp_t", NULL, fileNames2)),
                            Hyoid_vel_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Hyoid_vel_t", NULL, fileNames2)),
                            Mandible_ang_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Mandible_ang_t", NULL, fileNames2)),
                            Maxilla_ang_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Maxilla_ang_t", NULL, fileNames2)),
                            Maxilla_tip_disp_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Maxilla_tip_disp_t", NULL, fileNames2)),
                            Neurocranium_rotation_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Neurocranium_rotation_t", NULL, fileNames2)),
                            Neurocranium_angular_vel_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Neurocranium_angular_vel_t", NULL, fileNames2)),
                            Neurocranium_tip_disp_eye_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Neurocranium_tip_disp_eye_t", NULL, fileNames2)),
                            Maxilla_ang_cam_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Maxilla_ang_cam_t", NULL, fileNames2)),
                            Maxilla_tip_disp_cam_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Maxilla_tip_disp_cam_t", NULL, fileNames2)),
                            Hyoid_excursion_t_backup = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Hyoid_excursion_t_backup", NULL, fileNames2)),
                            Hyoid_excursion_vel_t_backup = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Hyoid_excursion_vel_t_backup", NULL, fileNames2)),
                            Hyoid_retraction_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Hyoid_retraction_t", NULL, fileNames2)),
                            Hyoid_nasal_dist_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Hyoid_nasal_dist_t", NULL, fileNames2)),
                            Ceratohyal_rotation_angle_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Ceratohyal_rotation_angle_t", NULL, fileNames2)),
                            Ceratohyal_lateral_angle_hypaxial_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Ceratohyal_lateral_angle_hypaxial_t", NULL, fileNames2)),
                            Ceratohyal_lateral_angle_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Ceratohyal_lateral_angle_t", NULL, fileNames2)),
                            Ceratohyal_half_width_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Ceratohyal_half_width_t", NULL, fileNames2)),
                            Buccal_width_expansion_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Buccal_width_expansion_t", NULL, fileNames2)),
                            V_buccal_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("V_buccal_t", NULL, fileNames2)),
                            V_buccal_vel_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("V_buccal_vel_t", NULL, fileNames2)),
                            U_flow_ff_predicted_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("U_flow_ff_predicted_t", NULL, fileNames2)),
                            U_ram_body_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("U_ram_body_t", NULL, fileNames2)),
                            A_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("A_t", NULL, fileNames2))
                            )




#Make a log of analyses that were skipped (or flagged) for each video, so I can check the landmarking quality
  # One row per skipped section per video
skippedVids <- data.frame(vidName = character(), section = character(), reason = character())







#Now open up each data set to record info about the strike

      # dataNow <- h5read(file=fileNames[i], #Read in the ith data set
      #                     name='/tracks',
      #                     index= list(NULL,NULL,NULL,NULL)) #This is here so I remember how to subset from h5read. Change any "NULL" to a number to get just that element number.
      #  #assign(fileNames2[i], dataNow) #This allows you to assign a value to a changing variable name

for(i in 1:nVids) {
  
  
  #The problem is that not all videos have the same naming format
    #Solution is to find where ACxTRC is in each name
    #Then separate between ones that have F2 before the animal name and those that don't
  
  splitStr <- strsplit(fileNames2[i],'_')
  
  #find where ACxTRC is
  stablePos <- which(splitStr[[1]]=='ACxTRC')
  
  #Extract data from filename
  if (splitStr[[1]][stablePos+2]=='F2'){
   finalData$date[i] <- splitStr[[1]][stablePos+1]
   finalData$Animal.ID[i] <- splitStr[[1]][stablePos+3]
   finalData$Event[i] <- splitStr[[1]][stablePos+4]
   finalData$eventNum[i] <- as.numeric(substring(finalData$Event[i],2))
  } else{
    finalData$date[i] <- splitStr[[1]][stablePos+1]
    finalData$Animal.ID[i] <- splitStr[[1]][stablePos+2]
    finalData$Event[i] <- splitStr[[1]][stablePos+3]
    finalData$eventNum[i] <- as.numeric(substring(finalData$Event[i],2))
  }
  
}
  
 

#clean up the spreadsheets so they can be combined

#make animal names consistent by adding an 'A' to the start of any names that don't have it (keep the A to make sure it isn't treated numerically)
fixRow <- which(substring(finalData$Animal.ID,1,1)!='A')
finalData$Animal.ID[fixRow] <- paste0('A',finalData$Animal.ID[fixRow])
#Also replace "+" and "-" in animal names with a period
finalData$Animal.ID <- gsub('\\+','\\.',finalData$Animal.ID)
finalData$Animal.ID <- gsub('\\-','\\.',finalData$Animal.ID)

#Do it all again for video data
vidData$Animal.ID <- substring(vidData$Animal.ID,4)
fixRow2 <- which(substring(vidData$Animal.ID,1,1)!='A')
vidData$Animal.ID[fixRow2] <- paste0('A',vidData$Animal.ID[fixRow2])

vidData$Animal.ID <- gsub('\\+','\\.',vidData$Animal.ID)
vidData$Animal.ID <- gsub('\\-','\\.',vidData$Animal.ID)


#Now find the matching row in video data and fill in final data frame with the info
#Important to include the date to differentiate between fish that have the same individual ID but come from different parents
#e.g.F1.2 A20.01 vs. F1.3 A20.01

print(head(vidData))
print(head(finalData))

for(sigh in 1:nrow(vidData)) {
  vidData$Date2[sigh] <- strsplit(vidData$Kine.filename,'_')[[sigh]][3]
}


for(check in 1:nrow(finalData)){
  dateNow <- finalData$date[check]
  animalNow <- finalData$Animal.ID[check]
  eventNow <- finalData$Event[check]
  
  vidRowNow <- which(vidData$Date2==dateNow & vidData$Animal.ID==animalNow & vidData$Event==eventNow)
  
  if(length(vidRowNow) != 1) {
    stop("Problem matching ", fileNames2[check], " to a single row from the video data sheet. It matched ", length(vidRowNow)," rows: ", paste(vidRowNow, collapse = ", "))
  }
  
  finalData$equivVidRow[check] <- vidRowNow
  
  finalData$daysNoFood[check] <- as.numeric(vidData$Days.without.food[vidRowNow])
  finalData$frameRate[check] <- vidData$Kine.frame.rate[vidRowNow]
  finalData$family[check] <- vidData$Family[vidRowNow]
  finalData$cal_lookup[check] <- substring(vidData$Kine.calibration[vidRowNow],1,nchar(vidData$Kine.calibration[vidRowNow])-7)
  finalData$preyNumber[check] <- vidData$Prey.number[vidRowNow]
  
  #Get the kinematic video calibration value
  calRowNow <- which(calValueTable$cal.File==finalData$cal_lookup[check])
  if(length(calRowNow)!=1){
    stop("Calibration lookup failed for ", fileNames2[check], ": '", finalData$cal_lookup[check], "' matched ", length(calRowNow), " rows")
  }
  finalData$cal_pix.mm[check] <- calValueTable$pix.mm[calRowNow]
  #Don't forget to correct to mm and seconds
}


  
#Make sure there are no duplicates of the same events in the spreadsheet
finalData$Unique.ID <- paste(finalData$family,finalData$Generation,finalData$Animal.ID,sep='_')
finalData$Event.ID <- paste(finalData$date,finalData$family,finalData$Generation,finalData$Animal.ID,finalData$Event,sep='_')

print(finalData$Event.ID[which(duplicated(finalData$Event.ID)==TRUE)]) #No duplicates, onto the analysis

print(length(unique(finalData$Unique.ID)))




# Get time of peak pressure and peak velocity from the PIV analysis for each video
  # Remove the suffix from the PIV video names, because it's not in the data sheet
  Pressure_vel_Data$video <- sub("_0{5,6}$", "", Pressure_vel_Data$video)
  
  # Columns to hold the PIV values (NA for videos without good PIV data)
  finalData$tPeakPressure <- NA_real_   # seconds
  finalData$tPeakVelocity <- NA_real_   # seconds
  finalData$peakPressure_Pa <- NA_real_
  finalData$peakVelocity_m.s <- NA_real_
  
  for(i in 1:nrow(finalData)){
    
    # PIV file name for this video, from its matching row in the data sheet
    pivNameNow <- vidData$PIV.filename[finalData$equivVidRow[i]]
    
    # Find that PIV file in the pressure/velocity results
    pressureRowNow <- which(Pressure_vel_Data$video == pivNameNow)
    
    # No PIV data for this video: leave the PIV columns as NA
    if(length(pressureRowNow) == 0) {
      next
    }
    
    # More than one match means something is wrong with the data
    if(length(pressureRowNow) > 1) {
      stop("PIV lookup for ", fileNames2[i], " matched ", length(pressureRowNow), " rows in the pressure/velocity data")
    }
    
    # Convert times from ms to seconds so they match timeSeriesNow
    finalData$tPeakPressure[i]    <- Pressure_vel_Data$pressure_peak_time_ms[pressureRowNow] / 1000
    finalData$tPeakVelocity[i]    <- Pressure_vel_Data$max_velocity_time_ms[pressureRowNow] / 1000
    finalData$peakPressure_Pa[i]  <- Pressure_vel_Data$pressure_peak_Pa[pressureRowNow]
    finalData$peakVelocity_m.s[i] <- Pressure_vel_Data$max_velocity_m_per_s[pressureRowNow]
  }
  
  # Get tmax and tmaxVel
  finalData$tmax <- finalData$tPeakPressure
  finalData$tmaxVel <- finalData$tPeakVelocity
  
  # Check: how many videos got PIV data, and were any PIV results not used?
  print(sum(!is.na(finalData$tmax)))
  print(setdiff(Pressure_vel_Data$video, vidData$PIV.filename[finalData$equivVidRow]))
  
  
  
  
  
###########################
# Anatomical measurements #
###########################

  # Read in the measurements from the scans for each individual
    # All lengths in mm, measured in a lateral (sagittal plane) view so they match what the camera sees
    scanData <- read.csv('PATH/TO/scan_measurements.csv')        ### CHANGE path
    
    # Columns in the scan file that identify each individual
    scanFamilyCol <- "Family"          ### CHANGE to the family column name in the scan file
    scanAnimalCol <- "Animal.ID"       ### CHANGE to the animal ID column name in the scan file
    
    # Measurements to pull from the scans
      # Left side = name it will have in finalData
      # Right side = column name in the scan file
      scanColumns <- c(
        # Mandible
        Mandible_length            = "Mandible_length",          ### CHANGE  jaw joint to LJ tip, anteroposterior leg
        Nasal_joint_length         = "Nasal_joint_length",       ### CHANGE  nasal tip to jaw joint, lateral view
        # Gape, buccal volume, neurocranium, hyoid
        W_jaw                      = "W_jaw",                    ### CHANGE
        W_head                     = "W_head",                   ### CHANGE
        D_head                     = "D_head",                   ### CHANGE
        Neurocranium_length        = "Neurocranium_length",      ### CHANGE
        Ceratohyal_length          = "Ceratohyal_length",        ### CHANGE
        D_head_backup              = "D_head_backup",            ### CHANGE  nasal tip to resting hyoid (urohyal) point, lateral view
        Nasal_eye_AP_length        = "Nasal_eye_AP_length",      ### CHANGE  AP distance, nasal tip to eye-level landmark (caudal base of parasphenoid)
        Hyoid_rest_AP              = "Hyoid_rest_AP",            ### CHANGE  AP distance, eye-level landmark to rostral ceratohyal tip. Positive if the tip is posterior
        Ceratohyal_AP_rest         = "Ceratohyal_AP_rest",       ### CHANGE  AP distance between rostral and caudal ceratohyal tips at rest
        Ceratohyal_DV_rest         = "Ceratohyal_DV_rest",       ### CHANGE  DV distance between rostral and caudal ceratohyal tips at rest. Positive if the rostral tip is ventral
        # Oral four-bar
        Mandible_input_length      = "Mandible_input_length",    ### CHANGE  input link
        Maxilla_length             = "Maxilla_length",           ### CHANGE  coupler link
        Nasal_length               = "Nasal_length",             ### CHANGE  output link
        Oral_fixed_length          = "Oral_fixed_length",        ### CHANGE  fixed link
        Oral_theta2_rest           = "Oral_theta2_rest",         ### CHANGE
        # Opercular four-bar
        Operculum_length           = "Operculum_length",         ### CHANGE  input link
        IOP_length                 = "IOP_length",               ### CHANGE  coupler part 1
        IOP_ligament_length        = "IOP_ligament_length",      ### CHANGE  coupler part 2
        RA_process_length          = "RA_process_length",        ### CHANGE  output link
        Opercular_fixed_length     = "Opercular_fixed_length",   ### CHANGE  fixed link
        Opercular_theta2_rest      = "Opercular_theta2_rest",    ### CHANGE
        Opercular_theta4_rest_scan = "Opercular_theta4_rest"     ### CHANGE  delete this line if theta4_rest isn't measured in the scans
      )
    
    # Check that every column exists in the scan file before going any further
    missingCols <- setdiff(c(scanFamilyCol, scanAnimalCol, scanColumns), names(scanData))
    if(length(missingCols) > 0) {
      stop("These columns are not in the scan file: ", paste(missingCols, collapse = ", "))
    }
    
    # Make the scan animal names match the format in finalData (A prefix, periods instead of + and -)
    scanData[[scanAnimalCol]] <- as.character(scanData[[scanAnimalCol]])
    fixRow3 <- which(substring(scanData[[scanAnimalCol]],1,1)!='A')
    scanData[[scanAnimalCol]][fixRow3] <- paste0('A',scanData[[scanAnimalCol]][fixRow3])
    scanData[[scanAnimalCol]] <- gsub('\\+','\\.',scanData[[scanAnimalCol]])
    scanData[[scanAnimalCol]] <- gsub('\\-','\\.',scanData[[scanAnimalCol]])
    
    # Make empty columns in finalData to hold the measurements
    for(colNow in names(scanColumns)) {
      finalData[[colNow]] <- NA_real_
    }
    finalData$mandClampFrames <- NA_real_   # frames where the mandible triangle was slightly impossible and had to be clamped
    
    
    for(i in 1:nrow(finalData)){
      
      # Match on family AND animal ID, because the same ID is reused across families
      scanRowNow <- which(scanData[[scanFamilyCol]] == finalData$family[i] & scanData[[scanAnimalCol]] == finalData$Animal.ID[i])
      
      # No scan for this individual: leave the measurements as NA
      if(length(scanRowNow) == 0) {
        next
      }
      
      if(length(scanRowNow) > 1) {
        stop("Scan lookup for ", fileNames2[i], " matched ", length(scanRowNow), " rows in the scan data")
      }
      
      for(colNow in names(scanColumns)) {
        finalData[[colNow]][i] <- as.numeric(scanData[[scanColumns[colNow]]][scanRowNow])
      }
    }
    
    
    # Checks
      # How many videos are missing each measurement?
    print(colSums(is.na(finalData[, names(scanColumns)])))
      # Which individuals have no scan match at all?
    print(unique(finalData$Unique.ID[is.na(finalData$Mandible_length)]))
  
  
  
  
  
    
    
    
###################################
# Four-bar resting angles and KT  #
###################################

        # Are the scan angles in degrees? If so, convert to radians
        scanAnglesInDegrees <- TRUE   ### CHANGE to FALSE if the scan angles are already in radians
        
        if(scanAnglesInDegrees) {
          angleCols <- intersect(c("Oral_theta2_rest", "Opercular_theta2_rest", "Opercular_theta4_rest_scan"), names(finalData))
          for(colNow in angleCols) {
            finalData[[colNow]] <- finalData[[colNow]] * pi / 180
          }
        }
        
        # Direction of opening in four-bar angles
          # Mandible_ang and Neurocranium_rotation are positive for opening/elevation
          # Use +1 if that increases theta2 in your four-bar convention, -1 if it decreases it
        oralOpenSign      <- 1   ### CHECK with one individual
        opercularOpenSign <- 1   ### CHECK with one individual
        
        # Input rotations used for static KT
        oralStaticAngle      <- pi/6    # 30 degrees
        opercularStaticAngle <- pi/18   # 10 degrees
        
        
        # IOP_link_length (opercular coupler link)
        finalData$IOP_link_length <- finalData$IOP_length + finalData$IOP_ligament_length
        
        
        # Link ratios
        finalData$Oral_input_ratio        <- finalData$Mandible_input_length / finalData$Oral_fixed_length
        finalData$Oral_coupler_ratio      <- finalData$Maxilla_length / finalData$Oral_fixed_length
        finalData$Oral_output_ratio       <- finalData$Nasal_length / finalData$Oral_fixed_length
        finalData$Opercular_input_ratio   <- finalData$Operculum_length / finalData$Opercular_fixed_length
        finalData$Opercular_coupler_ratio <- finalData$IOP_link_length / finalData$Opercular_fixed_length
        finalData$Opercular_output_ratio  <- finalData$RA_process_length / finalData$Opercular_fixed_length
        
        
        # Resting angles from the four-bar model
        finalData$Oral_theta3_rest <- fourBarTheta3(finalData$Oral_theta2_rest,
                                                    finalData$Mandible_input_length, finalData$Maxilla_length,
                                                    finalData$Nasal_length, finalData$Oral_fixed_length)
        
        finalData$Opercular_theta4_rest <- fourBarTheta4(finalData$Opercular_theta2_rest,
                                                         finalData$Operculum_length, finalData$IOP_link_length,
                                                         finalData$RA_process_length, finalData$Opercular_fixed_length)
        
        
        # Output angles at the static input rotations
        finalData$Oral_theta3_30 <- fourBarTheta3(finalData$Oral_theta2_rest + oralOpenSign * oralStaticAngle,
                                                  finalData$Mandible_input_length, finalData$Maxilla_length,
                                                  finalData$Nasal_length, finalData$Oral_fixed_length)
        
        finalData$Opercular_theta4_10 <- fourBarTheta4(finalData$Opercular_theta2_rest + opercularOpenSign * opercularStaticAngle,
                                                       finalData$Operculum_length, finalData$IOP_link_length,
                                                       finalData$RA_process_length, finalData$Opercular_fixed_length)
        
        
        # Static KT (average output rotation / input rotation, from rest to the static input angle)
        finalData$Oral_4bar_KT      <- wrapAngle(finalData$Oral_theta3_30 - finalData$Oral_theta3_rest) / (oralOpenSign * oralStaticAngle)
        finalData$Opercular_4bar_KT <- wrapAngle(finalData$Opercular_theta4_10 - finalData$Opercular_theta4_rest) / (opercularOpenSign * opercularStaticAngle)
        
        
        # Checks
          # Individuals whose linkage can't close (NaN). The link lengths or resting angle need checking
        print(unique(finalData$Unique.ID[is.nan(finalData$Oral_theta3_rest) | is.nan(finalData$Oral_theta3_30)]))
        print(unique(finalData$Unique.ID[is.nan(finalData$Opercular_theta4_rest) | is.nan(finalData$Opercular_theta4_10)]))
        
          # Calculated vs. scanned theta4_rest, in degrees. Should be close to 0
        if("Opercular_theta4_rest_scan" %in% names(finalData)) {
          finalData$Opercular_theta4_rest_diff <- wrapAngle(finalData$Opercular_theta4_rest - finalData$Opercular_theta4_rest_scan)
          print(summary(finalData$Opercular_theta4_rest_diff * 180 / pi))
        }
        
        
    
    
    
    
  

  
###########################
########  Analysis ########
###########################


#The plan:
    #I need to correct for two things: 1) jitter in point placement and, 2) misplaced points in individual frames
    
    #For the jitter I will use a kalman filter
      #?kalman_filter
    # For bad landmarks I will use distance between the same point in successive frames
      # Measure distance using a moving measure of median and MAD (mean absolute deviation)
      # This moving measurement is implimented in a Hampel filter


#####################
# Analysis settings #
#####################
  # Created before the loop so they survive the per-video reset and can be used after the loop

  # Smoothing
  lambdaRef        <- 5e-6
  lambdaRefDerived <- lambdaRef/5       # smoothing for derived series (R_t etc.). Should do less smoothing than the original pass
  spanRef          <- 73
  lambdaRefRam     <- lambdaRefDerived  ### CHANGE  smoothing for x_mouth and the nasal point in the Ram section. Acceleration is very sensitive to it
  ######################### RESET THESE BASED ON TEST VIDEO FITS #########################

  # Hyoid
  hyoidNegTolerance <- 0    ### CHANGE  how far below 0 (mm) depression can go before it counts as negative. Raise if tracking noise trips the flag
  betaMinFrames     <- 50   ### CHANGE  minimum number of eye-based strike frames needed to fit the pooled retraction-to-depression ratio

  # Buccal volume
  uFlowAreaFrac <- 0.5      ### CHANGE  big-gape mean of U_flow_ff_predicted only uses frames where A >= this fraction of A_Max

  # QC plots (one PNG per video)
  saveQCplots  <- TRUE                   ### CHANGE  FALSE to skip (they take time and disk space for ~2900 videos)
  qcPlotFolder <- 'PATH/TO/QC_plots'     ### CHANGE

  # finalData columns for each hyoid version. Names on the left match the values returned by hyoidOutputs()
  hyoidColsMain <- c(Time = "Time_hyoid", visible_frac = "Hyoid_visible_frac", vel_mean_visible = "Hyoid_vel_mean_visible",
                     vel_tmax = "Hyoid_vel_tmax", Max = "Hyoid_Max", tmax = "Hyoid_depression_tmax",
                     t_mean = "Hyoid_depression_t_mean", vel_mean = "Hyoid_vel_mean")

  hyoidColsBackup <- c(Time = "Time_hyoid_backup", visible_frac = "Hyoid_visible_frac_backup", vel_mean_visible = "Hyoid_excursion_vel_mean_visible_backup",
                       vel_tmax = "Hyoid_excursion_vel_tmax_backup", Max = "Hyoid_excursion_Max_backup", tmax = "Hyoid_excursion_tmax_backup",
                       t_mean = "Hyoid_excursion_t_mean_backup", vel_mean = "Hyoid_excursion_vel_mean_backup")
  
  
  
  
  
  
  


#Record everything that exists before the analysis loop
  # At the top of each loop iteration, anything NOT in this list is deleted, so no per-video variable can carry over to the next video
  # Anything that needs to survive across videos (finalData, finalTimeSeriesData, skippedVids, functions, settings) must be created ABOVE this line
preLoopVars <- c(ls(), "preLoopVars")

for(i in 1:nVids) {
  
  
  
        ###########################################
        # Clear per-video variables from the last video #
        ###########################################
          # Deletes everything created inside the loop on the previous iteration
          # If a section uses a variable before it's set for this video, R will stop with "object not found" instead of silently using the last video's value

          staleVars <- setdiff(ls(), c(preLoopVars, "i"))
          rm(list = staleVars)
          rm(staleVars)
  
  
        #Get the names of the landmarks from this file and check that they are in the same order as expected
            landmarkNames <- h5read(file=fileNames[i],
                                    name='/node_names')
            expectedNames <- c("UJ", "LJ", "Nasal", "VentMax", "Eye", "Hyoid", "OpercTab")
            
            if(!identical(as.vector(landmarkNames), expectedNames)){
              stop("The landmark names are different from expected in: ", fileNames[i])
            }
            
            #Read in the actual tracking data
            dataNow <- h5read(file=fileNames[i], #Read in the ith data set
                              name='/tracks',
                              index= list(NULL,NULL,NULL,NULL))
                
            #[a,b,x,n]
            #a is the frame number
            #b is the node number (names of the nodes are stored in the "node_names" data set)
            #x holds the x,y coordinates of all the points
            #n is the number of unique individuals that were tracked through the video (always 1 in this data set)
        
        
        
        ##### Hampel filter to find misplaced points and set to NA
          # Run it over x and y for each landmark, replace flagged points with NA
        
        hampFlagsNow <- 0
        
        for(node in 1:length(dataNow[1,,1,1])){
          
          xFlags <- applyHampel(dataNow[,node,1,1])
          yFlags <- applyHampel(dataNow[,node,2,1])
          
          allFlags <- union(xFlags, yFlags)
          
          dataNow[allFlags,node,,1] <- NA
          
          hampFlagsNow <- hampFlagsNow + length(allFlags)
        }
        
        finalData$hampelFlags[i] <- hampFlagsNow
        
        ##### Interpolate missing landmarks using na_interpolation
          # Only fill gaps <= 3 frames so it doesn't try to fill times when the landmark isn't visible
          # Only fill gaps between the first and last real landmarks.
        
        for(node in 1:length(dataNow[1,,1,1])){
          
          xPointsNow <- dataNow[,node,1,1]
          yPointsNow <- dataNow[,node,2,1]
          realFrames <- which(!is.na(xPointsNow) & !is.na(yPointsNow))
          
          if(length(realFrames)>=10) {
            firstReal <- min(realFrames)
            lastReal <- max(realFrames)
            
            interpolatedX <- na_interpolation(xPointsNow[firstReal:lastReal], option="linear", maxgap=3)
            interpolatedY <- na_interpolation(yPointsNow[firstReal:lastReal], option="linear", maxgap=3)
            
            dataNow[firstReal:lastReal,node,1,1] <- interpolatedX
            dataNow[firstReal:lastReal,node,2,1] <- interpolatedY
          }
        
        }
        
        ##### Smooth the data using smooth.spline
          # Helps keep first and second derivatives useful for velocity and acceleration
        
        
        ######################### RESET THESE BASED ON TEST VIDEO FITS #########################
        
        # Make an empty array the same shape as dataNow to hold the smoothed coordinates
          # dataNow stays unsmoothed so raw and smoothed can be compared
        smoothNow <- array(NA_real_, dim = dim(dataNow))
        
        for(node in 1:dim(dataNow)[2]){
          for(xy in 1:2){                          # 1 = x, 2 = y
            
            pointsNow  <- dataNow[,node,xy,1]
            realFrames <- which(!is.na(pointsNow))
            
            # Skip landmarks with too few points to fit
            if(length(realFrames) < 10) {
              next
            }
            
            # Scale lambda by the span of frames being fit, so smoothing is the same in every video
            spanNow   <- max(realFrames) - min(realFrames)
            lambdaNow <- lambdaRef * (spanRef / spanNow)^3
            
            splineFit <- smooth.spline(x = realFrames, y = pointsNow[realFrames], lambda = lambdaNow)
            
            # Only predict at frames that had data, so long gaps and leading/trailing NAs stay NA
            smoothNow[realFrames,node,xy,1] <- predict(splineFit, realFrames)$y
          }
        }
        
        
        
        ### Convert coordinates to mm
        
        smoothNowmm <- smoothNow / finalData$cal_pix.mm[i]
      
        
        ### Make a time variable that is in the units of seconds
          # Time 0 is the first frame
        timeSeriesNow <- 0:(length(smoothNowmm[,1,1,1]) - 1)/finalData$frameRate[i]
        
        
        
        
        
        #landmarkNames
        #LMs in order 1:UJ, 2:LJ, 3:Nasal, 4:Maxilla, 5:Eye, 6:Hyoid, 7:Operculum (these are identities, not necessarily the variable names)
        #[a,b,x,n]
            #a is the frame number
            #b is the node number (names of the nodes are stored in the "node_names" data set)
            #x holds the x,y coordinates of all the points
            #n is the number of unique individuals that were tracked through the video (always 1 in this data set)
        
        xUJ <- smoothNowmm[,1,1,1]
        yUJ <- smoothNowmm[,1,2,1]
        
        xLJ <- smoothNowmm[,2,1,1]
        yLJ <- smoothNowmm[,2,2,1]
        
        xNasal <- smoothNowmm[,3,1,1]
        yNasal <- smoothNowmm[,3,2,1]
        
        xMaxilla <- smoothNowmm[,4,1,1]
        yMaxilla <- smoothNowmm[,4,2,1]
        
        xEye <- smoothNowmm[,5,1,1]
        yEye <- smoothNowmm[,5,2,1]
        
        xHyoid <- smoothNowmm[,6,1,1]
        yHyoid <- smoothNowmm[,6,2,1]
        
        xOperculum <- smoothNowmm[,7,1,1]
        yOperculum <- smoothNowmm[,7,2,1]
        
        
        
        
        
        
        
        ####################
        # Strike direction #
        ####################
            # +1 = fish facing right, -1 = fish facing left
            # Uses the first frame with UJ, LJ and VentMax all present
              # VentMax sits behind the jaw tips, so:
              # Facing right: maxilla is to the LEFT of both UJ and LJ
              # Facing left:  maxilla is to the RIGHT of both UJ and LJ
            # Only x is used, so the y-down image axis doesn't matter here
            # Anything that depends on rotation direction or forward motion (neurocranium, maxilla angle, x_mouth/U_ram) multiplies by this
  
            dfTemp <- data.frame(xUJ=xUJ, xLJ=xLJ, xMaxilla=xMaxilla)
            firstDirFrame <- which(complete.cases(dfTemp))[1]
  
            if(is.na(firstDirFrame)) {
              stop("No frame has UJ, LJ and VentMax all present, so strike direction can't be found in: ", fileNames2[i])
            }
  
            if(xMaxilla[firstDirFrame] < min(xUJ[firstDirFrame], xLJ[firstDirFrame])) {
              strikeDirection <- 1
            } else if(xMaxilla[firstDirFrame] > max(xUJ[firstDirFrame], xLJ[firstDirFrame])) {
              strikeDirection <- -1
            } else {
              stop("VentMax is between UJ and LJ in x at frame ", firstDirFrame, ", so strike direction is ambiguous (check landmark placement) in: ", fileNames2[i])
            }
  
            finalData$strikeDirection[i] <- strikeDirection
          
          
          
        
        
        
        
        ########
        # Gape #
        ########
        
        # R(t)
          R_t <- sqrt((xLJ-xUJ)^2+(yLJ-yUJ)^2)/2     #Radius is half of gape, so it is halved
            #save it to the final array
          finalTimeSeriesData$R_t[1, 1:length(R_t), i] <- R_t
          
          
          
          
          
        ####################
        # Timing Variables #
        ####################
        
        # Tstart and Tend
        
          #Clear residual values from previous video
            peakFrame <- NA
            R_max <- NA
            R_min <- NA
            startThreshold <- NA
            startPos <- NA
            endPos <- NA
            Tstart <- NA
            Tend <- NA
            
        if(any(!is.na(R_t))) {
            
            peakFrame <- which.max(R_t)
            
            R_max <- R_t[peakFrame]
            R_min <- min(R_t[1:peakFrame], na.rm = TRUE)
            
            
            if(any(!is.na(R_t[1:peakFrame]))) {
            
                  #Find the gape threshold that defines the start of the strike
                startThreshold <- R_min + 0.2 * (R_max - R_min)
                
                startPos <- which(R_t>startThreshold)[1]
                Tstart <- timeSeriesNow[startPos]
                
                endPos <- which(R_t[peakFrame:length(R_t)]<startThreshold)[1] + peakFrame - 1 #Have to add the time of peak position because it found which frame after peak was below threshold
                Tend <- timeSeriesNow[endPos]
                
                if (!is.na(Tstart) && is.na(Tend)) {
                  stop("No strike end found (gape never drops below threshold after peak) in: ", fileNames2[i])
                }
            
            }
        }
            
            finalData$Tstart[i] <- Tstart
            finalData$Tend[i] <- Tend
        
        
        # Ttpg
          finalData$Ttpg[i] <- timeSeriesNow[peakFrame] - Tstart
          
        # Ttotal
          finalData$Ttotal[i] <- Tend - Tstart
          
        # R_t_mean
          R_t_mean <- NA
          if(!any(is.na(c(startPos,endPos)))) {
          R_t_mean <- mean(R_t[startPos:endPos], na.rm=TRUE)
          }
          finalData$R_t_mean[i] <- R_t_mean
      
            if(!is.na(finalData$tmax[i])) {
      
              goodFrames_gape <- which(!is.na(R_t))
                
                
                #Check if there are enough frames for fitting, and that tmax is within the frames that have gape data
                if(length(goodFrames_gape) >= 10 && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_gape)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_gape)]) {
                  
                    # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                    spanNow   <- max(goodFrames_gape) - min(goodFrames_gape)
                    lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3
                    
                    R_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_gape], y = R_t[goodFrames_gape], lambda = lambdaNow)
                    
                    # R_tmax (mm)
                    finalData$R_tmax[i] <- predict(R_t_smooth, x = finalData$tmax[i], deriv = 0)$y
                    
                    # R_vel_tmax (mm/s)
                    finalData$R_vel_tmax[i] <- predict(R_t_smooth, x = finalData$tmax[i], deriv = 1)$y
                }
          }
          
        
        # Relative_ttpg
        
        finalData$Relative_ttpg[i] <- finalData$Ttpg[i] / finalData$Ttotal[i]
        
        
        #Relative_tmax
        
        finalData$Relative_tmax[i] <- (finalData$tmax[i] - finalData$Tstart[i]) / finalData$Ttotal[i]
        
        
        
        ##### Record skipped analyses: Gape and timing
          # Strikes with no Tstart are only logged here, not again in every later section

          if(is.na(finalData$Tstart[i])) {
            if(all(is.na(R_t))) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Gape", "UJ and LJ never both present")
            } else {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Gape", "No strike start found from gape threshold")
            }
          } else if(!is.na(finalData$tmax[i]) && is.na(finalData$R_tmax[i])) {
            skippedVids <- addSkip(skippedVids, fileNames2[i], "Gape", "Values at tmax not calculated: tmax outside tracked frames or < 10 frames")
          }
        
        
        
        
        
        ##############
        # Gape area  #
        ##############

        # A(t) = (pi * W_jaw / 2) * R(t): ellipse with half-width W_jaw/2 (tooth row, from scans) and half-height R(t)
          # W_jaw is fixed for each fish, so every gape area value is R times a constant
          # Values at tmax, the mean and the max come straight from the R values, so the SEM equations are exact identities
          # Assumes the jaws don't widen during the strike
        # Variables
          # A_t:        gape area over time (mm^2)
          # A_tmax:     (pi * W_jaw / 2) * R_tmax (mm^2)
          # A_vel_tmax: (pi * W_jaw / 2) * R_vel_tmax (mm^2/s)
          # A_mean:     (pi * W_jaw / 2) * R_t_mean, mean over Tstart-Tend (mm^2)
          # A_Max:      (pi * W_jaw / 2) * largest R during the strike (mm^2)
        # NA if W_jaw is missing. Not logged, since missing scans aren't a landmarking problem

          gapeAreaConst <- pi * finalData$W_jaw[i] / 2

          if(!is.na(gapeAreaConst)) {

            # A_t (mm^2)
            A_t <- gapeAreaConst * R_t
              #save it to the final array
            finalTimeSeriesData$A_t[1, 1:length(A_t), i] <- A_t

            # A_tmax (mm^2) and A_vel_tmax (mm^2/s). NA if R wasn't calculated at tmax
            finalData$A_tmax[i]     <- gapeAreaConst * finalData$R_tmax[i]
            finalData$A_vel_tmax[i] <- gapeAreaConst * finalData$R_vel_tmax[i]

            # A_mean (mm^2). NA if there's no Tstart
            finalData$A_mean[i] <- gapeAreaConst * finalData$R_t_mean[i]

            # A_Max (mm^2)
            if(!is.na(finalData$Tstart[i])) {
              finalData$A_Max[i] <- gapeAreaConst * max(R_t[startPos:endPos], na.rm=TRUE)
            }
          }
        
        
        
        
        
        
        ##############  
        # Protrusion #
        ##############
        
        # Find protrusion in first frame with both nasal and UJ points. Set this as zero point for protrusion
        
          dfTemp <- data.frame(xUJ=xUJ, yUJ=yUJ, xNasal=xNasal, yNasal=yNasal)
          firstProtFrame <- which(complete.cases(dfTemp))[1]
          
          Protrusion_t_mean <- NA
          
          
          if(!is.na(firstProtFrame) && !is.na(finalData$Tstart[i]) && (timeSeriesNow[firstProtFrame])<=finalData$Tstart[i]) {
          
            # Protrusion_t
              
              Protrusion_firstProtFrame <- sqrt((xUJ[firstProtFrame]-xNasal[firstProtFrame])^2+(yUJ[firstProtFrame]-yNasal[firstProtFrame])^2)
              Protrusion_t <- sqrt((xUJ-xNasal)^2+(yUJ-yNasal)^2)-Protrusion_firstProtFrame
                #save it to the final array
              finalTimeSeriesData$Protrusion_t[1, 1:length(Protrusion_t), i] <- Protrusion_t
            
            # Protrusion_t_mean
              
              if(!any(is.na(c(startPos,endPos)))) {
              Protrusion_t_mean <- mean(Protrusion_t[startPos:endPos], na.rm=TRUE)
              }
              
              finalData$Protrusion_t_mean[i] <- Protrusion_t_mean
              
              
            # Protrusion_Max
              if(any(!is.na(Protrusion_t[startPos:endPos]))) {
                  Protrusion_Max <- max(Protrusion_t[startPos:endPos], na.rm=TRUE)
                  finalData$Protrusion_Max[i] <- Protrusion_Max
              }
            
                    if(!is.na(finalData$tmax[i])) {
              
                        goodFrames_prot <- which(!is.na(Protrusion_t))
                        
                        
                        #Check if there are enough frames for fitting, and that tmax is within the frames that have protrusion data
                        if(length(goodFrames_prot) >= 10 && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_prot)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_prot)]) {
                      
                            # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                            spanNow   <- max(goodFrames_prot) - min(goodFrames_prot)
                            lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3
                            
                            Protrusion_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_prot], y = Protrusion_t[goodFrames_prot], lambda = lambdaNow)
                            
                            # Protrusion_tmax (mm)
                            finalData$Protrusion_tmax[i] <- predict(Protrusion_t_smooth, x = finalData$tmax[i], deriv = 0)$y
                            
                            # Protrusion_vel_tmax (mm/s)
                            finalData$Protrusion_vel_tmax[i] <- predict(Protrusion_t_smooth, x = finalData$tmax[i], deriv = 1)$y
                        }
                    }
          }
        
        
                ##### Record skipped analyses: Protrusion

          if(!is.na(finalData$Tstart[i])) {
            if(is.na(firstProtFrame)) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Protrusion", "UJ and Nasal never both present")
            } else if(timeSeriesNow[firstProtFrame] > finalData$Tstart[i]) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Protrusion", "No rest frame: first frame with UJ and Nasal is after Tstart")
            } else if(!is.na(finalData$tmax[i]) && is.na(finalData$Protrusion_tmax[i])) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Protrusion", "Values at tmax not calculated: tmax outside tracked frames or < 10 frames")
            }
          }
          

        ############
        # Mandible #
        ############

        # Mandible rotation from the triangle Nasal - jaw joint - LJ
          # D = Nasal tip to jaw joint distance (Nasal_joint_length, from scans)
          # r = jaw joint to LJ tip (Mandible_length, from scans)
          # d = Nasal to LJ distance in each frame (from the video)
          # Law of cosines gives the angle at the jaw joint: phi = acos((D^2 + r^2 - d^2) / (2*D*r))
        # Variables
          # Mandible_tip_disp_t:   change in Nasal to LJ distance from rest (mm). Measured directly, no scan values needed
          # Mandible_ang_t:        change in phi from rest (radians, opening is positive)
          # Mandible_depression_t: Mandible_length * sin(Mandible_ang_t) (mm). Identity with the SEM equation

          dfTemp <- data.frame(xLJ=xLJ, yLJ=yLJ, xNasal=xNasal, yNasal=yNasal)
          firstMandFrame <- which(complete.cases(dfTemp))[1]


          if(!is.na(firstMandFrame) && !is.na(finalData$Tstart[i]) && (timeSeriesNow[firstMandFrame])<=finalData$Tstart[i]) {

            ##### Mandible_tip_disp_t

              nasalLJ_t <- sqrt((xLJ-xNasal)^2+(yLJ-yNasal)^2)
              Mandible_tip_disp_t <- nasalLJ_t - nasalLJ_t[firstMandFrame]
                #save it to the final array
              finalTimeSeriesData$Mandible_tip_disp_t[1, 1:length(Mandible_tip_disp_t), i] <- Mandible_tip_disp_t
              finalData$Mandible_tip_dist_rest[i] <- nasalLJ_t[firstMandFrame]
              
              
              
            ##### Mandible_tip_disp_Max (largest change in Nasal-LJ distance from rest during the strike, mm)
                # No scan values needed, so this exists even for fish without scans
  
                if(any(!is.na(Mandible_tip_disp_t[startPos:endPos]))) {
                  finalData$Mandible_tip_disp_Max[i] <- max(Mandible_tip_disp_t[startPos:endPos], na.rm=TRUE)
                }

            

                    if(!is.na(finalData$tmax[i])) {

                        goodFrames_mandTip <- which(!is.na(Mandible_tip_disp_t))

                        #Check if there are enough frames for fitting, and that tmax is within the frames that have mandible data
                        if(length(goodFrames_mandTip) >= 10 && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_mandTip)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_mandTip)]) {

                            # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                            spanNow   <- max(goodFrames_mandTip) - min(goodFrames_mandTip)
                            lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

                            Mandible_tip_disp_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_mandTip], y = Mandible_tip_disp_t[goodFrames_mandTip], lambda = lambdaNow)

                            # Mandible_tip_disp_tmax (mm)
                            finalData$Mandible_tip_disp_tmax[i] <- predict(Mandible_tip_disp_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                            # Mandible_tip_disp_vel_tmax (mm/s)
                            finalData$Mandible_tip_disp_vel_tmax[i] <- predict(Mandible_tip_disp_t_smooth, x = finalData$tmax[i], deriv = 1)$y
                        }
                    }


            ##### Mandible_ang_t and Mandible_depression_t (need scan values)

              jawD <- finalData$Nasal_joint_length[i]
              jawR <- finalData$Mandible_length[i]

              if(!is.na(jawD) && !is.na(jawR)) {

                cosPhi <- (jawD^2 + jawR^2 - nasalLJ_t^2) / (2*jawD*jawR)

                # The rest frame must be a possible triangle. If not, the scan values or the calibration are wrong
                if(abs(cosPhi[firstMandFrame]) > 1) {
                  stop("Nasal-LJ distance at rest (", round(nasalLJ_t[firstMandFrame],2), " mm) is impossible for the scanned jaw length and nasal-joint distance in: ", fileNames2[i])
                }

                # Tracking noise can push cosPhi slightly past +/-1 during the strike. Count those frames, then clamp them
                finalData$mandClampFrames[i] <- sum(abs(cosPhi) > 1, na.rm=TRUE)
                cosPhi <- pmin(pmax(cosPhi, -1), 1)

                phi_t <- acos(cosPhi)

                # Mandible_ang_t (radians)
                Mandible_ang_t <- phi_t - phi_t[firstMandFrame]
                  #save it to the final array
                finalTimeSeriesData$Mandible_ang_t[1, 1:length(Mandible_ang_t), i] <- Mandible_ang_t

                # Mandible_depression_t (mm)
                Mandible_depression_t <- jawR * sin(Mandible_ang_t)
                  #save it to the final array
                finalTimeSeriesData$Mandible_depression_t[1, 1:length(Mandible_depression_t), i] <- Mandible_depression_t


                if(any(!is.na(Mandible_ang_t[startPos:endPos]))) {

                  # Mandible_depression_t_mean
                    finalData$Mandible_depression_t_mean[i] <- mean(Mandible_depression_t[startPos:endPos], na.rm=TRUE)

                  # Mandible_Max (largest jaw rotation during the strike, radians)
                    finalData$Mandible_Max[i] <- max(Mandible_ang_t[startPos:endPos], na.rm=TRUE)
                    
                    
                  # Mandible_depression_Max (mm)
                    # Depression = Mandible_length * sin(angle) keeps increasing with angle up to 90 degrees, so its max is at Mandible_Max
                    # Calculated from Mandible_Max so it's an exact identity, not a separate max search
                    # A jaw opening past 90 degrees is impossible, so it means bad tracking, scan values or calibration
                    if(finalData$Mandible_Max[i] > pi/2) {
                      stop("Mandible_Max is over 90 degrees (", round(finalData$Mandible_Max[i]*180/pi, 1), " deg) in: ", fileNames2[i])
                    }
                    finalData$Mandible_depression_Max[i] <- jawR * sin(finalData$Mandible_Max[i])

                  
                }


                    if(!is.na(finalData$tmax[i])) {

                        goodFrames_mandAng <- which(!is.na(Mandible_ang_t))

                        #Check if there are enough frames for fitting, and that tmax is within the frames that have mandible data
                        if(length(goodFrames_mandAng) >= 10 && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_mandAng)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_mandAng)]) {

                            # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                            spanNow   <- max(goodFrames_mandAng) - min(goodFrames_mandAng)
                            lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

                            Mandible_ang_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_mandAng], y = Mandible_ang_t[goodFrames_mandAng], lambda = lambdaNow)

                            # Mandible_ang_tmax (radians)
                            finalData$Mandible_ang_tmax[i] <- predict(Mandible_ang_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                            # Mandible_angular_vel_tmax (radians/s)
                            finalData$Mandible_angular_vel_tmax[i] <- predict(Mandible_ang_t_smooth, x = finalData$tmax[i], deriv = 1)$y

                            # Depression and its velocity come from the angle, so the SEM equations are exact identities
                            # Mandible_depression_tmax (mm)
                            finalData$Mandible_depression_tmax[i] <- jawR * sin(finalData$Mandible_ang_tmax[i])

                            # Mandible_dep_vel_tmax (mm/s)
                            finalData$Mandible_dep_vel_tmax[i] <- jawR * finalData$Mandible_angular_vel_tmax[i] * cos(finalData$Mandible_ang_tmax[i])
                        }
                    }
              }
          }
          
          
          
         ##### Record skipped analyses: Mandible

          if(!is.na(finalData$Tstart[i])) {
            if(is.na(firstMandFrame)) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Mandible", "LJ and Nasal never both present")
            } else if(timeSeriesNow[firstMandFrame] > finalData$Tstart[i]) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Mandible", "No rest frame: first frame with LJ and Nasal is after Tstart")
            } else {

              # Tip displacement (no scan values needed)
              if(!is.na(finalData$tmax[i]) && is.na(finalData$Mandible_tip_disp_tmax[i])) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Mandible", "Tip displacement at tmax not calculated: tmax outside tracked frames or < 10 frames")
              }

              # Angle (only checked when scan values exist, since missing scans aren't a landmarking problem)
              if(!is.na(finalData$Nasal_joint_length[i]) && !is.na(finalData$Mandible_length[i])) {
                if(!is.na(finalData$tmax[i]) && is.na(finalData$Mandible_ang_tmax[i])) {
                  skippedVids <- addSkip(skippedVids, fileNames2[i], "Mandible", "Angle at tmax not calculated: tmax outside tracked frames or < 10 frames")
                }

                # Flag (not a skip): frames where the triangle was impossible and had to be clamped
                if(!is.na(finalData$mandClampFrames[i]) && finalData$mandClampFrames[i] > 0) {
                  skippedVids <- addSkip(skippedVids, fileNames2[i], "Mandible", paste0("FLAG: ", finalData$mandClampFrames[i], " frames clamped in the mandible triangle"))
                }
              }
            }
          }
          
          
          
          
        ######################
        # Oral 4-bar at tmax #
        ######################
          # Input angle at tmax = resting input angle + measured mandible rotation
        
          if(!is.na(finalData$Mandible_ang_tmax[i]) && !is.na(finalData$Oral_theta2_rest[i])) {
            
            oralIn  <- finalData$Mandible_input_length[i]
            oralCp  <- finalData$Maxilla_length[i]
            oralOut <- finalData$Nasal_length[i]
            oralFx  <- finalData$Oral_fixed_length[i]
            
            finalData$Oral_theta2_tmax[i] <- finalData$Oral_theta2_rest[i] + oralOpenSign * finalData$Mandible_ang_tmax[i]
            finalData$Oral_theta3_tmax[i] <- fourBarTheta3(finalData$Oral_theta2_tmax[i], oralIn, oralCp, oralOut, oralFx)
            finalData$Oral_theta4_tmax[i] <- fourBarTheta4(finalData$Oral_theta2_tmax[i], oralIn, oralCp, oralOut, oralFx)
            
            # Oral_KT_tmax (instantaneous d theta3 / d theta2)
            finalData$Oral_KT_tmax[i] <- fourBarKTinst(finalData$Oral_theta2_tmax[i], oralIn, oralCp, oralOut, oralFx, output = "theta3")
            
            # Maxilla_ang_predicted_tmax (radians, change in theta3 from rest)
            finalData$Maxilla_ang_predicted_tmax[i] <- wrapAngle(finalData$Oral_theta3_tmax[i] - finalData$Oral_theta3_rest[i])
          }

          
          
          
          
          
          
          
        ################
        # Neurocranium #
        ################

        # Neurocranium rotation from the angle of the Eye -> Nasal vector
          # On a rigid skull, the line between any two points on it rotates by exactly the skull's rotation,
            # no matter where the pivot (craniovertebral joint) is. So the eye works as a base point without knowing the pivot
          # Eye sliding ALONG the Eye-Nasal line doesn't change the angle
          # Eye sliding PERPENDICULAR to it does (error = slide / Eye-Nasal distance)
          # Whole-body pitch is included in this angle. Known source of error, assumed small because strikes are fast
        # Variables
          # Neurocranium_rotation_t:      change in Eye->Nasal angle from rest (radians, elevation is positive)
          # Neurocranium_tip_disp_t:      Neurocranium_length * Neurocranium_rotation_t (mm). Arc length, so rotation = tip_disp / Neurocranium_length exactly
          # Neurocranium_tip_disp_eye_t:  nasal tip displacement relative to the eye, perpendicular to the resting Eye-Nasal line (mm). No scan values needed
          # Neurocranium_angular_vel_t:   d(rotation)/dt from a spline (radians/s)
          # Neurocranium_linear_vel_t:    Neurocranium_length * Neurocranium_angular_vel_t (mm/s). Identical to the spline derivative of tip_disp
          # Neurocranium_eyeNasal_stretch_max: QC. Largest proportional change in Eye-Nasal distance during the strike
            # Should be 0 for a rigid skull. Only detects eye sliding along the line, so it's a proxy for sliding in general

          eyeStretchCutoff <- 0.20   ### CHANGE  lenient QC cutoff for flagging eye sliding (proportion of resting Eye-Nasal distance)

          dfTemp <- data.frame(xNasal=xNasal, yNasal=yNasal, xEye=xEye, yEye=yEye)
          firstNeuroFrame <- which(complete.cases(dfTemp))[1]


          if(!is.na(firstNeuroFrame) && !is.na(finalData$Tstart[i]) && (timeSeriesNow[firstNeuroFrame])<=finalData$Tstart[i]) {

            ##### Eye -> Nasal vector in each frame
              # y is flipped so up is positive (image y points down)

              vxNeuro <- xNasal - xEye
              vyNeuro <- -(yNasal - yEye)

              vxRest <- vxNeuro[firstNeuroFrame]
              vyRest <- vyNeuro[firstNeuroFrame]
              eyeNasalRest <- sqrt(vxRest^2 + vyRest^2)

              # The nasal must be in front of the eye in the direction the fish is facing. If not, landmarks are swapped or the direction is wrong
              if(sign(vxRest) != finalData$strikeDirection[i]) {
                stop("Nasal is not in front of the eye at rest (frame ", firstNeuroFrame, ") for strikeDirection = ", finalData$strikeDirection[i], ". Check Nasal/Eye/VentMax landmarks in: ", fileNames2[i])
              }


            ##### Neurocranium_rotation_t (radians)
              # Signed angle from the resting vector to the current vector, atan2(cross, dot). Positive = counterclockwise (with y up)
              # Elevation is counterclockwise for a fish facing right and clockwise for a fish facing left, so multiply by strikeDirection

              Neurocranium_rotation_t <- finalData$strikeDirection[i] * atan2(vxRest*vyNeuro - vyRest*vxNeuro,
                                                                              vxRest*vxNeuro + vyRest*vyNeuro)
                #save it to the final array
              finalTimeSeriesData$Neurocranium_rotation_t[1, 1:length(Neurocranium_rotation_t), i] <- Neurocranium_rotation_t


            ##### Neurocranium_tip_disp_eye_t (mm)
              # Unit vector perpendicular to the resting Eye-Nasal line, pointing dorsally for either facing direction

              xPerp <- -finalData$strikeDirection[i] * vyRest / eyeNasalRest
              yPerp <-  finalData$strikeDirection[i] * vxRest / eyeNasalRest

              Neurocranium_tip_disp_eye_t <- (vxNeuro - vxRest)*xPerp + (vyNeuro - vyRest)*yPerp
                #save it to the final array
              finalTimeSeriesData$Neurocranium_tip_disp_eye_t[1, 1:length(Neurocranium_tip_disp_eye_t), i] <- Neurocranium_tip_disp_eye_t


            ##### QC: Eye-Nasal stretch during the strike

              eyeNasalStretch_t <- sqrt(vxNeuro^2 + vyNeuro^2) / eyeNasalRest - 1

              if(any(!is.na(eyeNasalStretch_t[startPos:endPos]))) {
                finalData$Neurocranium_eyeNasal_stretch_max[i] <- max(abs(eyeNasalStretch_t[startPos:endPos]), na.rm=TRUE)
                finalData$Neurocranium_eyeQC_flag[i] <- finalData$Neurocranium_eyeNasal_stretch_max[i] > eyeStretchCutoff
              }


            ##### Neurocranium_Max and Time_cranial
              # Time_cranial is measured from Tstart, same as Ttpg and Time_hyoid

              if(any(!is.na(Neurocranium_rotation_t[startPos:endPos]))) {

                cranialPeakFrame <- startPos - 1 + which.max(Neurocranium_rotation_t[startPos:endPos])

                # Neurocranium_Max (largest cranial elevation during the strike, radians)
                finalData$Neurocranium_Max[i] <- Neurocranium_rotation_t[cranialPeakFrame]

                # Time_cranial (s)
                finalData$Time_cranial[i] <- timeSeriesNow[cranialPeakFrame] - finalData$Tstart[i]
              }


            ##### Neurocranium_tip_disp_t (mm, needs scan value)

              neuroL <- finalData$Neurocranium_length[i]

              if(!is.na(neuroL)) {
                Neurocranium_tip_disp_t <- neuroL * Neurocranium_rotation_t
                  #save it to the final array
                finalTimeSeriesData$Neurocranium_tip_disp_t[1, 1:length(Neurocranium_tip_disp_t), i] <- Neurocranium_tip_disp_t
              }


            ##### Spline on rotation: velocity time series, mean velocities, and values at tmax

              goodFrames_neuro <- which(!is.na(Neurocranium_rotation_t))

              if(length(goodFrames_neuro) >= 10) {

                  # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                  spanNow   <- max(goodFrames_neuro) - min(goodFrames_neuro)
                  lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

                  Neurocranium_rotation_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_neuro], y = Neurocranium_rotation_t[goodFrames_neuro], lambda = lambdaNow)

                  # Neurocranium_angular_vel_t (radians/s). Only at frames that had data
                  Neurocranium_angular_vel_t <- rep(NA_real_, length(Neurocranium_rotation_t))
                  Neurocranium_angular_vel_t[goodFrames_neuro] <- predict(Neurocranium_rotation_t_smooth, x = timeSeriesNow[goodFrames_neuro], deriv = 1)$y
                    #save it to the final array
                  finalTimeSeriesData$Neurocranium_angular_vel_t[1, 1:length(Neurocranium_angular_vel_t), i] <- Neurocranium_angular_vel_t

                  # Neurocranium_linear_vel_t (mm/s, needs scan value)
                  if(!is.na(neuroL)) {
                    Neurocranium_linear_vel_t <- neuroL * Neurocranium_angular_vel_t
                      #save it to the final array
                    finalTimeSeriesData$Neurocranium_linear_vel_t[1, 1:length(Neurocranium_linear_vel_t), i] <- Neurocranium_linear_vel_t
                  }


                  # Mean velocities from Tstart to peak cranial elevation
                    # Over the whole strike these would be ~0, because the cranium elevates and then returns
                    # Plain means. Muscle/MA scaling for <Neurocranium_Adj_vel> is done outside this script
                  if(!is.na(finalData$Time_cranial[i]) && any(!is.na(Neurocranium_angular_vel_t[startPos:cranialPeakFrame]))) {

                    # Neurocranium_angular_vel_mean (radians/s)
                    finalData$Neurocranium_angular_vel_mean[i] <- mean(Neurocranium_angular_vel_t[startPos:cranialPeakFrame], na.rm=TRUE)

                    # Neurocranium_linear_vel_mean (mm/s). NA if no scan value
                    finalData$Neurocranium_linear_vel_mean[i] <- neuroL * finalData$Neurocranium_angular_vel_mean[i]
                  }


                  if(!is.na(finalData$tmax[i])) {

                      #Check that tmax is within the frames that have neurocranium data
                      if(finalData$tmax[i] >= timeSeriesNow[min(goodFrames_neuro)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_neuro)]) {

                          # Neurocranium_rotation_tmax (radians)
                          finalData$Neurocranium_rotation_tmax[i] <- predict(Neurocranium_rotation_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                          # Neurocranium_angular_vel_tmax (radians/s)
                          finalData$Neurocranium_angular_vel_tmax[i] <- predict(Neurocranium_rotation_t_smooth, x = finalData$tmax[i], deriv = 1)$y

                          # Tip displacement and linear velocity come from the rotation spline, so the SEM equations are exact identities
                            # NA if no scan value
                          # Neurocranium_tip_disp_tmax (mm)
                          finalData$Neurocranium_tip_disp_tmax[i] <- neuroL * finalData$Neurocranium_rotation_tmax[i]

                          # Neurocranium_linear_vel_tmax (mm/s)
                          finalData$Neurocranium_linear_vel_tmax[i] <- neuroL * finalData$Neurocranium_angular_vel_tmax[i]


                          # Neurocranium_tip_disp_eye_tmax (mm). Separate spline, same frames and smoothing
                          Neurocranium_tip_disp_eye_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_neuro], y = Neurocranium_tip_disp_eye_t[goodFrames_neuro], lambda = lambdaNow)
                          finalData$Neurocranium_tip_disp_eye_tmax[i] <- predict(Neurocranium_tip_disp_eye_t_smooth, x = finalData$tmax[i], deriv = 0)$y
                      }
                  }
              }
          }
          
          
          
          
        ##### Record skipped analyses: Neurocranium

          if(!is.na(finalData$Tstart[i])) {
            if(is.na(firstNeuroFrame)) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Neurocranium", "Eye and Nasal never both present")
            } else if(timeSeriesNow[firstNeuroFrame] > finalData$Tstart[i]) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Neurocranium", "No rest frame: first frame with Eye and Nasal is after Tstart")
            } else {

              if(sum(!is.na(Neurocranium_rotation_t)) < 10) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Neurocranium", "Fewer than 10 frames with Eye and Nasal: no velocities or values at tmax")
              } else if(!is.na(finalData$tmax[i]) && is.na(finalData$Neurocranium_rotation_tmax[i])) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Neurocranium", "Values at tmax not calculated: tmax outside tracked frames")
              }

              # Flag (not a skip): eye slid relative to the nasal more than the QC cutoff
              if(isTRUE(finalData$Neurocranium_eyeQC_flag[i])) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Neurocranium", paste0("FLAG: Eye-Nasal stretch ", round(finalData$Neurocranium_eyeNasal_stretch_max[i], 3), " is above the QC cutoff"))
              }
            }
          }
          
          
          
          
        ###########
        # Maxilla #
        ###########

        # Maxilla rotation from the angle of the Nasal -> VentMax vector
          # The nasal tip is about where the dorsal end of the maxilla meets the nasal, so this vector approximates the maxilla (oral 4-bar coupler)
          # Signed angle from rest using atan2(cross, dot), multiplied by strikeDirection
            # Positive = ventral tip swinging anteriorly (counterclockwise for a fish facing right, clockwise for facing left)
            # This is the same rotation direction as cranial elevation, so Neurocranium_rotation_t can be subtracted directly
        # Two versions of each variable
          # Skull frame (Maxilla_ang_t, Maxilla_tip_disp_t): neurocranium rotation removed. Comparable to Maxilla_ang_predicted_tmax. Needs the eye
          # Camera frame (Maxilla_ang_cam_t, Maxilla_tip_disp_cam_t): includes neurocranium rotation. No eye needed, so available in more videos
        # Variables
          # Maxilla_ang_cam_t:       change in Nasal->VentMax angle from rest, camera frame (radians)
          # Maxilla_ang_t:           Maxilla_ang_cam_t - Neurocranium_rotation_t (radians)
          # Maxilla_tip_disp_cam_t:  distance VentMax has moved from its resting position relative to the nasal tip, camera frame (mm)
          # Maxilla_tip_disp_t:      same, after rotating out the neurocranium rotation (mm)
            # Displacements are distances, so they are always >= 0, and tracking noise biases them slightly positive near rest

          dfTemp <- data.frame(xNasal=xNasal, yNasal=yNasal, xMaxilla=xMaxilla, yMaxilla=yMaxilla)
          firstMaxFrame <- which(complete.cases(dfTemp))[1]


          if(!is.na(firstMaxFrame) && !is.na(finalData$Tstart[i]) && (timeSeriesNow[firstMaxFrame])<=finalData$Tstart[i]) {

            ##### Nasal -> VentMax vector in each frame
              # y is flipped so up is positive (image y points down)

              vxMax <- xMaxilla - xNasal
              vyMax <- -(yMaxilla - yNasal)

              vxMaxRest <- vxMax[firstMaxFrame]
              vyMaxRest <- vyMax[firstMaxFrame]

              # VentMax must be below the nasal tip at rest. If not, the landmarks are swapped or misplaced
              if(vyMaxRest >= 0) {
                stop("VentMax is not ventral to the nasal tip at rest (frame ", firstMaxFrame, "). Check Nasal/VentMax landmarks in: ", fileNames2[i])
              }


            ##### Maxilla_ang_cam_t (radians)

              Maxilla_ang_cam_t <- finalData$strikeDirection[i] * atan2(vxMaxRest*vyMax - vyMaxRest*vxMax,
                                                                        vxMaxRest*vxMax + vyMaxRest*vyMax)
                #save it to the final array
              finalTimeSeriesData$Maxilla_ang_cam_t[1, 1:length(Maxilla_ang_cam_t), i] <- Maxilla_ang_cam_t


            ##### Maxilla_tip_disp_cam_t (mm)

              Maxilla_tip_disp_cam_t <- sqrt((vxMax - vxMaxRest)^2 + (vyMax - vyMaxRest)^2)
                #save it to the final array
              finalTimeSeriesData$Maxilla_tip_disp_cam_t[1, 1:length(Maxilla_tip_disp_cam_t), i] <- Maxilla_tip_disp_cam_t


            ##### Maxilla_cam_Max (largest camera-frame maxilla rotation during the strike, radians)

              if(any(!is.na(Maxilla_ang_cam_t[startPos:endPos]))) {
                finalData$Maxilla_cam_Max[i] <- max(Maxilla_ang_cam_t[startPos:endPos], na.rm=TRUE)
              }


            ##### Camera-frame values at tmax

              if(!is.na(finalData$tmax[i])) {

                  goodFrames_maxCam <- which(!is.na(Maxilla_ang_cam_t))

                  #Check if there are enough frames for fitting, and that tmax is within the frames that have maxilla data
                  if(length(goodFrames_maxCam) >= 10 && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_maxCam)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_maxCam)]) {

                      # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                      spanNow   <- max(goodFrames_maxCam) - min(goodFrames_maxCam)
                      lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

                      Maxilla_ang_cam_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_maxCam], y = Maxilla_ang_cam_t[goodFrames_maxCam], lambda = lambdaNow)

                      # Maxilla_ang_cam_tmax (radians)
                      finalData$Maxilla_ang_cam_tmax[i] <- predict(Maxilla_ang_cam_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                      # Maxilla_angular_vel_cam_tmax (radians/s)
                      finalData$Maxilla_angular_vel_cam_tmax[i] <- predict(Maxilla_ang_cam_t_smooth, x = finalData$tmax[i], deriv = 1)$y


                      Maxilla_tip_disp_cam_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_maxCam], y = Maxilla_tip_disp_cam_t[goodFrames_maxCam], lambda = lambdaNow)

                      # Maxilla_tip_disp_cam_tmax (mm)
                      finalData$Maxilla_tip_disp_cam_tmax[i] <- predict(Maxilla_tip_disp_cam_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                      # Maxilla_linear_vel_cam_tmax (mm/s)
                      finalData$Maxilla_linear_vel_cam_tmax[i] <- predict(Maxilla_tip_disp_cam_t_smooth, x = finalData$tmax[i], deriv = 1)$y
                  }
              }


            ##### Skull-frame versions (need Neurocranium_rotation_t, which needs the eye)
              # Neurocranium_rotation_t only exists if the neurocranium section ran for THIS video (the loop reset deletes it otherwise)

              if(exists("Neurocranium_rotation_t", inherits = FALSE)) {

                ##### Maxilla_ang_t (radians)

                  Maxilla_ang_t <- Maxilla_ang_cam_t - Neurocranium_rotation_t
                    #save it to the final array
                  finalTimeSeriesData$Maxilla_ang_t[1, 1:length(Maxilla_ang_t), i] <- Maxilla_ang_t


                ##### Maxilla_tip_disp_t (mm)
                  # Rotate the Nasal -> VentMax vector back by the skull's rotation, so cranial elevation doesn't count as maxilla motion
                  # Skull rotation in camera coordinates (counterclockwise positive) is strikeDirection * Neurocranium_rotation_t
                  # The rest vector needs no rotation, because the skull is at rest in the rest frame

                  skullRot <- finalData$strikeDirection[i] * Neurocranium_rotation_t

                  vxMaxSkull <-  vxMax*cos(skullRot) + vyMax*sin(skullRot)
                  vyMaxSkull <- -vxMax*sin(skullRot) + vyMax*cos(skullRot)

                  Maxilla_tip_disp_t <- sqrt((vxMaxSkull - vxMaxRest)^2 + (vyMaxSkull - vyMaxRest)^2)
                    #save it to the final array
                  finalTimeSeriesData$Maxilla_tip_disp_t[1, 1:length(Maxilla_tip_disp_t), i] <- Maxilla_tip_disp_t


                ##### Maxilla_Max (largest skull-frame maxilla rotation during the strike, radians)

                  if(any(!is.na(Maxilla_ang_t[startPos:endPos]))) {
                    finalData$Maxilla_Max[i] <- max(Maxilla_ang_t[startPos:endPos], na.rm=TRUE)
                  }


                ##### Skull-frame values at tmax

                  if(!is.na(finalData$tmax[i])) {

                      goodFrames_max <- which(!is.na(Maxilla_ang_t))

                      #Check if there are enough frames for fitting, and that tmax is within the frames that have maxilla and eye data
                      if(length(goodFrames_max) >= 10 && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_max)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_max)]) {

                          # Same lambda scaling as the coordinate smoothing, based on the frame span being fit
                          spanNow   <- max(goodFrames_max) - min(goodFrames_max)
                          lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

                          Maxilla_ang_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_max], y = Maxilla_ang_t[goodFrames_max], lambda = lambdaNow)

                          # Maxilla_ang_tmax (radians)
                          finalData$Maxilla_ang_tmax[i] <- predict(Maxilla_ang_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                          # Maxilla_angular_vel_tmax (radians/s)
                          finalData$Maxilla_angular_vel_tmax[i] <- predict(Maxilla_ang_t_smooth, x = finalData$tmax[i], deriv = 1)$y


                          Maxilla_tip_disp_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_max], y = Maxilla_tip_disp_t[goodFrames_max], lambda = lambdaNow)

                          # Maxilla_tip_disp_tmax (mm)
                          finalData$Maxilla_tip_disp_tmax[i] <- predict(Maxilla_tip_disp_t_smooth, x = finalData$tmax[i], deriv = 0)$y

                          # Maxilla_linear_vel_tmax (mm/s)
                          finalData$Maxilla_linear_vel_tmax[i] <- predict(Maxilla_tip_disp_t_smooth, x = finalData$tmax[i], deriv = 1)$y
                      }
                  }
              }
          } 
          
          
          
        ##### Record skipped analyses: Maxilla

          if(!is.na(finalData$Tstart[i])) {
            if(is.na(firstMaxFrame)) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Maxilla", "VentMax and Nasal never both present")
            } else if(timeSeriesNow[firstMaxFrame] > finalData$Tstart[i]) {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Maxilla", "No rest frame: first frame with VentMax and Nasal is after Tstart")
            } else {

              # Camera frame
              if(!is.na(finalData$tmax[i]) && is.na(finalData$Maxilla_ang_cam_tmax[i])) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Maxilla", "Camera-frame values at tmax not calculated: tmax outside tracked frames or < 10 frames")
              }

              # Skull frame (needs the neurocranium rotation, so needs the eye)
              if(!exists("Maxilla_ang_t", inherits = FALSE)) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Maxilla", "Skull-frame maxilla skipped: no neurocranium rotation (eye missing)")
              } else if(!is.na(finalData$tmax[i]) && is.na(finalData$Maxilla_ang_tmax[i])) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Maxilla", "Skull-frame values at tmax not calculated: tmax outside frames with eye and maxilla, or < 10 frames")
              }
            }
          }
          
          
          
          
          
          
          
          
          
          
        #########
        # Hyoid #
        #########

        # Hyoid landmark = anterior tip of the urohyal, seen through the skin
          # Only visible mid-strike, so resting positions come from the scans
        # Main version (Hyoid_*): dorsoventral depth below the Eye-Nasal line, minus D_head
          # Videos with the eye: measured here, in the skull frame (line taken from the same frame, so neurocranium rotation is already removed)
          # Videos without the eye: calculated after the loop from the Nasal-Hyoid distance with the triangle and the pooled retraction ratio
          # Hyoid_depth_source records which ("eye" or "triangle")
        # Hyoid_retraction_t (eye videos): posterior movement along the Eye-Nasal line, minus Hyoid_rest_AP
        # Backup (Hyoid_excursion_*_backup): change in Nasal-Hyoid distance from D_head_backup
          # Not a depression: mixes retraction with dorsoventral motion. Kept as a separate phenotype
        # All per-version calculations are in hyoidOutputs() (functions section)

          # Saved for the after-loop pass (triangle depth and buccal volume), which works from the saved time series
          finalData$nFrames[i]    <- length(timeSeriesNow)
          finalData$startFrame[i] <- startPos
          finalData$endFrame[i]   <- endPos


          ##### Skull frame: depth below, and position along, the Eye-Nasal line (mm)
            # y is flipped so up is positive (image y points down)

            vxEN <- xNasal - xEye
            vyEN <- -(yNasal - yEye)
            lengthEN <- sqrt(vxEN^2 + vyEN^2)

            # Unit vector perpendicular to the Eye->Nasal line, pointing ventrally for either facing direction
            xVentral <-  finalData$strikeDirection[i] * vyEN / lengthEN
            yVentral <- -finalData$strikeDirection[i] * vxEN / lengthEN

            # Unit vector along the line, pointing posteriorly (from the nasal toward the eye)
            xPosterior <- -vxEN / lengthEN
            yPosterior <- -vyEN / lengthEN

            hyoidDepthSkull_t <- (xHyoid - xEye)*xVentral   - (yHyoid - yEye)*yVentral     # depth below the line
            hyoidAlongSkull_t <- (xHyoid - xEye)*xPosterior - (yHyoid - yEye)*yPosterior   # distance behind the eye, along the line


          ##### Hyoid_retraction_t (mm, eye videos only)
            # Positive = moved posteriorly from rest

            if(!is.na(finalData$Hyoid_rest_AP[i])) {
              Hyoid_retraction_t <- hyoidAlongSkull_t - finalData$Hyoid_rest_AP[i]
                #save it to the final array
              finalTimeSeriesData$Hyoid_retraction_t[1, 1:length(Hyoid_retraction_t), i] <- Hyoid_retraction_t
            }


          ##### Nasal-Hyoid distance (mm). Used for the excursion backup here, and for the triangle after the loop

            hyoidDistNasal_t <- sqrt((xHyoid - xNasal)^2 + (yHyoid - yNasal)^2)
              #save it to the final array
            finalTimeSeriesData$Hyoid_nasal_dist_t[1, 1:length(hyoidDistNasal_t), i] <- hyoidDistNasal_t


          # Settings for each version
          hyoidVersions <- list(
            list(section = "Hyoid",        missing = "Hyoid, Eye and Nasal never all present", raw_t = hyoidDepthSkull_t, rest = finalData$D_head[i],
                 cols = hyoidColsMain,   depressionSeries = "Hyoid_depression_t",        velSeries = "Hyoid_vel_t"),
            list(section = "Hyoid backup", missing = "Hyoid and Nasal never both present",   raw_t = hyoidDistNasal_t,  rest = finalData$D_head_backup[i],
                 cols = hyoidColsBackup, depressionSeries = "Hyoid_excursion_t_backup",  velSeries = "Hyoid_excursion_vel_t_backup")
          )


          # Strikes with no Tstart are already logged in the Gape section
          if(!is.na(finalData$Tstart[i])) {

            for(hv in hyoidVersions) {

              hyOut <- hyoidOutputs(hv$raw_t, hv$rest, timeSeriesNow, startPos, endPos, finalData$tmax[i], finalData$Tstart[i],
                                    lambdaRefDerived, spanRef, hyoidNegTolerance, hv$missing)

              # Save values and time series
              for(nm in names(hv$cols)) {
                finalData[[hv$cols[[nm]]]][i] <- hyOut$values[[nm]]
              }
              if(!is.null(hyOut$vel_t)) {
                finalTimeSeriesData[[hv$velSeries]][1, 1:length(hyOut$vel_t), i] <- hyOut$vel_t
              }
              if(!is.null(hyOut$depression_t)) {
                finalTimeSeriesData[[hv$depressionSeries]][1, 1:length(hyOut$depression_t), i] <- hyOut$depression_t
              }

              # Log
              for(msgNow in hyOut$log) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], hv$section, msgNow)
              }
            }

            # Eye-based if the skull-frame version had any strike frames. Otherwise the triangle is tried after the loop
            if(!is.na(finalData$Time_hyoid[i])) {
              finalData$Hyoid_depth_source[i] <- "eye"
            }
          }
          
        
          
          
          
          
        #######
        # Ram #
        #######

        # Strike axis = direction of the body's net movement over the strike, from the Nasal position at Tstart to Tend
          # Mouth velocity is projected onto this axis, so:
            # Up/down motion of the gape center from jaw opening (perpendicular to the strike path) is not counted as ram
            # A tilted strike path is measured at full speed
          # Falls back to the camera x-axis if the Nasal isn't tracked at Tstart and Tend, or the body barely moves
          # Signed so that moving toward the prey is positive
          # PIV velocity is a speed magnitude, assumed collinear with ram. The kinematic camera is mirrored relative to the PIV camera,
            # but that doesn't matter here: ram is signed by strikeDirection and the PIV value has no sign
          # Assumes the strike path is roughly straight over Tstart-Tend
        # x_mouth(t) = midpoint of UJ and LJ (Muller et al. gape center)
          # Includes jaw protrusion (along the strike path) as well as body motion, which is correct for the pressure equation
        # Variables
          # x_mouth_t:          position of the gape center along the strike axis (mm, increasing = moving toward the prey)
          # U_ram_t:            d(x_mouth)/dt from a spline (mm/s)
          # U_ram_mean:         (x_mouth(Tend) - x_mouth(Tstart)) / Ttotal. Exact identity for (1/T) * integral of U_ram
          # U_ram_mean_frames:  plain mean of U_ram_t over Tstart-Tend frames. Check, should be nearly equal to U_ram_mean
          # U_ram_body_*:       same, using the Nasal point projected onto the same axis. Separates whole-body motion from jaw-driven motion of the mouth
            # Nasal moves slightly with cranial elevation, so this still includes a small part of head motion
          # Cov_R_Uram:         (1/T) * integral of (R - <R>)(U_ram - <U_ram>) over Tstart-Tend (mm^2/s)
          # strikeAxisTilt:     angle of the strike axis above horizontal (radians, positive = strike path tilted upward)
          # strikeAxisFromNasal: TRUE if the axis came from the Nasal, FALSE if it fell back to the camera x-axis

          strikeAxisMinDisp <- 0.5   ### CHANGE  minimum net Nasal displacement over the strike (mm) for the strike axis to be meaningful


          # Strikes with no Tstart are already logged in the Gape section
          if(!is.na(finalData$Tstart[i])) {

            ##### Strike axis
              # y is flipped so up is positive (image y points down)

              # Default: camera x-axis, pointing toward the prey
              xAxis <- finalData$strikeDirection[i]
              yAxis <- 0
              finalData$strikeAxisFromNasal[i] <- FALSE

              nasalDispX <- xNasal[endPos] - xNasal[startPos]
              nasalDispY <- -(yNasal[endPos] - yNasal[startPos])

              if(is.na(nasalDispX) || is.na(nasalDispY)) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram", "FLAG: Nasal not tracked at Tstart or Tend, strike axis set to camera x-axis")

              } else if(sqrt(nasalDispX^2 + nasalDispY^2) < strikeAxisMinDisp) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram", paste0("FLAG: body moved less than ", strikeAxisMinDisp, " mm, strike axis set to camera x-axis"))

              } else if(nasalDispX * finalData$strikeDirection[i] == 0) {
                skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram", "FLAG: body moved straight up or down, strike axis set to camera x-axis")

              } else {
                nasalDispLength <- sqrt(nasalDispX^2 + nasalDispY^2)

                # Point the axis toward the prey. If the fish backed up over the strike, the axis still points forward and ram comes out negative
                axisSign <- sign(nasalDispX * finalData$strikeDirection[i])
                xAxis <- axisSign * nasalDispX / nasalDispLength
                yAxis <- axisSign * nasalDispY / nasalDispLength
                finalData$strikeAxisFromNasal[i] <- TRUE
              }

              # strikeAxisTilt (radians above horizontal, in the forward direction)
              finalData$strikeAxisTilt[i] <- atan2(yAxis, finalData$strikeDirection[i] * xAxis)


            ##### x_mouth_t (mm along the strike axis)

              xMouthAxis_t <- ((xUJ + xLJ) / 2) * xAxis + (-(yUJ + yLJ) / 2) * yAxis
                #save it to the final array
              finalTimeSeriesData$x_mouth_t[1, 1:length(xMouthAxis_t), i] <- xMouthAxis_t

              goodFrames_mouth <- which(!is.na(xMouthAxis_t))


            ##### Mouth ram

              if(length(goodFrames_mouth) >= 10) {

                  # Same lambda scaling as the coordinate smoothing, based on the frame span being fit. Uses the ram-specific lambda
                  spanNow   <- max(goodFrames_mouth) - min(goodFrames_mouth)
                  lambdaNow <- lambdaRefRam * (spanRef / spanNow)^3

                  x_mouth_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_mouth], y = xMouthAxis_t[goodFrames_mouth], lambda = lambdaNow)

                  # U_ram_t (mm/s). Only at frames that had data
                  U_ram_t <- rep(NA_real_, length(xMouthAxis_t))
                  U_ram_t[goodFrames_mouth] <- predict(x_mouth_t_smooth, x = timeSeriesNow[goodFrames_mouth], deriv = 1)$y
                    #save it to the final array
                  finalTimeSeriesData$U_ram_t[1, 1:length(U_ram_t), i] <- U_ram_t


                  # U_ram_mean (mm/s)
                    # Uses spline positions at Tstart and Tend, so it exactly equals the time-average of the spline's velocity
                    # UJ and LJ always exist at Tstart and Tend, because those frames are defined by gape
                  finalData$U_ram_mean[i] <- (predict(x_mouth_t_smooth, x = timeSeriesNow[endPos], deriv = 0)$y -
                                              predict(x_mouth_t_smooth, x = timeSeriesNow[startPos], deriv = 0)$y) / finalData$Ttotal[i]

                  # U_ram_mean_frames (mm/s)
                  finalData$U_ram_mean_frames[i] <- mean(U_ram_t[startPos:endPos], na.rm=TRUE)


                  # Cov_R_Uram (mm^2/s)
                    # Population covariance (divide by n) to match (1/T) * integral. Means are taken over the same frames as the covariance
                  covFrames <- intersect(startPos:endPos, which(!is.na(R_t) & !is.na(U_ram_t)))

                  if(length(covFrames) >= 2) {
                    finalData$Cov_R_Uram[i] <- mean((R_t[covFrames] - mean(R_t[covFrames])) * (U_ram_t[covFrames] - mean(U_ram_t[covFrames])))
                  }


                  # Cov_R_Uflowefmeas and Cov_R_Uflowff: add once the PIV flow time series are imported   ### CHANGE
                    # Use the same frames as Cov_R_Uram, so Cov_R_Uflowff = Cov_R_Uflowefmeas - Cov_R_Uram is an exact identity
                    # PIV is in m/s, so convert to mm/s (* 1000) first
                  # covFramesPIV <- intersect(covFrames, which(!is.na(U_flow_ef_meas_t)))
                  # finalData$Cov_R_Uflowefmeas[i] <- mean((R_t[covFramesPIV] - mean(R_t[covFramesPIV])) * (U_flow_ef_meas_t[covFramesPIV] - mean(U_flow_ef_meas_t[covFramesPIV])))
                  # finalData$Cov_R_Uflowff[i]     <- finalData$Cov_R_Uflowefmeas[i] - finalData$Cov_R_Uram[i]


                  # Values at tmax
                    # x_mouth has the same frames as gape, so a tmax outside the frames is already logged in the Gape section
                  if(!is.na(finalData$tmax[i]) && finalData$tmax[i] >= timeSeriesNow[min(goodFrames_mouth)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_mouth)]) {

                      # U_ram_tmax (mm/s)
                      finalData$U_ram_tmax[i] <- predict(x_mouth_t_smooth, x = finalData$tmax[i], deriv = 1)$y

                      # U_ram_acc_tmax (mm/s^2)
                      finalData$U_ram_acc_tmax[i] <- predict(x_mouth_t_smooth, x = finalData$tmax[i], deriv = 2)$y
                  }

              } else {
                  skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram", "Fewer than 10 frames with UJ and LJ: no ram values")
              }


            ##### Body ram (Nasal point, projected onto the same strike axis)

              xNasalAxis_t <- xNasal * xAxis + (-yNasal) * yAxis
              goodFrames_nasalRam <- which(!is.na(xNasalAxis_t))

              if(length(goodFrames_nasalRam) >= 10) {

                  # Same lambda scaling, using the ram-specific lambda
                  spanNow   <- max(goodFrames_nasalRam) - min(goodFrames_nasalRam)
                  lambdaNow <- lambdaRefRam * (spanRef / spanNow)^3

                  x_nasal_t_smooth <- smooth.spline(x = timeSeriesNow[goodFrames_nasalRam], y = xNasalAxis_t[goodFrames_nasalRam], lambda = lambdaNow)

                  # U_ram_body_t (mm/s). Only at frames that had data
                  U_ram_body_t <- rep(NA_real_, length(xNasalAxis_t))
                  U_ram_body_t[goodFrames_nasalRam] <- predict(x_nasal_t_smooth, x = timeSeriesNow[goodFrames_nasalRam], deriv = 1)$y
                    #save it to the final array
                  finalTimeSeriesData$U_ram_body_t[1, 1:length(U_ram_body_t), i] <- U_ram_body_t

                  # U_ram_body_mean (mm/s)
                    # Only if the nasal is tracked at both Tstart and Tend. Predicting the spline outside its frames would extrapolate
                  if(timeSeriesNow[startPos] >= timeSeriesNow[min(goodFrames_nasalRam)] && timeSeriesNow[endPos] <= timeSeriesNow[max(goodFrames_nasalRam)]) {
                    finalData$U_ram_body_mean[i] <- (predict(x_nasal_t_smooth, x = timeSeriesNow[endPos], deriv = 0)$y -
                                                     predict(x_nasal_t_smooth, x = timeSeriesNow[startPos], deriv = 0)$y) / finalData$Ttotal[i]
                  } else {
                    skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram body", "Nasal not tracked at Tstart or Tend: U_ram_body_mean not calculated")
                  }

                  # U_ram_body_mean_frames (mm/s)
                  if(any(!is.na(U_ram_body_t[startPos:endPos]))) {
                    finalData$U_ram_body_mean_frames[i] <- mean(U_ram_body_t[startPos:endPos], na.rm=TRUE)
                  }

                  # Values at tmax
                  if(!is.na(finalData$tmax[i])) {
                    if(finalData$tmax[i] >= timeSeriesNow[min(goodFrames_nasalRam)] && finalData$tmax[i] <= timeSeriesNow[max(goodFrames_nasalRam)]) {

                        # U_ram_body_tmax (mm/s)
                        finalData$U_ram_body_tmax[i] <- predict(x_nasal_t_smooth, x = finalData$tmax[i], deriv = 1)$y

                        # U_ram_body_acc_tmax (mm/s^2)
                        finalData$U_ram_body_acc_tmax[i] <- predict(x_nasal_t_smooth, x = finalData$tmax[i], deriv = 2)$y

                    } else {
                        skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram body", "Values at tmax not calculated: tmax outside tracked Nasal frames")
                    }
                  }

              } else {
                  skippedVids <- addSkip(skippedVids, fileNames2[i], "Ram body", "Fewer than 10 frames with Nasal: no body ram values")
              }
          }
          
          
          
          
          
          
            
        ###########################
        # Opercular 4-bar at tmax #
        ###########################
          # Input angle at tmax = resting input angle + measured neurocranium rotation (proxy for opercular rotation)
        
          if(!is.na(finalData$Neurocranium_rotation_tmax[i]) && !is.na(finalData$Opercular_theta2_rest[i])) {
            
            operIn  <- finalData$Operculum_length[i]
            operCp  <- finalData$IOP_link_length[i]
            operOut <- finalData$RA_process_length[i]
            operFx  <- finalData$Opercular_fixed_length[i]
            
            finalData$Opercular_theta2_tmax[i] <- finalData$Opercular_theta2_rest[i] + opercularOpenSign * finalData$Neurocranium_rotation_tmax[i]
            finalData$Opercular_theta3_tmax[i] <- fourBarTheta3(finalData$Opercular_theta2_tmax[i], operIn, operCp, operOut, operFx)
            finalData$Opercular_theta4_tmax[i] <- fourBarTheta4(finalData$Opercular_theta2_tmax[i], operIn, operCp, operOut, operFx)
            
            # Opercular_KT_tmax (instantaneous d theta4 / d theta2)
            finalData$Opercular_KT_tmax[i] <- fourBarKTinst(finalData$Opercular_theta2_tmax[i], operIn, operCp, operOut, operFx, output = "theta4")
            
            # Mandible_ang_predicted_tmax (radians, change in theta4 from rest)
            finalData$Mandible_ang_predicted_tmax[i] <- wrapAngle(finalData$Opercular_theta4_tmax[i] - finalData$Opercular_theta4_rest[i])
          }
          





        
        
}
  
  









        
        
###################################################
# Hyoid without the eye (triangle) and buccal volume #
###################################################
  # Runs after the loop because:
    # The triangle needs the pooled retraction ratio (beta), which needs every eye-based video first
    # Buccal volume needs hyoid depression and retraction from both eye-based and triangle videos
  # Works from the saved time series. Frame numbers come from nFrames, startFrame and endFrame saved in the loop


  ##### Pooled retraction-to-depression ratio (beta)
    # Straight line through zero, Hyoid_retraction = beta * Hyoid_depression, fit on strike frames of eye-based videos
    # Pooled on purpose, so fish with extreme retraction are underestimated rather than exaggerated
    # Only videos whose depression passed the negative check (Hyoid_Max not NA)

    betaH <- c()
    betaR <- c()
    betaPhase <- c()   # opening (before peak depression) or closing (after)

    eyeRows <- which(finalData$Hyoid_depth_source == "eye" & !is.na(finalData$Hyoid_Max))

    for(i in eyeRows) {
      framesNow <- finalData$startFrame[i]:finalData$endFrame[i]
      hNow <- finalTimeSeriesData$Hyoid_depression_t[1, framesNow, i]
      rNow <- finalTimeSeriesData$Hyoid_retraction_t[1, framesNow, i]
      keepNow <- !is.na(hNow) & !is.na(rNow)

      peakFrameNow <- finalData$startFrame[i] + round(finalData$Time_hyoid[i] * finalData$frameRate[i])

      betaH     <- c(betaH, hNow[keepNow])
      betaR     <- c(betaR, rNow[keepNow])
      betaPhase <- c(betaPhase, ifelse(framesNow[keepNow] <= peakFrameNow, "opening", "closing"))
    }

    if(length(betaH) >= betaMinFrames) {

      hyoidBeta <- sum(betaR * betaH) / sum(betaH^2)
      betaResid <- betaR - hyoidBeta * betaH
      betaR2    <- 1 - sum(betaResid^2) / sum((betaR - mean(betaR))^2)   # centered R2, a stricter test than the uncentered one for a line through zero

    } else if(!is.null(testVids)) {

      # Test mode: too few eye-based frames in the test videos, so use the fixed test value
      hyoidBeta <- testBeta
      betaResid <- numeric(0)
      betaPhase <- character(0)
      betaR2    <- NA_real_
      print(paste0("TEST MODE: only ", length(betaH), " eye-based frames, so using testBeta = ", testBeta))

    } else {
      stop("Only ", length(betaH), " eye-based strike frames with both hyoid depression and retraction. Need at least betaMinFrames = ", betaMinFrames, " to fit the retraction ratio")
    }


  ##### Triangle depth for videos without the eye
    # Tracked Nasal-Hyoid distance d satisfies d^2 = (a + beta*h)^2 + (D_head + h)^2
      # h = depression, a = Nasal_eye_AP_length + Hyoid_rest_AP (resting AP distance from the nasal tip to the hyoid)
    # Quadratic in h: (1 + beta^2) h^2 + 2(a*beta + D_head) h + (a^2 + D_head^2 - d^2) = 0
      # The + root is the physical one (it's 0 at rest)
    # Retraction = beta * h

    quadA <- 1 + hyoidBeta^2

    triRows <- which(is.na(finalData$Hyoid_depth_source) & !is.na(finalData$Tstart))

    for(i in triRows) {

      aNow     <- finalData$Nasal_eye_AP_length[i] + finalData$Hyoid_rest_AP[i]
      dHeadNow <- finalData$D_head[i]

      # Missing scan values aren't logged, since they aren't a landmarking problem
      if(is.na(aNow) || is.na(dHeadNow)) {
        next
      }

      nNow    <- finalData$nFrames[i]
      timeNow <- (0:(nNow - 1)) / finalData$frameRate[i]
      dNow    <- finalTimeSeriesData$Hyoid_nasal_dist_t[1, 1:nNow, i]

      quadB <- 2 * (aNow * hyoidBeta + dHeadNow)
      quadC <- aNow^2 + dHeadNow^2 - dNow^2
      discNow <- quadB^2 - 4 * quadA * quadC

      # No real solution means the tracked distance is shorter than the closest the hyoid can get to the nasal tip in this model
      if(any(discNow < 0, na.rm = TRUE)) {
        stop("Triangle has no solution (Nasal-Hyoid distance too short for the scan geometry) in: ", fileNames2[i])
      }

      hTri     <- (-quadB + sqrt(discNow)) / (2 * quadA)
      depthTri <- dHeadNow + hTri

      hyOut <- hyoidOutputs(depthTri, dHeadNow, timeNow, finalData$startFrame[i], finalData$endFrame[i], finalData$tmax[i], finalData$Tstart[i],
                            lambdaRefDerived, spanRef, hyoidNegTolerance, "Hyoid and Nasal never both present")

      # Save into the main hyoid columns
      for(nm in names(hyoidColsMain)) {
        finalData[[hyoidColsMain[[nm]]]][i] <- hyOut$values[[nm]]
      }
      if(!is.null(hyOut$vel_t)) {
        finalTimeSeriesData$Hyoid_vel_t[1, 1:nNow, i] <- hyOut$vel_t
      }
      if(!is.null(hyOut$depression_t)) {
        finalTimeSeriesData$Hyoid_depression_t[1, 1:nNow, i] <- hyOut$depression_t
        finalTimeSeriesData$Hyoid_retraction_t[1, 1:nNow, i] <- hyoidBeta * hyOut$depression_t
      }

      for(msgNow in hyOut$log) {
        skippedVids <- addSkip(skippedVids, fileNames2[i], "Hyoid triangle", msgNow)
      }

      finalData$Hyoid_depth_source[i] <- "triangle"
    }


  ##### Buccal volume
    # Ceratohyal bar of fixed length: anterior end moves with the urohyal (retraction r, depression h), posterior end only moves laterally
      # Ceratohyal_half_width(t) = sqrt(L^2 - (Ceratohyal_AP_rest - r)^2 - (Ceratohyal_DV_rest + h)^2)
      # Ceratohyal_lateral_angle_hypaxial(t) = max(0, asin(half_width / L) - asin((W_head/2) / L))
      # Ceratohyal_lateral_angle(t) = Neurocranium_rotation(t) + hypaxial angle
      # Buccal_width_expansion(t) = max(0, 2L[sin(rest angle + lateral angle) - sin(rest angle)])
    # Width can't drop below rest, so both angles are floored at 0. Frames where each floor applies are counted
    # V_buccal(t) = frustum: front = A(t), back = pi(W_head + width expansion)/4 * (D_head + h), length = Neurocranium_length + Protrusion(t)
    # Over Tstart-Tend only. Hyoid depression, retraction and neurocranium rotation are 0 at Tstart/Tend where missing, with gaps filled linearly
      # Protrusion and gape area gaps are filled linearly, and missing ends take the nearest value
    # No eye (triangle videos): neurocranium rotation is set to 0, recorded in V_buccal_neuroIncluded
    # Only videos whose hyoid depression passed the negative check (Hyoid_Max not NA)

    bvRows <- which(!is.na(finalData$Hyoid_Max) & !is.na(finalData$Hyoid_depth_source))

    for(i in bvRows) {

      cerL   <- finalData$Ceratohyal_length[i]
      cerAP0 <- finalData$Ceratohyal_AP_rest[i]
      cerDV0 <- finalData$Ceratohyal_DV_rest[i]
      wHead  <- finalData$W_head[i]
      dHead  <- finalData$D_head[i]
      ncL    <- finalData$Neurocranium_length[i]

      # Missing scan values aren't logged, since they aren't a landmarking problem
      if(any(is.na(c(cerL, cerAP0, cerDV0, wHead, dHead, ncL)))) {
        next
      }

      # Resting geometry must be possible. If not, a scan measurement is wrong for this fish
      if(cerL^2 < cerAP0^2 + cerDV0^2) {
        stop("Ceratohyal_length is shorter than its resting AP and DV span combined for ", finalData$Unique.ID[i], ". Check the ceratohyal scan measurements")
      }
      if(wHead/2 > cerL) {
        stop("W_head/2 is longer than Ceratohyal_length for ", finalData$Unique.ID[i], ". Check the scan measurements")
      }

      cerPsi0 <- asin((wHead/2) / cerL)   # resting lateral angle (radians)

      framesNow <- finalData$startFrame[i]:finalData$endFrame[i]
      nStrike   <- length(framesNow)
      timeNow   <- (framesNow - 1) / finalData$frameRate[i]


      ##### Inputs over Tstart-Tend

        # Hyoid depression and retraction: 0 at Tstart/Tend where missing, gaps filled linearly
        hNow <- finalTimeSeriesData$Hyoid_depression_t[1, framesNow, i]
        rNow <- finalTimeSeriesData$Hyoid_retraction_t[1, framesNow, i]

        if(all(is.na(rNow))) {
          next   # no Hyoid_rest_AP scan value
        }

        if(is.na(hNow[1]))       { hNow[1] <- 0 }
        if(is.na(hNow[nStrike])) { hNow[nStrike] <- 0 }
        if(is.na(rNow[1]))       { rNow[1] <- 0 }
        if(is.na(rNow[nStrike])) { rNow[nStrike] <- 0 }
        hNow <- na_interpolation(hNow, option = "linear")
        rNow <- na_interpolation(rNow, option = "linear")

        # Neurocranium rotation: same filling. Set to 0 if there's none (no eye)
        neuroNow <- finalTimeSeriesData$Neurocranium_rotation_t[1, framesNow, i]

        if(all(is.na(neuroNow))) {
          neuroNow <- rep(0, nStrike)
          finalData$V_buccal_neuroIncluded[i] <- FALSE
        } else {
          if(is.na(neuroNow[1]))       { neuroNow[1] <- 0 }
          if(is.na(neuroNow[nStrike])) { neuroNow[nStrike] <- 0 }
          neuroNow <- na_interpolation(neuroNow, option = "linear")
          finalData$V_buccal_neuroIncluded[i] <- TRUE
        }

        # Protrusion and gape area: gaps filled linearly, missing ends take the nearest value
        protNow <- finalTimeSeriesData$Protrusion_t[1, framesNow, i]
        areaNow <- finalTimeSeriesData$A_t[1, framesNow, i]

        if(sum(!is.na(protNow)) < 2) {
          skippedVids <- addSkip(skippedVids, fileNames2[i], "Buccal volume", "Fewer than 2 strike frames with protrusion: no buccal volume")
          next
        }
        if(sum(!is.na(areaNow)) < 2) {
          next   # no W_jaw scan value
        }
        protNow <- na_interpolation(protNow, option = "linear")
        areaNow <- na_interpolation(areaNow, option = "linear")


      ##### Ceratohyal angles and buccal width (radians, mm)

        if(any(hNow / cerL > 1)) {
          stop("Hyoid depression is longer than Ceratohyal_length (max ", round(max(hNow), 2), " mm) in: ", fileNames2[i], ". Check the D_head zero or Ceratohyal_length")
        }

        # Ceratohyal_rotation_angle_t (phenotype only, doesn't feed the width)
        Ceratohyal_rotation_angle_t <- asin(hNow / cerL)

        # Half-width. A negative square means the bar can't reach, the extreme case of narrowing, so it's floored below
        halfWidthSq <- cerL^2 - (cerAP0 - rNow)^2 - (cerDV0 + hNow)^2
        halfWidth   <- sqrt(pmax(halfWidthSq, 0))

        # Hypaxial lateral angle, floored at 0 (posterior ends can't move medially past rest)
        latHypNow <- asin(halfWidth / cerL) - cerPsi0
        hypFloor  <- latHypNow < 0
        Ceratohyal_lateral_angle_hypaxial_t <- pmax(latHypNow, 0)

        # Effective half-width after the floor
        Ceratohyal_half_width_t <- pmax(halfWidth, wHead/2)

        # Total lateral angle
        Ceratohyal_lateral_angle_t <- neuroNow + Ceratohyal_lateral_angle_hypaxial_t

        # Buccal width expansion, floored at 0 (catches small negative neurocranium rotation)
        widthExpNow <- 2 * cerL * (sin(cerPsi0 + Ceratohyal_lateral_angle_t) - sin(cerPsi0))
        totalFloor  <- widthExpNow < 0
        Buccal_width_expansion_t <- pmax(widthExpNow, 0)

        finalData$Width_hypaxial_floor_frames[i] <- sum(hypFloor)
        finalData$Width_total_floor_frames[i]    <- sum(totalFloor)


      ##### V_buccal_t (mm^3)

        areaBackNow <- pi * (wHead + Buccal_width_expansion_t) / 4 * (dHead + hNow)
        V_buccal_t  <- ((ncL + protNow) / 3) * (areaNow + areaBackNow + sqrt(areaNow * areaBackNow))

        finalData$V_buccal_Max[i]               <- max(V_buccal_t)
        finalData$Buccal_width_expansion_Max[i] <- max(Buccal_width_expansion_t)


      ##### Save time series (strike frames only)

        finalTimeSeriesData$Ceratohyal_rotation_angle_t[1, framesNow, i]         <- Ceratohyal_rotation_angle_t
        finalTimeSeriesData$Ceratohyal_lateral_angle_hypaxial_t[1, framesNow, i] <- Ceratohyal_lateral_angle_hypaxial_t
        finalTimeSeriesData$Ceratohyal_lateral_angle_t[1, framesNow, i]          <- Ceratohyal_lateral_angle_t
        finalTimeSeriesData$Ceratohyal_half_width_t[1, framesNow, i]             <- Ceratohyal_half_width_t
        finalTimeSeriesData$Buccal_width_expansion_t[1, framesNow, i]            <- Buccal_width_expansion_t
        finalTimeSeriesData$V_buccal_t[1, framesNow, i]                          <- V_buccal_t


      ##### Spline on V_buccal: rate of change, values at tmax, predicted flow

        if(nStrike >= 10) {

          # Same lambda scaling as the other derived series
          spanNow   <- nStrike - 1
          lambdaNow <- lambdaRefDerived * (spanRef / spanNow)^3

          V_buccal_t_smooth <- smooth.spline(x = timeNow, y = V_buccal_t, lambda = lambdaNow)

          # V_buccal_vel_t (mm^3/s)
          V_buccal_vel_t <- predict(V_buccal_t_smooth, x = timeNow, deriv = 1)$y
          finalTimeSeriesData$V_buccal_vel_t[1, framesNow, i] <- V_buccal_vel_t

          # U_flow_ff_predicted_t (mm/s) = V_buccal_vel / A
          U_flow_ff_predicted_t <- V_buccal_vel_t / areaNow
          finalTimeSeriesData$U_flow_ff_predicted_t[1, framesNow, i] <- U_flow_ff_predicted_t

          # U_flow_ff_predicted_mean (mm/s): all strike frames
          finalData$U_flow_ff_predicted_mean[i] <- mean(U_flow_ff_predicted_t)

          # U_flow_ff_predicted_mean_bigGape (mm/s): only frames with A >= uFlowAreaFrac * A_Max, so small gapes can't dominate
          bigGapeFrames <- which(areaNow >= uFlowAreaFrac * finalData$A_Max[i])
          if(length(bigGapeFrames) > 0) {
            finalData$U_flow_ff_predicted_mean_bigGape[i] <- mean(U_flow_ff_predicted_t[bigGapeFrames])
          }

          # Values at tmax
          if(!is.na(finalData$tmax[i])) {
            if(finalData$tmax[i] >= timeNow[1] && finalData$tmax[i] <= timeNow[nStrike]) {

              # V_buccal_vel_tmax (mm^3/s)
              finalData$V_buccal_vel_tmax[i] <- predict(V_buccal_t_smooth, x = finalData$tmax[i], deriv = 1)$y

              # V_buccal_acc_tmax (mm^3/s^2)
              finalData$V_buccal_acc_tmax[i] <- predict(V_buccal_t_smooth, x = finalData$tmax[i], deriv = 2)$y

            } else {
              skippedVids <- addSkip(skippedVids, fileNames2[i], "Buccal volume", "Values at tmax not calculated: tmax outside Tstart-Tend")
            }
          }

        } else {
          skippedVids <- addSkip(skippedVids, fileNames2[i], "Buccal volume", "Fewer than 10 strike frames: no buccal volume rates")
        }
    }


  ##### Checks

    # Retraction ratio fit
  print(hyoidBeta)
  print(length(unique(eyeRows)))          # eye-based videos used
  print(length(betaH))                     # frames used
  print(betaR2)                            # how well a line through zero fits. Low = retraction isn't proportional to depression
  print(tapply(betaResid, betaPhase, mean))   # mean residual on opening vs closing. A clear difference = different paths, so one beta is a compromise

    # How many videos used each depth source
  print(table(finalData$Hyoid_depth_source, useNA = "ifany"))

    # Scan consistency, one row per fish
  scanRowsNow <- !duplicated(finalData$Unique.ID)
      # Resting nasal-to-hyoid distance from the triangle inputs vs. D_head_backup. Should be close (they differ only by the ceratohyal stand-in for depth)
  print(summary(sqrt((finalData$Nasal_eye_AP_length + finalData$Hyoid_rest_AP)^2 + finalData$D_head^2)[scanRowsNow] - finalData$D_head_backup[scanRowsNow]))
      # Resting ceratohyal half-width from the bar vs. W_head/2. Should be close to 0
  print(summary(sqrt(finalData$Ceratohyal_length^2 - finalData$Ceratohyal_AP_rest^2 - finalData$Ceratohyal_DV_rest^2)[scanRowsNow] - finalData$W_head[scanRowsNow]/2))

    # Eye-based vs. triangle videos: depression, width and volume should be in the same range
  print(tapply(finalData$Hyoid_Max, finalData$Hyoid_depth_source, summary))
  print(tapply(finalData$Buccal_width_expansion_Max, finalData$Hyoid_depth_source, summary))
  print(tapply(finalData$V_buccal_Max, finalData$Hyoid_depth_source, summary))

    # How often width had to be floored at 0, by depth source
  strikeFramesN <- finalData$endFrame - finalData$startFrame + 1
      # Share of videos with any floored frames
  print(tapply(finalData$Width_hypaxial_floor_frames > 0, finalData$Hyoid_depth_source, mean, na.rm = TRUE))
  print(tapply(finalData$Width_total_floor_frames > 0, finalData$Hyoid_depth_source, mean, na.rm = TRUE))
      # Fraction of strike frames floored
  print(tapply(finalData$Width_hypaxial_floor_frames / strikeFramesN, finalData$Hyoid_depth_source, summary))
  print(tapply(finalData$Width_total_floor_frames / strikeFramesN, finalData$Hyoid_depth_source, summary))

    # Are small gapes dominating U_flow_ff_predicted?   ### CHECK
      # Ratio far from 1, or a weak correlation, means the frames near Tstart/Tend (small A) are driving the full-strike mean
  print(summary(finalData$U_flow_ff_predicted_mean / finalData$U_flow_ff_predicted_mean_bigGape))
  print(corCheck(finalData$U_flow_ff_predicted_mean, finalData$U_flow_ff_predicted_mean_bigGape))     
        
        
        
        
        
  
        
#########################################
# Range of motion (per individual fish) #
#########################################
  # ROM = largest per-video _Max value across all of an individual's videos
  # Every video from the same individual gets the same ROM value
  # NA if none of that individual's videos have a value

  # Left side = ROM column, right side = per-video column it comes from
  romPairs <- c(Protrusion_ROM          = "Protrusion_Max",
                Mandible_ROM            = "Mandible_Max",
                Mandible_tip_disp_ROM   = "Mandible_tip_disp_Max",
                Mandible_depression_ROM = "Mandible_depression_Max",
                Hyoid_ROM               = "Hyoid_Max",
                Neurocranium_ROM        = "Neurocranium_Max",
                Hyoid_excursion_ROM_backup = "Hyoid_excursion_Max_backup")


  for(romNow in names(romPairs)) {

    maxColNow <- romPairs[romNow]

    for(idNow in unique(finalData$Unique.ID)) {

      rowsNow   <- which(finalData$Unique.ID == idNow)
      valuesNow <- finalData[[maxColNow]][rowsNow]

      if(any(!is.na(valuesNow))) {
        finalData[[romNow]][rowsNow] <- max(valuesNow, na.rm=TRUE)
      } else {
        finalData[[romNow]][rowsNow] <- NA_real_
      }
    }
  }


  # Checks
    # How many individuals have each ROM?
  print(sapply(names(romPairs), function(romNow) length(unique(finalData$Unique.ID[!is.na(finalData[[romNow]])]))))

    # Mandible_depression_ROM should equal Mandible_length * sin(Mandible_ROM) exactly, because each fish has one Mandible_length
  print(all.equal(finalData$Mandible_depression_ROM, finalData$Mandible_length * sin(finalData$Mandible_ROM)))
  
    # Maxilla sign check: most videos should have positive maxilla rotation at tmax (ventral tip swung anteriorly)
  print(sum(finalData$Maxilla_ang_tmax < 0, na.rm=TRUE))
  print(sum(finalData$Maxilla_ang_tmax > 0, na.rm=TRUE))

    # Measured vs. predicted maxilla angle. Should be positively correlated. If negative, oralOpenSign or the theta3 sign convention is flipped   ### CHECK
  print(corCheck(finalData$Maxilla_ang_tmax, finalData$Maxilla_ang_predicted_tmax))
  
  
    # Skipped analyses: how many videos per section and reason
  print(table(skippedVids$section))
  print(table(skippedVids$reason))

    # Videos with the most problems (best ones to look at first when checking landmarks)
  print(head(sort(table(skippedVids$vidName), decreasing = TRUE), 20))

    # Save the log so I can work through it alongside the videos
  write.csv(skippedVids, 'PATH/TO/skippedVids.csv', row.names = FALSE)   ### CHANGE path
  
  
  
  print(summary(finalData$strikeAxisTilt * 180/pi)) # shows how tilted strikes actually are.
  
  print(table(finalData$strikeAxisFromNasal)) # shows how often the axis fell back to the camera x-axis.
                          # If few strikes fall back and most tilts are only a few degrees, the projection is changing little compared with x-only ram, which is fine.
  
  
  
  
  
    
############################     CALCULATED VARIABLES      ############################
  
  
  
  ###########################
  # Kinetic synchronization #
  ###########################
      # SD of the three peak times divided by strike duration. Low = more synchronized
      # All three times are measured from Tstart (gape first passes 20% of its range): Time_hyoid, Ttpg, Time_cranial
      # sd() divides by n-1. With 3 values that's always sqrt(3/2) times the population SD, so it doesn't change any comparison
      # NA if any of the three times is missing
      # Kinetic_Synchronization_backup uses Time_hyoid_backup (Nasal-Hyoid distance) instead of the skull-frame Time_hyoid
        # Time_cranial needs the eye either way, so the backup only adds videos where the eye wasn't visible in the same frames as the hyoid
    
      # Kinetic_Synchronization
      peakTimes <- cbind(finalData$Time_hyoid, finalData$Ttpg, finalData$Time_cranial)
      finalData$Kinetic_Synchronization <- apply(peakTimes, 1, sd) / finalData$Ttotal
    
      # Kinetic_Synchronization_backup
      peakTimesBackup <- cbind(finalData$Time_hyoid_backup, finalData$Ttpg, finalData$Time_cranial)
      finalData$Kinetic_Synchronization_backup <- apply(peakTimesBackup, 1, sd) / finalData$Ttotal
    
    
      # Order of the peaks (s). SD loses the order, so these keep it
        # Positive = that peak came after peak gape
      finalData$Lag_hyoid_gape   <- finalData$Time_hyoid - finalData$Ttpg
      finalData$Lag_cranial_gape <- finalData$Time_cranial - finalData$Ttpg
    
    
      # Checks
        # How many videos have each version?
      sum(!is.na(finalData$Kinetic_Synchronization))
      sum(!is.na(finalData$Kinetic_Synchronization_backup))
    
        # Videos where both exist should be close. A low correlation means the two hyoid peak times disagree
      cor(finalData$Kinetic_Synchronization, finalData$Kinetic_Synchronization_backup, use = "complete.obs")
    
        # Typical order of the peaks
      summary(finalData$Lag_hyoid_gape)
      summary(finalData$Lag_cranial_gape)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
###########
# Save    #
###########
  # finalData: one row per video. Analyze in a separate script
  write.csv(finalData, 'PATH/TO/kinematicsAll_sleapOutput.csv', row.names = FALSE)   ### CHANGE path

  # finalTimeSeriesData is a list of 3D arrays, so it can't go in a CSV. Save it as an R object and read it back with readRDS()
  saveRDS(finalTimeSeriesData, 'PATH/TO/kinematicsTimeSeries.rds')   ### CHANGE path




