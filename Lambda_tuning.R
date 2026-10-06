
#############################################################################################
# Lambda tuning for the kinematics scripts                                                   #
#############################################################################################
  # Fits a few test videos over a range of smoothing values and saves plots and tables,
    # so lambdaRef, lambdaRefDerived and lambdaRefRam can be chosen and pasted into
    # ACxTRC_Kinematic_Measurements_Sept2026.R and ACxTRC_Kinematic_Measurements_JawScansOnly.R
  # Uses the same Hampel filter, gap filling and span scaling as the kinematics scripts
  # Works in pixels and frames. smooth.spline's lambda doesn't depend on the units of x or y,
    # so the chosen values apply directly to the kinematics scripts (which work in mm and seconds)

  # Three settings, tuned in this order:
    # 1. lambdaRef:        raw landmark coordinates (first pass)
    # 2. lambdaRefDerived: derived series (gape here; also protrusion, angles, hyoid in the main scripts)
    # 3. lambdaRefRam:     mouth position for ram speed and acceleration
  # Each stage uses the CURRENT value of the earlier stages, so after changing lambdaRef, rerun before judging the others

  # What to look for (per stage, across the multipliers):
    # Residual SD (raw - smoothed) close to the tracking noise. Much lower = still fitting noise, much higher = removing real motion
    # Lag-1 autocorrelation of the residuals near 0
      # Clearly negative = the fit is still chasing frame-to-frame noise (undersmoothing)
      # Clearly positive = smoothing is removing real motion (oversmoothing)
      # Only meaningful for stage 1. Stages 2 and 3 start from already-smoothed coordinates, so their residuals
        # are correlated by construction: judge those stages by the peaks, velocities and accelerations instead
    # Peak values that stay within a few % of the lighter fits (too much smoothing flattens peaks)
    # Velocities and accelerations that stop changing much between neighbouring multipliers (the plateau)
    # Ram mean speed from end positions vs. the frame mean of the velocity should agree within a few %
  # GCV: smooth.spline's own choice of lambda, converted to the reference span. Tends to undersmooth when
    # tracking errors are correlated between frames, so treat it as a lower bound

#This is how you install the rhdf5 package
#install.packages("BiocManager")
#BiocManager::install("rhdf5")

require(rhdf5)
require(imputeTS)

rm(list=ls())


##############
#  Settings  #
##############

  # Folder with the SLEAP .h5 files
  path <- 'N:/TRCxAC/Hybrid feeding project/HSV/KINE/SLEAP landmarks/F2'

  # Where to save the plots and tables
  pathOutput <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/R/Grant Update Results/Lambda_tuning'

  # Test videos: numbers (positions in the folder's file list) or names (without ".analysis.h5")
    # Pick ~8-10 uncropped videos: short and long, fast and slow, clean and messy tracking, both facing directions
  testVids <- seq(from=1, to=301, by=25)   ### CHANGE

  # Current values (copy from the kinematics script)
  lambdaRef        <- 5e-6                ### CHANGE
  lambdaRefDerived <- lambdaRef/5         ### CHANGE
  lambdaRefRam     <- lambdaRefDerived    ### CHANGE
  spanRef          <- 73

  # Multipliers to try on each lambda
  lambdaMults <- c(0.1, 0.3, 1, 3, 10)   ### CHANGE  1 = the current value

  # Landmarks used for the coordinate stage (OpercTab isn't used in the kinematics)
  expectedNames <- c("UJ", "LJ", "Nasal", "VentMax", "Eye", "Hyoid", "OpercTab")
  coordLandmarks <- c("UJ", "LJ", "Nasal", "VentMax", "Eye", "Hyoid")


#####################
#  Define functions #
#####################

    # Hampel filter (same as the kinematics scripts): frames more than threshold * MAD from the rolling median
    applyHampel <- function(x, halfWindow=3, threshold=3, madFloor=1) {
      flaggedFrames <- c()
      for(j in seq_along(x)){
        if(is.na(x[j])) {
          next
        }
        minFrameNow <- max(1, j - halfWindow)
        maxFrameNow <- min(length(x), j + halfWindow)
        pointsNow <- x[minFrameNow:maxFrameNow]
        medianNow <- median(pointsNow, na.rm=TRUE)
        MADnow <- mad(pointsNow, na.rm = TRUE)
        if(MADnow<madFloor) {MADnow=madFloor}
        if(abs(x[j] - medianNow) > threshold * MADnow) {
          flaggedFrames <- c(flaggedFrames,j)
        }
      }
      return(flaggedFrames)
    }


    # Smoothing spline with the span-scaled lambda, as in the kinematics scripts
      # Returns the fit, or NULL if there are fewer than 10 frames
    fitScaled <- function(frames, y, lambdaRefNow) {
      good <- which(!is.na(y))
      if(length(good) < 10) {
        return(NULL)
      }
      spanNow <- max(frames[good]) - min(frames[good])
      smooth.spline(x = frames[good], y = y[good], lambda = lambdaRefNow * (spanRef / spanNow)^3)
    }


    # GCV lambda, converted to the reference span so it's comparable with lambdaRef
    gcvRef <- function(frames, y) {
      good <- which(!is.na(y))
      if(length(good) < 10) {
        return(NA_real_)
      }
      spanNow <- max(frames[good]) - min(frames[good])
      fitNow <- smooth.spline(x = frames[good], y = y[good])
      fitNow$lambda * (spanNow / spanRef)^3
    }


    # Lag-1 autocorrelation, ignoring NAs
    lag1 <- function(r) {
      ok <- which(!is.na(r[-length(r)]) & !is.na(r[-1]))
      if(length(ok) < 5) {
        return(NA_real_)
      }
      cor(r[-length(r)][ok], r[-1][ok])
    }


    # Hampel + gap filling for one video (same settings as the kinematics scripts)
    cleanTracks <- function(dataNow) {
      for(node in 1:dim(dataNow)[2]){
        xFlags <- applyHampel(dataNow[,node,1,1])
        yFlags <- applyHampel(dataNow[,node,2,1])
        dataNow[union(xFlags, yFlags),node,,1] <- NA
      }
      for(node in 1:dim(dataNow)[2]){
        xPointsNow <- dataNow[,node,1,1]
        yPointsNow <- dataNow[,node,2,1]
        realFrames <- which(!is.na(xPointsNow) & !is.na(yPointsNow))
        if(length(realFrames)>=10) {
          firstReal <- min(realFrames)
          lastReal <- max(realFrames)
          dataNow[firstReal:lastReal,node,1,1] <- na_interpolation(xPointsNow[firstReal:lastReal], option="linear", maxgap=3)
          dataNow[firstReal:lastReal,node,2,1] <- na_interpolation(yPointsNow[firstReal:lastReal], option="linear", maxgap=3)
        }
      }
      dataNow
    }


    # Smoothed coordinates for one landmark (x and y) at a given lambdaRef, predicted only at frames with data
    smoothLandmark <- function(dataNow, node, lambdaRefNow) {
      nF <- dim(dataNow)[1]
      out <- matrix(NA_real_, nF, 2)
      for(xy in 1:2) {
        fitNow <- fitScaled(1:nF, dataNow[,node,xy,1], lambdaRefNow)
        if(!is.null(fitNow)) {
          good <- which(!is.na(dataNow[,node,xy,1]))
          out[good, xy] <- predict(fitNow, good)$y
        }
      }
      out
    }


###################
#  Find the files #
###################

  fileNames  <- list.files(path, pattern="\\.h5$", all.files=TRUE, no..=TRUE)
  fileNames2 <- gsub("\\.analysis\\.h5$", "", fileNames)

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

  dir.create(pathOutput, showWarnings = FALSE, recursive = TRUE)

  multCols <- hcl.colors(length(lambdaMults), "Viridis")
  results  <- data.frame()
  gcvTable <- data.frame()


#########################
#  Loop over test videos #
#########################

for(v in testKeep) {

  vidName <- fileNames2[v]
  print(paste0("Tuning on: ", vidName))

  landmarkNames <- h5read(file=file.path(path, fileNames[v]), name='/node_names')
  if(!identical(as.vector(landmarkNames), expectedNames)){
    stop("The landmark names are different from expected in: ", fileNames[v])
  }

  dataNow <- h5read(file=file.path(path, fileNames[v]), name='/tracks', index=list(NULL,NULL,NULL,NULL))
  dataNow <- cleanTracks(dataNow)
  nF <- dim(dataNow)[1]
  frames <- 1:nF


  ##### Tracking noise (px, per coordinate)
    # The Eye-Nasal distance is fixed on a rigid skull, so its frame-to-frame scatter is tracking noise
    # The distance carries the noise of two landmarks (sqrt(2) * per-coordinate noise), and frame-to-frame differences
      # double the variance again (another sqrt(2)), so per-coordinate noise = SD of the differences / 2
    # MAD instead of SD, so real motion (and eye sliding) doesn't inflate it
  eyeNasal <- sqrt((dataNow[,3,1,1] - dataNow[,5,1,1])^2 + (dataNow[,3,2,1] - dataNow[,5,2,1])^2)
  noisePx  <- if(sum(!is.na(eyeNasal)) >= 10) mad(diff(eyeNasal), na.rm=TRUE) / 2 else NA_real_


  ##### GCV lambdas for each landmark coordinate
  for(lmNow in coordLandmarks) {
    node <- match(lmNow, expectedNames)
    for(xy in 1:2) {
      gcvTable <- rbind(gcvTable, data.frame(video = vidName, landmark = lmNow, coord = c("x","y")[xy],
                                             lambdaRef_GCV = gcvRef(frames, dataNow[,node,xy,1])))
    }
  }


  ##### Stage 1: lambdaRef (raw coordinates)
    # Residual SD and lag-1 autocorrelation, pooled over the landmarks and both coordinates
  coordFits <- list()

  for(m in seq_along(lambdaMults)) {

    lambdaNow <- lambdaRef * lambdaMults[m]
    residsAll <- c()
    lagsAll   <- c()

    for(lmNow in coordLandmarks) {
      node <- match(lmNow, expectedNames)
      smoothNow <- smoothLandmark(dataNow, node, lambdaNow)
      if(lmNow == "LJ") {
        coordFits[[m]] <- smoothNow
      }
      for(xy in 1:2) {
        residNow <- dataNow[,node,xy,1] - smoothNow[,xy]
        residsAll <- c(residsAll, residNow)
        lagsAll   <- c(lagsAll, lag1(residNow))
      }
    }

    results <- rbind(results, data.frame(video = vidName, stage = "lambdaRef", mult = lambdaMults[m], lambda = lambdaNow,
                                         noisePx = noisePx, residSD = sd(residsAll, na.rm=TRUE), residLag1 = mean(lagsAll, na.rm=TRUE),
                                         peakR = NA, maxRvel = NA, maxRamAcc = NA, ramMeanDiffPct = NA))
  }


  ##### Coordinates at the CURRENT lambdaRef, for the derived and ram stages
  ujNow <- smoothLandmark(dataNow, 1, lambdaRef)
  ljNow <- smoothLandmark(dataNow, 2, lambdaRef)

  R_t <- sqrt((ljNow[,1] - ujNow[,1])^2 + (ljNow[,2] - ujNow[,2])^2) / 2   # gape radius (px)


  ##### Strike window, same threshold rule as the kinematics scripts (20% of the gape range)
  startPos <- NA
  endPos   <- NA
  if(any(!is.na(R_t))) {
    peakFrame <- which.max(R_t)
    R_min <- min(R_t[1:peakFrame], na.rm = TRUE)
    startThreshold <- R_min + 0.2 * (R_t[peakFrame] - R_min)
    startPos <- which(R_t > startThreshold)[1]
    endPos   <- which(R_t[peakFrame:nF] < startThreshold)[1] + peakFrame - 1
  }
  if(is.na(startPos) || is.na(endPos)) {
    print(paste0("No strike window found, skipping the derived and ram stages for: ", vidName))
    next
  }
  strikeFrames <- startPos:endPos


  ##### Stage 2: lambdaRefDerived (gape radius)
    # Peak gape and max opening velocity (px/frame) within the strike
  RFits <- list()
  RVels <- list()

  for(m in seq_along(lambdaMults)) {

    lambdaNow <- lambdaRefDerived * lambdaMults[m]
    fitNow <- fitScaled(frames, R_t, lambdaNow)
    if(is.null(fitNow)) next

    good <- which(!is.na(R_t))
    RFits[[m]] <- rep(NA_real_, nF)
    RVels[[m]] <- rep(NA_real_, nF)
    RFits[[m]][good] <- predict(fitNow, good)$y
    RVels[[m]][good] <- predict(fitNow, good, deriv = 1)$y

    results <- rbind(results, data.frame(video = vidName, stage = "lambdaRefDerived", mult = lambdaMults[m], lambda = lambdaNow,
                                         noisePx = noisePx, residSD = sd(R_t - RFits[[m]], na.rm=TRUE), residLag1 = lag1(R_t - RFits[[m]]),
                                         peakR = max(RFits[[m]][strikeFrames], na.rm=TRUE), maxRvel = max(RVels[[m]][strikeFrames], na.rm=TRUE),
                                         maxRamAcc = NA, ramMeanDiffPct = NA))
  }


  ##### Stage 3: lambdaRefRam (mouth position along the camera x-axis)
    # Camera x is enough for tuning (the kinematics scripts project onto the strike axis, which is nearly horizontal)
    # Max |acceleration| within the strike (px/frame^2), and the ram mean check:
      # (x at Tend - x at Tstart) / strike frames vs. the frame mean of the velocity
  xMouth <- (ujNow[,1] + ljNow[,1]) / 2
  ramVels <- list()
  ramAccs <- list()

  for(m in seq_along(lambdaMults)) {

    lambdaNow <- lambdaRefRam * lambdaMults[m]
    fitNow <- fitScaled(frames, xMouth, lambdaNow)
    if(is.null(fitNow)) next

    good <- which(!is.na(xMouth))
    ramVels[[m]] <- rep(NA_real_, nF)
    ramAccs[[m]] <- rep(NA_real_, nF)
    ramVels[[m]][good] <- predict(fitNow, good, deriv = 1)$y
    ramAccs[[m]][good] <- predict(fitNow, good, deriv = 2)$y

    meanFromEnds   <- (predict(fitNow, endPos)$y - predict(fitNow, startPos)$y) / (endPos - startPos)
    meanFromFrames <- mean(ramVels[[m]][strikeFrames], na.rm=TRUE)

    results <- rbind(results, data.frame(video = vidName, stage = "lambdaRefRam", mult = lambdaMults[m], lambda = lambdaNow,
                                         noisePx = noisePx, residSD = sd(xMouth - predict(fitNow, frames)$y, na.rm=TRUE),
                                         residLag1 = lag1(xMouth - predict(fitNow, frames)$y),
                                         peakR = NA, maxRvel = NA, maxRamAcc = max(abs(ramAccs[[m]][strikeFrames]), na.rm=TRUE),
                                         ramMeanDiffPct = 100 * (meanFromFrames - meanFromEnds) / abs(meanFromEnds)))
  }


  ##### Plot for this video
    # Row 1: LJ y raw vs. fits (lambdaRef), stage 1 residual SD and autocorrelation, gape fits, gape velocity
    # Row 2: zoom on the strike for LJ y, mouth x-velocity, mouth acceleration, legend
  png(file.path(pathOutput, paste0("lambda_", vidName, ".png")), width = 2000, height = 1000, res = 110)
  par(mfrow = c(2, 4), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 2, 0))

  # LJ y, whole video
  plot(frames, dataNow[,2,2,1], pch = 16, cex = 0.4, col = "grey60", xlab = "Frame", ylab = "LJ y (px)", main = "Stage 1: LJ y, raw vs fits")
  for(m in seq_along(coordFits)) lines(frames, coordFits[[m]][,2], col = multCols[m], lwd = 1.5)
  abline(v = c(startPos, endPos), lty = 3, col = "grey40")

  # Stage 1 summary (wider right margin for the second axis)
  par(mar = c(4, 4.5, 2.5, 4.5))
  s1 <- results[results$video == vidName & results$stage == "lambdaRef", ]
  plot(s1$mult, s1$residSD, log = "x", type = "b", pch = 16, xlab = "Multiplier on lambdaRef", ylab = "Residual SD (px)",
       main = "Stage 1: residual SD (black) vs noise (red)", ylim = range(c(s1$residSD, noisePx, 0), na.rm=TRUE))
  abline(h = noisePx, col = "red", lty = 2)
  par(new = TRUE)
  plot(s1$mult, s1$residLag1, log = "x", type = "b", pch = 1, lty = 3, col = "blue", axes = FALSE, xlab = "", ylab = "", ylim = c(-1, 1))
  axis(4, col = "blue", col.axis = "blue")
  abline(h = 0, col = "blue", lty = 3)
  mtext("Residual lag-1 autocorrelation (blue)", side = 4, line = 2.5, col = "blue", cex = 0.7)
  par(mar = c(4, 4.5, 2.5, 1))

  # Gape radius fits
  plot(frames, R_t, pch = 16, cex = 0.4, col = "grey60", xlim = range(strikeFrames) + c(-10, 10), xlab = "Frame", ylab = "R (px)", main = "Stage 2: gape radius fits")
  for(m in seq_along(RFits)) if(!is.null(RFits[[m]])) lines(frames, RFits[[m]], col = multCols[m], lwd = 1.5)

  # Gape velocity
  velRange <- range(unlist(lapply(RVels, function(z) z[strikeFrames])), na.rm=TRUE)
  plot(NA, xlim = range(strikeFrames), ylim = velRange, xlab = "Frame", ylab = "dR/dt (px/frame)", main = "Stage 2: gape velocity")
  for(m in seq_along(RVels)) if(!is.null(RVels[[m]])) lines(frames, RVels[[m]], col = multCols[m], lwd = 1.5)
  abline(h = 0, col = "grey80")

  # LJ y, strike zoom
  plot(frames, dataNow[,2,2,1], pch = 16, cex = 0.6, col = "grey60", xlim = range(strikeFrames) + c(-10, 10),
       ylim = range(dataNow[strikeFrames,2,2,1], na.rm=TRUE), xlab = "Frame", ylab = "LJ y (px)", main = "Stage 1: LJ y, strike zoom")
  for(m in seq_along(coordFits)) lines(frames, coordFits[[m]][,2], col = multCols[m], lwd = 1.5)

  # Mouth velocity
  velRange <- range(unlist(lapply(ramVels, function(z) z[strikeFrames])), na.rm=TRUE)
  plot(NA, xlim = range(strikeFrames) + c(-10, 10), ylim = velRange, xlab = "Frame", ylab = "Mouth x velocity (px/frame)", main = "Stage 3: ram speed (camera x)")
  for(m in seq_along(ramVels)) if(!is.null(ramVels[[m]])) lines(frames, ramVels[[m]], col = multCols[m], lwd = 1.5)
  abline(v = c(startPos, endPos), lty = 3, col = "grey40")

  # Mouth acceleration
  accRange <- range(unlist(lapply(ramAccs, function(z) z[strikeFrames])), na.rm=TRUE)
  plot(NA, xlim = range(strikeFrames) + c(-10, 10), ylim = accRange, xlab = "Frame", ylab = "Mouth x acceleration (px/frame^2)", main = "Stage 3: ram acceleration")
  for(m in seq_along(ramAccs)) if(!is.null(ramAccs[[m]])) lines(frames, ramAccs[[m]], col = multCols[m], lwd = 1.5)
  abline(v = c(startPos, endPos), lty = 3, col = "grey40")
  abline(h = 0, col = "grey80")

  # Legend
  plot.new()
  legend("center", legend = paste0("x ", lambdaMults, ifelse(lambdaMults == 1, "  (current)", "")), col = multCols, lwd = 3,
         title = "Multiplier on each lambda", bty = "n", cex = 1.3)

  mtext(paste0(vidName, "     |     grey dotted = strike window (Tstart, Tend)"), outer = TRUE, cex = 1.1)
  dev.off()
}


#############
#  Summary  #
#############

  write.csv(results,  file.path(pathOutput, "lambda_tuning_results.csv"), row.names = FALSE)
  write.csv(gcvTable, file.path(pathOutput, "lambda_tuning_GCV.csv"), row.names = FALSE)

  # Each value relative to the same video's value at the current lambda (mult = 1), then the median across videos
    # Look for the range of multipliers where peaks and velocities stay near 1 (the plateau)
  relCols <- c("residSD", "residLag1", "peakR", "maxRvel", "maxRamAcc")
  relTable <- results
  for(colNow in setdiff(relCols, "residLag1")) {
    baseNow <- ave(ifelse(results$mult == 1, results[[colNow]], NA), results$video, results$stage, FUN = function(z) z[!is.na(z)][1])
    relTable[[colNow]] <- results[[colNow]] / baseNow
  }

  for(stageNow in c("lambdaRef", "lambdaRefDerived", "lambdaRefRam")) {
    print(paste0("===== ", stageNow, ": median across videos (values relative to the current lambda, except residLag1 and ramMeanDiffPct) ====="))
    tabNow <- relTable[relTable$stage == stageNow, ]
    print(aggregate(cbind(residSD, residLag1, peakR, maxRvel, maxRamAcc, ramMeanDiffPct) ~ mult, data = tabNow,
                    FUN = function(z) round(median(z, na.rm=TRUE), 3), na.action = na.pass))
  }

  # Tracking noise vs. stage 1 residual SD at the current lambda (px). These should be similar
  print("Tracking noise (Eye-Nasal) vs. stage 1 residual SD at the current lambdaRef (px):")
  print(results[results$stage == "lambdaRef" & results$mult == 1, c("video", "noisePx", "residSD")])

  # GCV suggestion for lambdaRef (median over landmarks and videos). A lower bound, not a target
  print(paste0("GCV lambdaRef (median): ", signif(median(gcvTable$lambdaRef_GCV, na.rm=TRUE), 3),
               "   current lambdaRef: ", lambdaRef))

  print(paste0("Plots and tables saved to: ", pathOutput))
