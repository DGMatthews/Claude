
#############################################################################################
# Check the smoothing in 3_PIV_Flow_Measurements.R                                          #
#############################################################################################
  # For a random sample of videos, smooths the raw velocity and pressure traces the same way script 3 does,
    # at the current lambda and at weaker / stronger multiples of it, then:
    # 1. Saves one PNG per video: raw points, the smoothed curves, and the rate of change (acceleration / dP/dt) underneath
    # 2. Prints two diagnostics per trace type and lambda multiple, as medians across the sampled videos
      # resid_lag1: lag-1 autocorrelation of the residuals (raw - smoothed)
        # Near 0 = residuals look like noise (good)
        # Clearly positive = smoothing is removing real signal (too strong). Clearly negative = the curve is following the noise (too weak)
      # peak_kept: smoothed peak / peak of a light 5-frame running median of the raw trace
        # Near 1 = the peak survives. Well below 1 = the peak is being flattened (too strong)
    # Too weak shows up in the plots: the rate of change (bottom panel) is jagged instead of a smooth curve
  # Uses the settings and functions at the top of script 3, so it always matches what script 3 does

rm(list=ls())

# Full path to script 3 (or just the file name if R's working directory is the folder with the scripts)
script3 <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/R/3_PIV_Flow_Measurements.R'   ### CHANGE
if(!file.exists(script3)) {
  stop("Can't find script 3 at '", script3, "'. Set script3 to its full path (working directory now: ", getwd(), ")")
}


##### Settings and functions from script 3 (everything above "Read the data")
  # Script 3 starts with rm(list=ls()), so the settings for this check are set after it
  src3 <- readLines(script3)
  eval(parse(text = src3[1:(grep("#  Read the data", src3, fixed = TRUE) - 2)]))


##### Settings for this check
  nCheck      <- 20                         ### CHANGE  number of videos to check
  lambdaMults <- c(0.2, 1, 5)               ### CHANGE  multiples of lambdaRefFlow to compare (1 = current)
  plotFolder  <- file.path(pathData, 'Flow_smoothing_check')   ### CHANGE  folder for the PNGs (made if needed). pathData comes from script 3
  set.seed(1)


##### Videos to check: kinematic videos with a filtered PIV result and a strike
  kin <- read.csv(kinFinalData, stringsAsFactors = FALSE)
  pivFiltered <- read_excel(pivFilteredFile)
  pivFiltered$videoKey <- pivKey(pivFiltered$video)

  # PIV name for each kinematic video. Kinematics output from before PIV.filename was saved: from the data sheet, as script 3 does
  if(!("PIV.filename" %in% names(kin))) {
    vidDataTemp <- read.csv(dataSheet, stringsAsFactors = FALSE)
    vidData     <- vidDataTemp[which(vidDataTemp$Generation == 'F2' & vidDataTemp$Event != '-'), ]   # same filter as the kinematics script
    kin$PIV.filename <- as.character(vidData$PIV.filename[kin$equivVidRow])
    print("Kinematics output has no PIV.filename column (made before 25 Oct): PIV names taken from the data sheet. Rerun script 2 to update it")
  }

  kin$pivRow <- match(pivKey(kin$PIV.filename), pivFiltered$videoKey)
  candidates <- which(!is.na(kin$pivRow) & !is.na(kin$Tstart) & !is.na(kin$nFrames))
  checkRows  <- candidates[sample(length(candidates), min(nCheck, length(candidates)))]

  dir.create(plotFolder, showWarnings = FALSE)
  diag <- data.frame()


##### Smooth and plot
  for(i in checkRows) {

    pivVideo  <- pivFiltered$video[kin$pivRow[i]]
    nNow      <- kin$nFrames[i]
    frameRate <- kin$frameRate[i]
    timeNow   <- (0:(nNow - 1)) / frameRate

    traces <- list()
    fileTraces <- file.path(pathRawTraces, paste0(pivVideo, suffixPeakTraces))
    if(file.exists(fileTraces)) {
      tr <- read.csv(fileTraces)
      traces$staticRaw_velocity <- list(y = toKinFrames(tr$velocity_time_ms, tr$velocity_m_per_s * 1000, nNow, frameRate), madFloor = madFloorVel, unit = "mm/s")
      traces$pressure           <- list(y = toKinFrames(tr$pressure_time_ms, tr$pressure_Pa, nNow, frameRate), madFloor = madFloorPressure, unit = "Pa")
    }
    fileFront <- file.path(pathMouthFront, paste0(pivVideo, suffixMouthFront))
    if(file.exists(fileFront)) {
      fr <- read.csv(fileFront)
      traces$mouthFront_velocity <- list(y = toKinFrames(fr$time_ms, fr$mouth_front_velocity_raw_m_per_s * 1000, nNow, frameRate), madFloor = madFloorVel, unit = "mm/s")
    }
    traces <- traces[sapply(traces, function(x) sum(!is.na(x$y)) >= 10)]
    if(length(traces) == 0) {
      next
    }

    png(file.path(plotFolder, paste0(kin$vidName[i], ".png")), width = 400 * length(traces), height = 700)
    par(mfcol = c(2, length(traces)), mar = c(4, 4, 3, 1))
    colsNow <- c("blue", "black", "red")[seq_along(lambdaMults)]

    for(trNm in names(traces)) {
      yRaw <- traces[[trNm]]$y
      fits <- lapply(lambdaMults, function(m) smoothTrace(yRaw, timeNow, traces[[trNm]]$madFloor, lambdaRefFlow * m))

      # Value panel
      plot(timeNow, yRaw, pch = 16, cex = 0.5, col = "grey50", main = paste(trNm), xlab = "Time (s)", ylab = traces[[trNm]]$unit)
      for(k in seq_along(fits)) lines(timeNow, fits[[k]]$y, col = colsNow[k], lwd = 1.5)
      abline(v = c(kin$Tstart[i], kin$Tend[i]), lty = 3)
      abline(v = kin$tmax[i], col = "red", lty = 2)
      abline(v = kin$tmaxVel[i], col = "blue", lty = 2)
      legend("topright", legend = paste0("lambda x", lambdaMults), col = colsNow, lwd = 1.5, bty = "n", cex = 0.8)

      # Rate-of-change panel (acceleration for velocity, dP/dt for pressure)
      dyRange <- range(unlist(lapply(fits, function(f) f$dy)), na.rm = TRUE)
      plot(NA, xlim = range(timeNow), ylim = dyRange, main = paste0(trNm, ": rate of change"), xlab = "Time (s)", ylab = paste0(traces[[trNm]]$unit, "/s"))
      for(k in seq_along(fits)) lines(timeNow, fits[[k]]$dy, col = colsNow[k], lwd = 1.5)
      abline(h = 0, col = "grey")

      # Diagnostics
      runMed <- stats::runmed(ifelse(is.na(yRaw), median(yRaw, na.rm = TRUE), yRaw), 5)
      peakRaw <- if(grepl("pressure", trNm)) min(runMed) else max(runMed)
      for(k in seq_along(fits)) {
        res <- yRaw - fits[[k]]$y
        res <- res[!is.na(res)]
        peakSm <- if(grepl("pressure", trNm)) min(fits[[k]]$y, na.rm = TRUE) else max(fits[[k]]$y, na.rm = TRUE)
        diag <- rbind(diag, data.frame(trace = trNm, lambdaMult = lambdaMults[k],
                                       resid_lag1 = if(length(res) > 3) cor(head(res, -1), tail(res, -1)) else NA,
                                       peak_kept  = peakSm / peakRaw))
      }
    }
    dev.off()
  }


##### Summary across the sampled videos
  print(aggregate(cbind(resid_lag1, peak_kept) ~ trace + lambdaMult, data = diag, FUN = median))
  print(paste0("Plots saved in: ", normalizePath(plotFolder)))
