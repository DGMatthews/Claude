
#############################################################################################
# Jaw width (W_jaw) from 3D Slicer landmarks                                                 #
#############################################################################################
  # Reads the mandible and premaxilla landmark files (.mrk.json from 3D Slicer) for each F2 individual,
    # calculates the width of the tooth row, and saves one row per fish for the kinematics scripts
  # Midline (sagittal) plane = the plane through 3 midline landmarks:
    # Mandible landmark 1, Mandible landmark 2, Premaxilla landmark 1
  # Half-width = shortest (perpendicular) distance from Mandible landmark 3 to that plane
  # W_jaw = 2 * half-width (assumes the left and right sides of the jaw are symmetric)
  # All values in mm (the Slicer coordinates are already in mm)

  # Output: jaw_widths_F2.csv, one row per fish, with columns:
    # Family, Animal.ID: used to match the kinematic videos (ACxTRC_Kinematic_Measurements_JawScansOnly.R)
    # W_jaw (mm), halfWidth (mm)
    # planeSine: how far the 3 midline landmarks are from lying on a line (1 = right angle, 0 = on a line). Low = unreliable plane
    # specimen, mandibleFile, premaxFile: so any odd value can be traced back to its files

#This script needs the jsonlite package to read the Slicer files
#install.packages("jsonlite")

require(jsonlite)

rm(list=ls())


##############
#  Settings  #
##############

  # Folders with the landmark files
  pathMandible <- 'D:/CT scans/Analysis/F2 Integration/AM-Jaw complex/Landmarks/Mandible/Combined'
  pathPremax   <- 'D:/CT scans/Analysis/F2 Integration/AM-Jaw complex/Landmarks/Premax/Combined'

  # Where to save the results
  pathOutput   <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/R/Grant Update Results'

  # Landmark numbers
  mandMidline1 <- 1   # Mandible landmark on the midline
  mandMidline2 <- 2   # Mandible landmark on the midline
  pmxMidline   <- 1   # Premaxilla landmark on the midline
  mandLateral  <- 3   # Mandible landmark whose distance to the midline plane is the half-width

  # The filenames hold the family as a number (e.g. "1.2" in F2_1.2_22.02_Mandible.mrk.json)
    # The kinematics script matches on the Family column of the feeding data sheet, so the family has to be written the same way
  familyPrefix <- "F"   ### CHECK  "F" turns "1.2" into "F1.2". Set to "" if the data sheet writes families as "1.2"

  # If the 3 midline landmarks are almost on a line, the plane through them can tilt freely and the width is meaningless
    # planeSine = sine of the angle at Mandible landmark 1 between the other two midline landmarks
  minPlaneSine <- 0.2   ### CHANGE  stop if any fish is below this (0.2 is about 12 degrees)


#####################
#  Define functions #
#####################

    # Read one landmark from a Slicer .mrk.json file
      # Points are found by their id, not their order in the file, so extra or reordered points don't matter
      # Some files contain leftover points copied from other specimens (e.g. labels "OC_10-2") that reuse the same ids
        # If an id appears more than once, only the point whose label starts with this file's specimen name is used
    readLandmark <- function(filePath, landmarkNum, specimenLabel) {

      controlPoints <- fromJSON(filePath)$markups$controlPoints[[1]]

      rowsNow <- which(controlPoints$id == as.character(landmarkNum))

      if(length(rowsNow) > 1) {
        rowsNow <- rowsNow[startsWith(controlPoints$label[rowsNow], specimenLabel)]
      }

      if(length(rowsNow) == 0) {
        stop("Landmark ", landmarkNum, " is missing (or only belongs to another specimen) in: ", filePath)
      }
      if(length(rowsNow) > 1) {
        stop("Landmark ", landmarkNum, " appears ", length(rowsNow), " times with this specimen's label in: ", filePath)
      }
      if(controlPoints$positionStatus[rowsNow] != "defined") {
        stop("Landmark ", landmarkNum, " is not placed (positionStatus = '", controlPoints$positionStatus[rowsNow], "') in: ", filePath)
      }

      unlist(controlPoints$position[rowsNow])
    }


    # Cross product of two 3D vectors
    crossProduct <- function(u, v) {
      c(u[2]*v[3] - u[3]*v[2],
        u[3]*v[1] - u[1]*v[3],
        u[1]*v[2] - u[2]*v[1])
    }


    # Specimen name from a filename, e.g. "F2_1.2_22.02_Mandible.mrk.json" -> "F2_1.2_22.02"
    specimenFromFile <- function(fileName, boneName) {
      sub(paste0("_", boneName, "(\\.mrk)?\\.json$"), "", fileName)
    }


###################
#  Find the files #
###################

  mandFiles <- list.files(pathMandible, pattern = "\\.json$")
  pmxFiles  <- list.files(pathPremax,   pattern = "\\.json$")

  if(length(mandFiles) == 0) {
    stop("No .json files in the mandible folder: ", pathMandible)
  }
  if(length(pmxFiles) == 0) {
    stop("No .json files in the premaxilla folder: ", pathPremax)
  }

  # Every file must end in _Mandible.mrk.json / _Premax.mrk.json, otherwise the specimen name can't be read reliably
  badMand <- mandFiles[!grepl("_Mandible(\\.mrk)?\\.json$", mandFiles)]
  badPmx  <- pmxFiles[!grepl("_Premax(\\.mrk)?\\.json$", pmxFiles)]
  if(length(badMand) > 0 || length(badPmx) > 0) {
    stop("These files don't end in _Mandible.mrk.json or _Premax.mrk.json: ", paste(c(badMand, badPmx), collapse = ", "))
  }

  mandSpecimens <- specimenFromFile(mandFiles, "Mandible")
  pmxSpecimens  <- specimenFromFile(pmxFiles, "Premax")

  if(any(duplicated(mandSpecimens)) || any(duplicated(pmxSpecimens))) {
    stop("More than one file for the same specimen: ", paste(unique(c(mandSpecimens[duplicated(mandSpecimens)], pmxSpecimens[duplicated(pmxSpecimens)])), collapse = ", "))
  }

  # Only F2 specimens (filenames starting with "F2_"). Anything else in the folders is listed and left out
  notF2 <- union(mandSpecimens[!startsWith(mandSpecimens, "F2_")], pmxSpecimens[!startsWith(pmxSpecimens, "F2_")])
  if(length(notF2) > 0) {
    print("Not F2, left out:")
    print(notF2)
  }

  # Specimens with only one of the two bones can't be measured
  onlyMand <- setdiff(mandSpecimens, pmxSpecimens)
  onlyPmx  <- setdiff(pmxSpecimens, mandSpecimens)
  print("Mandible landmarks but no premaxilla landmarks (left out):")
  print(onlyMand)
  print("Premaxilla landmarks but no mandible landmarks (left out):")
  print(onlyPmx)

  specimens <- sort(intersect(mandSpecimens, pmxSpecimens))
  specimens <- specimens[startsWith(specimens, "F2_")]
  nSpec <- length(specimens)

  print(paste0("Measuring ", nSpec, " F2 specimens"))


##################
#  Measure width #
##################

  jawWidths <- data.frame(specimen     = specimens,
                          Family       = rep(NA_character_, nSpec),
                          Animal.ID    = rep(NA_character_, nSpec),
                          W_jaw        = rep(NA_real_, nSpec),
                          halfWidth    = rep(NA_real_, nSpec),
                          planeSine    = rep(NA_real_, nSpec),
                          mandibleFile = mandFiles[match(specimens, mandSpecimens)],
                          premaxFile   = pmxFiles[match(specimens, pmxSpecimens)])

  for(i in 1:nSpec) {

    ##### Family and individual from the specimen name, e.g. "F2_1.2_22.02" -> family "1.2", individual "22.02"
      nameParts <- strsplit(specimens[i], "_")[[1]]

      if(length(nameParts) != 3) {
        stop("Can't read family and individual from '", specimens[i], "'. Expected F2_family_individual")
      }

      jawWidths$Family[i]    <- paste0(familyPrefix, nameParts[2])
      jawWidths$Animal.ID[i] <- nameParts[3]


    ##### Landmarks (mm)

      mandFileNow <- file.path(pathMandible, jawWidths$mandibleFile[i])
      pmxFileNow  <- file.path(pathPremax,   jawWidths$premaxFile[i])

      midA    <- readLandmark(mandFileNow, mandMidline1, specimenLabel = paste0(specimens[i], "_Mandible"))
      midB    <- readLandmark(mandFileNow, mandMidline2, specimenLabel = paste0(specimens[i], "_Mandible"))
      midC    <- readLandmark(pmxFileNow,  pmxMidline,   specimenLabel = paste0(specimens[i], "_Premax"))
      lateral <- readLandmark(mandFileNow, mandLateral,  specimenLabel = paste0(specimens[i], "_Mandible"))


    ##### Midline plane through the 3 midline landmarks
      # Its normal (perpendicular) direction is the cross product of two vectors lying in the plane

      inPlane1 <- midB - midA
      inPlane2 <- midC - midA
      normalNow <- crossProduct(inPlane1, inPlane2)

      # planeSine = |normal| / (|inPlane1| * |inPlane2|) = sine of the angle between the two in-plane vectors
      jawWidths$planeSine[i] <- sqrt(sum(normalNow^2)) / (sqrt(sum(inPlane1^2)) * sqrt(sum(inPlane2^2)))

      if(jawWidths$planeSine[i] < minPlaneSine) {
        stop("The 3 midline landmarks are almost on a line (planeSine = ", round(jawWidths$planeSine[i], 3),
             "), so the midline plane isn't defined. Check the landmarks for: ", specimens[i])
      }

      unitNormal <- normalNow / sqrt(sum(normalNow^2))


    ##### Half-width = perpendicular distance from the lateral landmark to the plane

      jawWidths$halfWidth[i] <- abs(sum((lateral - midA) * unitNormal))
      jawWidths$W_jaw[i]     <- 2 * jawWidths$halfWidth[i]
  }


############
#  Checks  #
############

  # Range of widths (mm)
  print(summary(jawWidths$W_jaw))

  # How well-defined the midline planes are. Values near minPlaneSine deserve a look
  print(summary(jawWidths$planeSine))

  # Most extreme widths, to check against the scans
    # Far from the median (more than 3 MADs) = worth opening in Slicer
  widthMedian <- median(jawWidths$W_jaw)
  widthMAD    <- mad(jawWidths$W_jaw)
  print("Widths more than 3 MADs from the median:")
  print(jawWidths[abs(jawWidths$W_jaw - widthMedian) > 3 * widthMAD, c("specimen", "W_jaw", "planeSine")])

  # Each family-individual pair must be unique, or the kinematics script can't match fish to videos
  if(any(duplicated(paste(jawWidths$Family, jawWidths$Animal.ID)))) {
    stop("The same family and individual appear more than once: ",
         paste(jawWidths$specimen[duplicated(paste(jawWidths$Family, jawWidths$Animal.ID))], collapse = ", "))
  }


##########
#  Save  #
##########

  dir.create(pathOutput, showWarnings = FALSE, recursive = TRUE)
  write.csv(jawWidths, file.path(pathOutput, "jaw_widths_F2.csv"), row.names = FALSE)

  print(paste0("Saved ", nSpec, " jaw widths to: ", file.path(pathOutput, "jaw_widths_F2.csv")))
