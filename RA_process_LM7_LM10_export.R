#############################################################################################
# One-off: distance between landmark 10 (tip of retroarticular process – L)                 #
# and landmark 7 (ventral point mandibular joint – L), for comparison with a previous analysis #
#############################################################################################
  # Reads the same hard-tissue landmark files as 1_Scan_measurements_from_landmarks.R (one .mrk.json per fish)
  # and writes one row per specimen with the 3D distance between the two landmarks (mm, Slicer coordinates are already in mm)
  # For reference it also writes the 3D distance from landmark 10 to the jaw joint centre (midpoint of landmarks 6 and 7),
    # which is what RA_process_length_3D in scan_measurements.csv uses

#This script needs the jsonlite package to read the Slicer files
#install.packages("jsonlite")

require(jsonlite)

rm(list=ls())


##############
#  Settings  #
##############

  # Folder with the hard-tissue landmark files (searched including subfolders), as in script 1
  pathHard <- 'D:/Hybrid feeding project/CT Landmarks/Hard'   ### CHANGE

  # Where to save the results
  outFile <- 'C:/Users/Dave/Documents/RealDocuments/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data/RA_process_LM7_LM10.csv'   ### CHANGE

  # Only F2 specimens (as script 1)? FALSE keeps every file, e.g. if the previous analysis included other fish
  onlyF2 <- FALSE   ### CHANGE

  # Landmark numbers (position in the file), from List_of_landmarks.docx
  lmJointDorsL <- 6    # Dorsal point mandibular joint – L
  lmJointVentL <- 7    # Ventral point mandibular joint – L
  lmRaTip      <- 10   # Tip of retroarticular process – L


#####################
#  Define functions #
#####################

  # Read all landmarks from a Slicer .mrk.json file, as an n x 3 matrix (rows = landmark number)
    # Points not marked "defined" (skipped or missing in Slicer) are returned as NA
  readLandmarks <- function(filePath) {
    controlPoints <- fromJSON(filePath)$markups$controlPoints[[1]]
    pts <- do.call(rbind, controlPoints$position)
    if(!is.null(controlPoints$positionStatus)) {
      pts[controlPoints$positionStatus != "defined", ] <- NA
    }
    pts
  }

  # Specimen name from a filename, e.g. "F2_1.2_11.19_HardLM.mrk.json" -> "F2_1.2_11.19" (a leading "ACxTRC_" is removed, as in script 1)
  specimenFromFile <- function(fileName) {
    sub("^ACxTRC_", "", sub("_HardLM(\\.mrk)?\\.json$", "", basename(fileName)))
  }

  dist3D <- function(a, b) sqrt(sum((b - a)^2))


###################
#  Find the files #
###################

  hardFiles <- list.files(pathHard, pattern = "_HardLM(\\.mrk)?\\.json$", recursive = TRUE, full.names = TRUE)
  if(length(hardFiles) == 0) {
    stop("No *_HardLM.mrk.json files in: ", pathHard)
  }
  specimens <- specimenFromFile(hardFiles)

  if(onlyF2) {
    keep <- startsWith(specimens, "F2_")
    hardFiles <- hardFiles[keep]
    specimens <- specimens[keep]
  }
  print(paste0("Measuring ", length(hardFiles), " specimens"))


##################
#  Measurements  #
##################

  out <- do.call(rbind, lapply(seq_along(hardFiles), function(i) {
    # Family and individual from names like "F2_1.2_11.19" (NA if the name has another form)
    nameParts <- strsplit(specimens[i], "_")[[1]]
    familyNow <- if(length(nameParts) == 3) nameParts[2] else NA_character_
    animalNow <- if(length(nameParts) == 3) nameParts[3] else NA_character_

    pts <- tryCatch(readLandmarks(hardFiles[i]), error = function(e) NULL)
    noteNow <- ""
    if(is.null(pts) || nrow(pts) < lmRaTip) {
      noteNow <- "could not read the file, or fewer than 10 landmarks"
      lm7 <- lm6 <- lm10 <- rep(NA_real_, 3)
    } else {
      lm6  <- pts[lmJointDorsL, ]
      lm7  <- pts[lmJointVentL, ]
      lm10 <- pts[lmRaTip, ]
      missingNow <- c(6, 7, 10)[c(any(is.na(lm6)), any(is.na(lm7)), any(is.na(lm10)))]
      if(length(missingNow) > 0) noteNow <- paste0("missing landmark(s): ", paste(missingNow, collapse = ";"))
    }

    data.frame(specimen = specimens[i], Family = familyNow, Animal.ID = animalNow,
               LM7_LM10_dist_3D = dist3D(lm7, lm10),
               jointCentre_LM10_dist_3D = dist3D((lm6 + lm7) / 2, lm10),
               LM7_x = lm7[1], LM7_y = lm7[2], LM7_z = lm7[3],
               LM10_x = lm10[1], LM10_y = lm10[2], LM10_z = lm10[3],
               note = noteNow, hardFile = hardFiles[i], stringsAsFactors = FALSE)
  }))

  print(summary(out$LM7_LM10_dist_3D))
  if(any(out$note != "")) {
    print("Specimens with problems:")
    print(out[out$note != "", c("specimen", "note")], row.names = FALSE)
  }

  write.csv(out, outFile, row.names = FALSE)
  print(paste0("Saved: ", outFile))
