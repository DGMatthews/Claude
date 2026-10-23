
#############################################################################################
# Scan measurements from 3D Slicer landmarks                                                 #
#############################################################################################
  # Reads the hard-tissue landmarks (one .mrk.json per fish, 27 points in the order of List_of_landmarks.docx)
    # and the soft-tissue cross-sectional areas (soft_tissue_areas.csv from the Slicer area script),
    # calculates every scan measurement, and saves one row per fish for ACxTRC_Kinematic_Measurements_Sept2026.R
  # All lengths in mm, angles in radians, areas in mm^2 (the Slicer coordinates are already in mm)

  # Skull frame (same baseline the video uses):
    # Left-right axis = average direction of the 4 left-to-right landmark pairs (tooth row, jaw joints, nasals, caudal ceratohyals)
    # Anteroposterior (AP) axis = nasal tip (L) to caudal base of parasphenoid (the scan version of the Eye-Nasal line), with any left-right part removed. Positive = posterior
    # Dorsoventral (DV) axis = perpendicular to both. Positive = ventral
    # Side view = the skull's midline plane (AP and DV axes), which is what the camera sees

  # Four-bar links and resting angles are taken in the side view, so each linkage closes exactly at rest
    # The true 3D lengths are saved as well (*_3D), for genetic mapping
  # Four-bar angles are counterclockwise-positive with the fish facing left (anterior to the left, dorsal up)
    # With this convention, jaw opening and opercular rotation both increase theta2 (oralOpenSign = opercularOpenSign = +1 in the kinematics script)

  # Output: scan_measurements.csv, one row per fish, with columns:
    # specimen, Family, Animal.ID: used to match the kinematic videos
    # Every measurement in List_of_landmarks.docx, plus the four-bar resting angles
    # QC columns (end in _check or _diff_deg), and missingLandmarks / hardFile so any odd value can be traced back

#This script needs the jsonlite package to read the Slicer files
#install.packages("jsonlite")

require(jsonlite)

rm(list=ls())


##############
#  Settings  #
##############

  # Folder with the hard-tissue landmark files (searched including subfolders)
  pathHard <- 'D:/Hybrid feeding project/CT Landmarks/Hard'

  # Summary CSV written by the Slicer area script (slicer_add_areas.py)
  pathSoftSummary <- 'D:/Hybrid feeding project/CT Landmarks/Soft_withArea/soft_tissue_areas.csv'   ### CHANGE if you overwrote the originals instead

  # Where to save the results
  pathOutput <- 'D:/Hybrid feeding project/CT Landmarks'   ### CHANGE
  outFile    <- "scan_measurements.csv"

  # The filenames hold the family as a number (e.g. "1.2" in F2_1.2_11.19_HardLM.mrk.json)
    # The kinematics script matches on the Family column of the feeding data sheet, so the family has to be written the same way
  familyPrefix <- "F"   ### CHECK  "F" turns "1.2" into "F1.2". Set to "" if the data sheet writes families as "1.2"

  # Scans where the fish isn't at rest, so its resting positions can't be measured
    # Mouth open: the resting four-bar angles and the input link (joint to the ventral maxilla tip, which slides along the mandible as the mouth opens) are set to NA
      # The kinematics script then leaves that fish's four-bar outputs as NA instead of starting from a wrong rest position
    # Hyoid not at rest (depressed or abducted): the hyoid zero points and W_head are set to NA, so its hyoid depression and buccal volume are NA
    # Bone lengths are kept, since they don't depend on posture
  mouthOpenScans <- c("F2_1.3_8-11.08")   ### CHANGE  add any other scans with the mouth open
  hyoidNotAtRestScans <- c()              ### CHANGE  add scans where the hyoid is visibly depressed or the ceratohyals spread

  mouthOpenCols <- c("Oral_theta2_rest", "Opercular_theta2_rest", "Opercular_theta4_rest",
                     "Oral_theta4_rest_diff_deg", "Opercular_theta4_rest_diff_deg",
                     "Mandible_input_length", "Mandible_input_length_3D")
  hyoidRestCols <- c("D_head", "D_head_backup", "D_head_backup_3D", "Hyoid_rest_AP",
                     "Ceratohyal_AP_rest", "Ceratohyal_DV_rest", "Ceratohyal_half_width_check", "W_head")

  # Landmark numbers (position in the file), from List_of_landmarks.docx
  lmNames <- c(mandTip      = 1,    # Mandible tip
               pmxTip       = 2,    # Anterior tip of premaxilla – L
               pmxAscArm    = 3,    # Tip of ascending arm of premaxilla – L
               toothL       = 4,    # Widest point of tooth row – L
               toothR       = 5,    # Widest point of tooth row – R
               jointDorsL   = 6,    # Dorsal point mandibular joint – L
               jointVentL   = 7,    # Ventral point mandibular joint – L
               jointDorsR   = 8,    # Dorsal point mandibular joint – R
               jointVentR   = 9,    # Ventral point mandibular joint – R
               raTip        = 10,   # Tip of retroarticular process – L
               maxVent      = 11,   # Ventral maxilla tip – L
               maxDors      = 12,   # Dorsal maxilla tip – L
               palatine     = 13,   # Anterior tip of palatine – L
               nasalL       = 14,   # Anterior tip of nasal – L
               nasalR       = 15,   # Anterior tip of nasal – R
               neuroTip     = 16,   # Anterior tip of neurocranium – L
               iopTip       = 17,   # Anterior tip of IOP – L
               iopSopSym    = 18,   # Symphysis of interoperculum and suboperculum – L
               hyoOpJoint   = 19,   # Hyomandibular-opercular joint – L
               crestBase    = 20,   # Base of supraoccipital crest – L
               crestTip     = 21,   # Dorsal tip of supraoccipital crest
               urohyalRost  = 22,   # Rostral tip urohyal – L
               urohyalCaud  = 23,   # Caudal tip urohyal – L
               ceratoRost   = 24,   # Rostral tip ceratohyal – L (= tracked Hyoid point)
               ceratoCaudL  = 25,   # Caudal tip of ceratohyal – L
               ceratoCaudR  = 26,   # Caudal tip of ceratohyal – R
               parasph      = 27)   # Caudal base of parasphenoid (= eye-level landmark)
  nLandmarks <- 27

  # Soft-tissue files are named <specimen>_<code>.mrk.json. Code = name in the output
  softCodes <- c(Epaxial   = "Epaxial_CSA_crest",        ### CHECK  Epaxial   = at the highest point of the supraoccipital crest
                 Epaxial_2 = "Epaxial_CSA_flat",         ### CHECK  Epaxial_2 = posterior-most section where the neurocranium flattens out
                 Hypaxial  = "Hypaxial_CSA",
                 SH        = "Sternohyoideus_CSA",
                 IOPL      = "IOP_ligament_CSA",
                 PH        = "Protractor_hyoideus_CSA")

  # Checks
  maxPairAngle   <- 30     ### CHANGE  stop if any left-right pair points more than this many degrees away from the average left-right axis (sides swapped or a misplaced point)
  maxOffPlane    <- 0.1    ### CHANGE  mm. Cross-sections that stray further than this from one plane are listed
  maxAreaMismatch <- 0.02  ### CHANGE  Slicer's own area and the planar area may differ by this fraction before the file is listed


#####################
#  Define functions #
#####################

    # Read all landmarks from a Slicer .mrk.json file, as a 27 x 3 matrix (rows = landmark number)
      # Landmark ids and labels aren't used, only the order of the points
      # Points not marked "defined" (skipped or missing in Slicer) are returned as NA
    readLandmarks <- function(filePath) {

      controlPoints <- fromJSON(filePath)$markups$controlPoints[[1]]

      if(nrow(controlPoints) != nLandmarks) {
        stop(nrow(controlPoints), " points instead of ", nLandmarks, " in: ", filePath)
      }

      pts <- do.call(rbind, controlPoints$position)
      if(!is.null(controlPoints$positionStatus)) {
        pts[controlPoints$positionStatus != "defined", ] <- NA
      }
      pts
    }


    # Specimen name from a filename, e.g. "F2_1.2_11.19_HardLM.mrk.json" -> "F2_1.2_11.19"
      # A leading "ACxTRC_" is removed, as in the other scripts
    specimenFromFile <- function(fileName, suffix) {
      specimenNow <- sub(paste0("_", suffix, "(\\.mrk)?\\.json$"), "", basename(fileName))
      sub("^ACxTRC_", "", specimenNow)
    }


    # Vector helpers
    vecLength    <- function(v) sqrt(sum(v^2))
    unitVec      <- function(v) v / vecLength(v)
    crossProduct <- function(u, v) {
      c(u[2]*v[3] - u[3]*v[2],
        u[3]*v[1] - u[1]*v[3],
        u[1]*v[2] - u[2]*v[1])
    }


    # Four-bar output angle theta4 for a given input angle theta2 (radians)
      # Same function (Freudenstein equation, + solution) as the kinematics script. Used here only to check that the measured linkage matches it
    fourBarTheta4 <- function(theta2, inputLength, couplerLength, outputLength, fixedLength) {
      K1 <- fixedLength / inputLength
      K2 <- fixedLength / outputLength
      K3 <- (inputLength^2 - couplerLength^2 + outputLength^2 + fixedLength^2) / (2 * inputLength * outputLength)

      A <- K1 - cos(theta2)
      B <- -sin(theta2)
      C <- K2 * cos(theta2) - K3

      2 * atan((B + sqrt(A^2 + B^2 - C^2)) / (A + C))
    }


    # Wrap an angle (radians) into -pi to pi
    wrapAngle <- function(angle) {
      ((angle + pi) %% (2*pi)) - pi
    }


###################
#  Find the files #
###################

  hardFiles <- list.files(pathHard, pattern = "_HardLM(\\.mrk)?\\.json$", recursive = TRUE, full.names = TRUE)

  if(length(hardFiles) == 0) {
    stop("No *_HardLM.mrk.json files in: ", pathHard)
  }

  # Anything else in the folder is listed, in case a file is misnamed
  otherHard <- setdiff(list.files(pathHard, pattern = "\\.json$", recursive = TRUE, full.names = TRUE), hardFiles)
  if(length(otherHard) > 0) {
    print("Other .json files in the hard-tissue folder (left out):")
    print(otherHard)
  }

  hardSpecimens <- specimenFromFile(hardFiles, "HardLM")

  if(any(duplicated(hardSpecimens))) {
    stop("More than one hard-tissue file for the same specimen: ", paste(unique(hardSpecimens[duplicated(hardSpecimens)]), collapse = ", "))
  }

  # Only F2 specimens. Anything else is listed and left out
  notF2 <- hardSpecimens[!startsWith(hardSpecimens, "F2_")]
  if(length(notF2) > 0) {
    print("Not F2, left out:")
    print(notF2)
  }
  hardFiles     <- hardFiles[startsWith(hardSpecimens, "F2_")]
  hardSpecimens <- hardSpecimens[startsWith(hardSpecimens, "F2_")]

  nSpec <- length(hardSpecimens)
  print(paste0("Measuring ", nSpec, " F2 specimens"))


###########################
#  Hard-tissue measures   #
###########################

  # Output columns, in order
  hardCols <- c(
    # Widths
    "W_jaw", "W_jaw_joint", "W_nasal", "W_head",
    # Mandible triangle (side-view lengths, by the List_of_landmarks.docx formulas)
    "Mandible_length", "Nasal_QMJ_dist",
    # Hyoid and buccal volume
    "D_head", "D_head_backup", "Nasal_eye_AP_length", "Hyoid_rest_AP",
    "Ceratohyal_length", "Ceratohyal_AP_rest", "Ceratohyal_DV_rest",
    # Other bones
    "Neurocranium_length", "Supraoccipital_crest_depth", "Urohyal_length", "Length_premaxilla_asc_arm",
    # Oral four-bar (side view)
    "Mandible_input_length", "Maxilla_length", "Maxilla_palatine_length", "Oral_fixed_length", "Oral_theta2_rest",
    # Opercular four-bar (side view)
    "Operculum_length", "IOP_length", "IOP_ligament_length", "RA_process_length", "Opercular_fixed_length",
    "Opercular_theta2_rest", "Opercular_theta4_rest",
    # True 3D sizes, for genetic mapping
    "Mandible_length_3D", "Nasal_QMJ_dist_3D", "D_head_backup_3D",
    "Mandible_input_length_3D", "Maxilla_length_3D", "Maxilla_palatine_length_3D", "Oral_fixed_length_3D",
    "Operculum_length_3D", "IOP_length_3D", "IOP_ligament_length_3D", "RA_process_length_3D", "Opercular_fixed_length_3D",
    # Checks
    "pairAngle_max_check", "Ceratohyal_half_width_check", "Oral_theta4_rest_diff_deg", "Opercular_theta4_rest_diff_deg")

  scanMeas <- data.frame(specimen  = hardSpecimens,
                         Family    = rep(NA_character_, nSpec),
                         Animal.ID = rep(NA_character_, nSpec))
  for(colNow in hardCols) {
    scanMeas[[colNow]] <- rep(NA_real_, nSpec)
  }
  scanMeas$missingLandmarks <- rep("", nSpec)
  scanMeas$hardFile         <- hardFiles


  for(i in 1:nSpec) {

    ##### Family and individual from the specimen name, e.g. "F2_1.2_11.19" -> family "1.2", individual "11.19"
      nameParts <- strsplit(hardSpecimens[i], "_")[[1]]

      if(length(nameParts) != 3) {
        stop("Can't read family and individual from '", hardSpecimens[i], "'. Expected F2_family_individual")
      }

      scanMeas$Family[i]    <- paste0(familyPrefix, nameParts[2])
      scanMeas$Animal.ID[i] <- nameParts[3]


    ##### Landmarks (mm)

      pts <- readLandmarks(hardFiles[i])
      lmk <- function(name) pts[lmNames[[name]], ]

      missingNow <- which(is.na(pts[, 1]))
      scanMeas$missingLandmarks[i] <- paste(missingNow, collapse = ";")

      # Jaw joint centres = average of the dorsal and ventral joint points
      jointL <- (lmk("jointDorsL") + lmk("jointVentL")) / 2
      jointR <- (lmk("jointDorsR") + lmk("jointVentR")) / 2


    ##### Skull frame
      # The frame needs the left-right pairs, nasal tip, parasphenoid and supraoccipital crest
      if(any(is.na(c(lmk("toothL"), lmk("toothR"), jointL, jointR, lmk("nasalL"), lmk("nasalR"), lmk("ceratoCaudL"), lmk("ceratoCaudR"),
                     lmk("parasph"), lmk("crestBase"), lmk("crestTip"))))) {
        print(paste0(hardSpecimens[i], ": a skull-frame landmark is missing, so the fish is left out"))
        next
      }

      # Left-right axis: average of the 4 left-to-right directions
      pairDirs <- rbind(unitVec(lmk("toothR")      - lmk("toothL")),
                        unitVec(jointR            - jointL),
                        unitVec(lmk("nasalR")      - lmk("nasalL")),
                        unitVec(lmk("ceratoCaudR") - lmk("ceratoCaudL")))
      lrAxis <- unitVec(colMeans(pairDirs))

      # Every pair should point roughly the same way. If one doesn't, its sides are probably swapped
      pairAngles <- acos(pmin(1, pairDirs %*% lrAxis)) * 180 / pi
      scanMeas$pairAngle_max_check[i] <- max(pairAngles)
      if(max(pairAngles) > maxPairAngle) {
        stop(hardSpecimens[i], ": a left-right pair is ", round(max(pairAngles)), " degrees off the left-right axis ",
             "(tooth row, joint, nasal, ceratohyal = ", paste(round(pairAngles), collapse = ", "), "). Check for swapped sides")
      }

      # AP axis: nasal tip to parasphenoid, without its left-right part. Positive = posterior
      apAxis <- lmk("parasph") - lmk("nasalL")
      apAxis <- unitVec(apAxis - sum(apAxis * lrAxis) * lrAxis)

      # DV axis: perpendicular to both. Positive = ventral
      dvAxis <- crossProduct(lrAxis, apAxis)

      # The supraoccipital crest tip must be dorsal to its base. If not, the axes (or the coordinate system) are flipped
      if(sum((lmk("crestTip") - lmk("crestBase")) * dvAxis) >= 0) {
        stop(hardSpecimens[i], ": the supraoccipital crest tip isn't dorsal to its base in the skull frame. Check the landmarks")
      }


    ##### Distance helpers for this fish
      dist3D   <- function(a, b) vecLength(b - a)
      # Side view: distance with the left-right part removed
      distSide <- function(a, b) {
        v <- b - a
        vecLength(v - sum(v * lrAxis) * lrAxis)
      }
      apComp <- function(a, b) sum((b - a) * apAxis)   # positive if b is posterior to a
      dvComp <- function(a, b) sum((b - a) * dvAxis)   # positive if b is ventral to a

      # Side-view 2D coordinates: x = posterior (+AP), y = dorsal (-DV), so the fish faces left
        # Angles in this plane are counterclockwise-positive
      side2D <- function(p) c(sum(p * apAxis), -sum(p * dvAxis))
      # Signed angle (radians) from direction u to direction v, both 2D
      signedAngle <- function(u, v) atan2(u[1]*v[2] - u[2]*v[1], sum(u * v))
      # Side-view direction from a to b
      sideDir <- function(a, b) side2D(b) - side2D(a)


    ##### Widths
      scanMeas$W_jaw[i]       <- dist3D(lmk("toothL"), lmk("toothR"))
      scanMeas$W_jaw_joint[i] <- dist3D(jointL, jointR)
      scanMeas$W_nasal[i]     <- dist3D(lmk("nasalL"), lmk("nasalR"))
      scanMeas$W_head[i]      <- dist3D(lmk("ceratoCaudL"), lmk("ceratoCaudR"))


    ##### Mandible triangle
      # Mandible_length = sqrt(hypotenuse^2 - (W_jaw_joint/2)^2): mandible tip (midline) to left joint, without the joint's distance from the midline
      scanMeas$Mandible_length_3D[i] <- dist3D(lmk("mandTip"), jointL)
      scanMeas$Mandible_length[i]    <- sqrt(scanMeas$Mandible_length_3D[i]^2 - (scanMeas$W_jaw_joint[i] / 2)^2)

      # Nasal_QMJ_dist = sqrt(hypotenuse^2 - ((W_jaw_joint - W_nasal)/2)^2): left nasal tip to left joint, without the side-to-side gap between them
      scanMeas$Nasal_QMJ_dist_3D[i] <- dist3D(lmk("nasalL"), jointL)
      scanMeas$Nasal_QMJ_dist[i]    <- sqrt(scanMeas$Nasal_QMJ_dist_3D[i]^2 - ((scanMeas$W_jaw_joint[i] - scanMeas$W_nasal[i]) / 2)^2)


    ##### Hyoid and buccal volume

      # D_head: parasphenoid base to the line through the two caudal ceratohyal tips (perpendicular distance)
      lineDir  <- unitVec(lmk("ceratoCaudR") - lmk("ceratoCaudL"))
      toPoint  <- lmk("parasph") - lmk("ceratoCaudL")
      scanMeas$D_head[i] <- vecLength(toPoint - sum(toPoint * lineDir) * lineDir)

      # D_head_backup: nasal tip to the tracked Hyoid point (rostral ceratohyal tip) at rest, side view like the video
      scanMeas$D_head_backup[i]    <- distSide(lmk("nasalL"), lmk("ceratoRost"))
      scanMeas$D_head_backup_3D[i] <- dist3D(lmk("nasalL"), lmk("ceratoRost"))

      # AP distances along the Eye-Nasal baseline
      scanMeas$Nasal_eye_AP_length[i] <- apComp(lmk("nasalL"), lmk("parasph"))
      scanMeas$Hyoid_rest_AP[i]       <- apComp(lmk("parasph"), lmk("ceratoRost"))   # positive if the rostral ceratohyal tip is posterior to the parasphenoid

      # Ceratohyal bar at rest
      scanMeas$Ceratohyal_length[i]  <- dist3D(lmk("ceratoRost"), lmk("ceratoCaudL"))
      scanMeas$Ceratohyal_AP_rest[i] <- apComp(lmk("ceratoRost"), lmk("ceratoCaudL"))   # positive if the caudal tip is posterior
      scanMeas$Ceratohyal_DV_rest[i] <- dvComp(lmk("ceratoCaudL"), lmk("ceratoRost"))   # positive if the rostral tip is ventral

      # Check: the bar model's resting half-width should be close to W_head/2 (ratio near 1)
        # NaN if the scan values can't form the bar at rest (the kinematics script stops on these)
      scanMeas$Ceratohyal_half_width_check[i] <- sqrt(scanMeas$Ceratohyal_length[i]^2 - scanMeas$Ceratohyal_AP_rest[i]^2 - scanMeas$Ceratohyal_DV_rest[i]^2) /
                                                 (scanMeas$W_head[i] / 2)


    ##### Other bones (3D)
      scanMeas$Neurocranium_length[i]        <- dist3D(lmk("neuroTip"), lmk("crestBase"))
      scanMeas$Supraoccipital_crest_depth[i] <- dist3D(lmk("crestBase"), lmk("crestTip"))
      scanMeas$Urohyal_length[i]             <- dist3D(lmk("urohyalRost"), lmk("urohyalCaud"))
      scanMeas$Length_premaxilla_asc_arm[i]  <- dist3D(lmk("pmxTip"), lmk("pmxAscArm"))


    ##### Oral four-bar
      # Input = mandible (joint to ventral maxilla tip), coupler = maxilla, output = dorsal maxilla tip to palatine, fixed = joint to palatine
      scanMeas$Mandible_input_length[i]   <- distSide(jointL, lmk("maxVent"))
      scanMeas$Maxilla_length[i]          <- distSide(lmk("maxVent"), lmk("maxDors"))
      scanMeas$Maxilla_palatine_length[i] <- distSide(lmk("maxDors"), lmk("palatine"))
      scanMeas$Oral_fixed_length[i]       <- distSide(jointL, lmk("palatine"))

      scanMeas$Mandible_input_length_3D[i]   <- dist3D(jointL, lmk("maxVent"))
      scanMeas$Maxilla_length_3D[i]          <- dist3D(lmk("maxVent"), lmk("maxDors"))
      scanMeas$Maxilla_palatine_length_3D[i] <- dist3D(lmk("maxDors"), lmk("palatine"))
      scanMeas$Oral_fixed_length_3D[i]       <- dist3D(jointL, lmk("palatine"))

      # theta2 = angle at the joint from the fixed link (joint -> palatine) to the input link (joint -> ventral maxilla tip)
      oralFixedDir <- sideDir(jointL, lmk("palatine"))
      scanMeas$Oral_theta2_rest[i] <- signedAngle(oralFixedDir, sideDir(jointL, lmk("maxVent")))

      # Check: measured theta4 (fixed link direction to output link, at the palatine) vs. the Freudenstein + solution
        # Should be ~0. A large value means this linkage sits in the other assembly, and the kinematics script's four-bar would be wrong for it
      oralTheta4Meas <- signedAngle(oralFixedDir, sideDir(lmk("palatine"), lmk("maxDors")))
      oralTheta4Pred <- fourBarTheta4(scanMeas$Oral_theta2_rest[i], scanMeas$Mandible_input_length[i], scanMeas$Maxilla_length[i],
                                      scanMeas$Maxilla_palatine_length[i], scanMeas$Oral_fixed_length[i])
      scanMeas$Oral_theta4_rest_diff_deg[i] <- wrapAngle(oralTheta4Pred - oralTheta4Meas) * 180 / pi


    ##### Opercular four-bar
      # Input = operculum (hyomandibular-opercular joint to IOP-SOP symphysis), coupler = IOP + IOP ligament,
        # output = retroarticular process (jaw joint to RA tip), fixed = hyomandibular-opercular joint to jaw joint
      scanMeas$Operculum_length[i]       <- distSide(lmk("hyoOpJoint"), lmk("iopSopSym"))
      scanMeas$IOP_length[i]             <- distSide(lmk("iopTip"), lmk("iopSopSym"))
      scanMeas$IOP_ligament_length[i]    <- distSide(lmk("raTip"), lmk("iopTip"))
      scanMeas$RA_process_length[i]      <- distSide(jointL, lmk("raTip"))
      scanMeas$Opercular_fixed_length[i] <- distSide(jointL, lmk("hyoOpJoint"))

      scanMeas$Operculum_length_3D[i]       <- dist3D(lmk("hyoOpJoint"), lmk("iopSopSym"))
      scanMeas$IOP_length_3D[i]             <- dist3D(lmk("iopTip"), lmk("iopSopSym"))
      scanMeas$IOP_ligament_length_3D[i]    <- dist3D(lmk("raTip"), lmk("iopTip"))
      scanMeas$RA_process_length_3D[i]      <- dist3D(jointL, lmk("raTip"))
      scanMeas$Opercular_fixed_length_3D[i] <- dist3D(jointL, lmk("hyoOpJoint"))

      # theta2 = angle at the hyomandibular-opercular joint from the fixed link (-> jaw joint) to the input link (-> symphysis)
      # theta4 = angle at the jaw joint from the same fixed-link direction to the output link (-> RA tip)
      opFixedDir <- sideDir(lmk("hyoOpJoint"), jointL)
      scanMeas$Opercular_theta2_rest[i] <- signedAngle(opFixedDir, sideDir(lmk("hyoOpJoint"), lmk("iopSopSym")))
      scanMeas$Opercular_theta4_rest[i] <- signedAngle(opFixedDir, sideDir(jointL, lmk("raTip")))

      # Check, as for the oral four-bar. Here it isn't exactly 0 even when correct, because the coupler is IOP + ligament
        # (bent at the IOP tip) rather than the straight symphysis-to-RA distance. The kinematics script reports the same difference
      opTheta4Pred <- fourBarTheta4(scanMeas$Opercular_theta2_rest[i], scanMeas$Operculum_length[i],
                                    scanMeas$IOP_length[i] + scanMeas$IOP_ligament_length[i],
                                    scanMeas$RA_process_length[i], scanMeas$Opercular_fixed_length[i])
      scanMeas$Opercular_theta4_rest_diff_deg[i] <- wrapAngle(opTheta4Pred - scanMeas$Opercular_theta4_rest[i]) * 180 / pi
  }


  ##### Scans that aren't at rest
    # Names that don't match any scan are listed, in case of a typo
  notFound <- setdiff(c(mouthOpenScans, hyoidNotAtRestScans), scanMeas$specimen)
  if(length(notFound) > 0) {
    print("Scans listed as not at rest that aren't in the hard-tissue folder (check the names):")
    print(notFound)
  }

  scanMeas[scanMeas$specimen %in% mouthOpenScans, mouthOpenCols]      <- NA
  scanMeas[scanMeas$specimen %in% hyoidNotAtRestScans, hyoidRestCols] <- NA

  scanMeas$notAtRest <- ifelse(scanMeas$specimen %in% mouthOpenScans,
                               ifelse(scanMeas$specimen %in% hyoidNotAtRestScans, "mouth open; hyoid", "mouth open"),
                               ifelse(scanMeas$specimen %in% hyoidNotAtRestScans, "hyoid", ""))


###########################
#  Soft-tissue areas      #
###########################

  softSummary <- read.csv(pathSoftSummary, stringsAsFactors = FALSE)

  # Specimen and muscle code from the filename, e.g. "F2_1.2_11.01\F2_1.2_11.01_Epaxial_2.mrk.json" -> "F2_1.2_11.01", "Epaxial_2"
  softBase <- sub("^ACxTRC_", "", basename(gsub("\\\\", "/", softSummary$file)))
  softPattern <- "^(F2_[^_]+_[^_]+)_(.+?)(\\.mrk)?\\.json$"
  isCSA <- grepl(softPattern, softBase) & sub(softPattern, "\\2", softBase) %in% names(softCodes)

  if(any(!isCSA)) {
    print("Soft-tissue files that aren't one of the CSA codes (left out):")
    print(softSummary$file[!isCSA])
  }

  softSummary <- softSummary[isCSA, ]
  softSummary$specimen <- sub(softPattern, "\\1", softBase[isCSA])
  softSummary$code     <- sub(softPattern, "\\2", softBase[isCSA])

  if(any(duplicated(paste(softSummary$specimen, softSummary$code)))) {
    stop("More than one area file for the same specimen and muscle: ",
         paste(unique(paste(softSummary$specimen, softSummary$code)[duplicated(paste(softSummary$specimen, softSummary$code))]), collapse = ", "))
  }

  softSummary$area_mm2       <- suppressWarnings(as.numeric(softSummary$area_mm2))
  softSummary$slicerArea     <- suppressWarnings(as.numeric(softSummary$slicerArea))
  softSummary$maxOffPlane_mm <- suppressWarnings(as.numeric(softSummary$maxOffPlane_mm))

  # Files with no area (errors in the Slicer script)
  noArea <- is.na(softSummary$area_mm2)
  if(any(noArea)) {
    print("CSA files with no area (left as NA):")
    print(softSummary[noArea, c("file", "note")])
  }

  # Check: Slicer's own area vs. the planar area. Slicer usually reports cm2
  unitFactor <- ifelse(softSummary$slicerUnits == "cm2", 100, ifelse(softSummary$slicerUnits == "mm2", 1, NA))
  areaMismatch <- abs(softSummary$slicerArea * unitFactor / softSummary$area_mm2 - 1)
  if(any(areaMismatch > maxAreaMismatch, na.rm = TRUE)) {
    print("CSA files where Slicer's area and the planar area disagree (check the curve for a crossing):")
    print(softSummary[which(areaMismatch > maxAreaMismatch), c("file", "slicerArea", "slicerUnits", "area_mm2")])
  }

  # Check: cross-sections that aren't flat
  if(any(softSummary$maxOffPlane_mm > maxOffPlane, na.rm = TRUE)) {
    print(paste0("CSA files more than ", maxOffPlane, " mm off one plane:"))
    print(softSummary[which(softSummary$maxOffPlane_mm > maxOffPlane), c("file", "maxOffPlane_mm")])
  }

  # One column per muscle, one row per specimen
  softSpecimens <- sort(unique(softSummary$specimen))
  softMeas <- data.frame(specimen = softSpecimens)
  for(codeNow in names(softCodes)) {
    rowsNow <- softSummary[softSummary$code == codeNow, ]
    softMeas[[softCodes[[codeNow]]]] <- rowsNow$area_mm2[match(softSpecimens, rowsNow$specimen)]
  }


###########################
#  Combine and check      #
###########################

  # Fish with only one of the two scan types are kept, with NA for the other
  onlyHard <- setdiff(scanMeas$specimen, softMeas$specimen)
  onlySoft <- setdiff(softMeas$specimen, scanMeas$specimen)
  print("Hard-tissue landmarks but no soft-tissue areas:")
  print(onlyHard)
  print("Soft-tissue areas but no hard-tissue landmarks (Family and Animal.ID filled in from the name):")
  print(onlySoft)

  scanOut <- merge(scanMeas, softMeas, by = "specimen", all = TRUE)

  # Family and individual for the soft-only fish
  noID <- is.na(scanOut$Family)
  scanOut$Family[noID]    <- paste0(familyPrefix, sapply(strsplit(scanOut$specimen[noID], "_"), `[`, 2))
  scanOut$Animal.ID[noID] <- sapply(strsplit(scanOut$specimen[noID], "_"), `[`, 3)

  # Each family-individual pair must be unique, or the kinematics script can't match fish to videos
  if(any(duplicated(paste(scanOut$Family, scanOut$Animal.ID)))) {
    stop("The same family and individual appear more than once: ",
         paste(scanOut$specimen[duplicated(paste(scanOut$Family, scanOut$Animal.ID))], collapse = ", "))
  }


  # Missing landmarks
  if(any(scanOut$missingLandmarks != "", na.rm = TRUE)) {
    print("Fish with missing landmarks (measurements that need them are NA):")
    print(scanOut[which(scanOut$missingLandmarks != ""), c("specimen", "missingLandmarks")])
  }

  # Range of every measurement
  measCols <- c(hardCols, unname(softCodes))
  print(summary(scanOut[, measCols]))

  # Values far from the median (more than 4 MADs) = worth opening in Slicer
  print("Values more than 4 MADs from the median:")
  for(colNow in setdiff(measCols, c("pairAngle_max_check", "Oral_theta4_rest_diff_deg", "Opercular_theta4_rest_diff_deg"))) {
    valNow <- scanOut[[colNow]]
    if(sum(!is.na(valNow)) < 3 || mad(valNow, na.rm = TRUE) == 0) {
      next
    }
    farNow <- which(abs(valNow - median(valNow, na.rm = TRUE)) > 4 * mad(valNow, na.rm = TRUE))
    if(length(farNow) > 0) {
      print(paste0(colNow, ": ", paste(scanOut$specimen[farNow], collapse = ", ")))
    }
  }

  # Ceratohyal bar check: ratios far from 1 mean the ceratohyal tips aren't where the bar model assumes
  print("Ceratohyal half-width / (W_head/2):")
  print(summary(scanOut$Ceratohyal_half_width_check))

  # Four-bar checks (degrees)
  print("Oral four-bar: predicted - measured theta4 at rest (should be ~0):")
  print(summary(scanOut$Oral_theta4_rest_diff_deg))
  print(scanOut$specimen[which(abs(scanOut$Oral_theta4_rest_diff_deg) > 1)])
  print("Opercular four-bar: predicted - measured theta4 at rest (small, from the bend at the IOP tip):")
  print(summary(scanOut$Opercular_theta4_rest_diff_deg))

  # Linkages that can't close at rest
  print("Fish whose four-bar can't close at rest (NaN):")
  print(scanOut$specimen[which(is.nan(scanOut$Oral_theta4_rest_diff_deg) | is.nan(scanOut$Opercular_theta4_rest_diff_deg))])


##########
#  Save  #
##########

  dir.create(pathOutput, showWarnings = FALSE, recursive = TRUE)
  write.csv(scanOut, file.path(pathOutput, outFile), row.names = FALSE)

  print(paste0("Saved ", nrow(scanOut), " fish to: ", file.path(pathOutput, outFile)))
