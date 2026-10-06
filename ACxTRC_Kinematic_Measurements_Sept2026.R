
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

#First view the structure of the file
h5ls(file=fileNames[1])


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
                            Neurocranium_tip_disp_eye_t = array(data=NA, dim = c(1,longestVid, nVids), dimnames = list("Neurocranium_tip_disp_eye_t", NULL, fileNames2))
                            )




#Make a list of videos where an analysis was skipped so I can check the landmarking quality
skippedVids <- vector()







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

head(vidData)
head(finalData)

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

finalData$Event.ID[which(duplicated(finalData$Event.ID)==TRUE)] #No duplicates, onto the analysis

length(unique(finalData$Unique.ID))




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
  sum(!is.na(finalData$tmax))
  setdiff(Pressure_vel_Data$video, vidData$PIV.filename[finalData$equivVidRow])
  
  
  
  
  
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
        Hyoid_rest                 = "Hyoid_rest",               ### CHANGE  resting hyoid position proxy
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
    colSums(is.na(finalData[, names(scanColumns)]))
      # Which individuals have no scan match at all?
    unique(finalData$Unique.ID[is.na(finalData$Mandible_length)])
  
  
  
  
  
    
    
    
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
        unique(finalData$Unique.ID[is.nan(finalData$Oral_theta3_rest) | is.nan(finalData$Oral_theta3_30)])
        unique(finalData$Unique.ID[is.nan(finalData$Opercular_theta4_rest) | is.nan(finalData$Opercular_theta4_10)])
        
          # Calculated vs. scanned theta4_rest, in degrees. Should be close to 0
        if("Opercular_theta4_rest_scan" %in% names(finalData)) {
          finalData$Opercular_theta4_rest_diff <- wrapAngle(finalData$Opercular_theta4_rest - finalData$Opercular_theta4_rest_scan)
          summary(finalData$Opercular_theta4_rest_diff * 180 / pi)
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
        
        lambdaRef <- 5e-6
          lambdaRefDerived <- lambdaRef/5   # smoothing for derived series (R_t etc.). Should do less smoothing than the original pass.
        spanRef <- 73
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
          
          
          
          
          
          
                    
                        # Time_hyoid
                        # Time_cranial
                        # U_ram_mean

                        # Hyoid_depression_t_mean

                        
                        # Maxilla_tip_disp_tmax
                        # Maxilla_linear_vel_tmax
                        # Maxilla_ang_tmax
                        # Maxilla_angular_vel_tmax
                        # Hyoid_depression_tmax
                        # Hyoid_vel_tmax
                        # U_ram_tmax
                        # U_ram_acc_tmax
                        # Kinetic_Synchronization

                        # Protrusion_ROM
                        # Mandible_ROM
                        # Hyoid_ROM
                        # Neurocranium_ROM



        # finalTimeSeriesData: x_mouth_t, Hyoid_depression_t













































        #Run through each frame and find the gape (distance from upper to lower jaw), save the max value and time of max
        tempVect <- rep(NA, length(dataNow[,1,1,1]))
        for (k in 1:length(dataNow[,1,1,1])) {
          xUJnow <- dataNow[k,1,1,1]
          yUJnow <- dataNow[k,1,2,1]

          xLJnow <- dataNow[k,2,1,1]
          yLJnow <- dataNow[k,2,2,1]

          gapeNowpix <- sqrt((xUJnow-xLJnow)^2+(yUJnow-yLJnow)^2)

          gapeNow <- gapeNowpix/finalData$cal_pix.mm[i]

          if(is.na(gapeNow)==FALSE) {
            if(gapeNow>=0) {
              tempVect[k] <- gapeNow
            }
          }
        }


        #This runs a polynomial fit on gape, then saves out the max gape from the fitted values

        if (sum(!is.na(tempVect))>10) {

          gapeNow <- na.omit(tempVect)
          framesNow <- seq(from=1, to=longestVid, by=1)[which(!is.na(tempVect))]


          #Create the first polynomial fit, and then get the positions that this line represents
          polyFit <- lm(gapeNow ~ poly(framesNow, 6,raw = TRUE))
          errorFitLine <- predict(polyFit)

          #plot(gapeNow~framesNow)
          #lines(errorFitLine~framesNow)

          #Save the polynomial fit data to gapeDF
          for(b in 1:length(errorFitLine)){
            gapeDF[framesNow[b],1,i] <- errorFitLine[b]
          }

          timeMaxGape <- which(gapeDF[,1,i]==max(gapeDF[,1,i], na.rm = TRUE))

          #Save the max best fit gape
          finalData$maxGape[i] <- max(gapeDF[,1,i], na.rm=TRUE)
          finalData$timeMaxGape[i] <- timeMaxGape


          #Now find the gape velocity
          #Only do this if the maximum gape falls outside the first 10 frames in the video
          if (timeMaxGape>9){
            #Make an expression to represent the best fit polynomial from lm()
            gapeBestFit <- makeLMintoFunction(pc = polyFit$coefficients, TRUE)

            # Find the derivative of the best fit polynomial
            gapeDeriv <- D(gapeBestFit, name = "x")

            #Find max gape velocity from frame 5 until the max gape is reached
            gapeVelocityNow <- eval(gapeDeriv, envir= list(x=2:timeMaxGape))
            #Divide frame rate calibration by 1000 so it's in units of mm/milliseconds
            finalData$gapeVel[i] <- max(gapeVelocityNow)*(finalData$frameRate[i]/1000)

            #Plot the gape velocity
            #plot(gapeVelocityNow)

          }
        } else {
          finalData$maxGape[i] <- NA
          finalData$timeMaxGape[i] <- NA
          finalData$gapeVel[i] <- NA
          skippedVids <- c(skippedVids,i)
          print(paste('Skipped video (gape): ', i))
        }


        ##############
        # Protrusion #

        #Run through each frame and find the protrusion (distance from nasal bone to upper jaw), save the max value and time of max
        tempVect <- rep(NA, length(dataNow[,1,1,1]))
        for (k in 1:length(dataNow[,1,1,1])) {
          xUJnow <- dataNow[k,1,1,1]
          yUJnow <- dataNow[k,1,2,1]

          xNasalnow <- dataNow[k,3,1,1]
          yNasalnow <- dataNow[k,3,2,1]

          protNowpix <- sqrt((xUJnow-xNasalnow)^2+(yUJnow-yNasalnow)^2)

          protNow <- protNowpix/finalData$cal_pix.mm[i]

          if(is.na(protNow)==FALSE) {
            if(protNow>=0) {
              tempVect[k] <- protNow
            }
          }
        }


        #This runs a polynomial fit on protrusion, then saves out the max protrusion from the fitted values

        if (sum(!is.na(tempVect))>10) {

          protNow <- na.omit(tempVect)
          framesNow <- seq(from=1, to=longestVid, by=1)[which(!is.na(tempVect))]


          #Create the first polynomial fit, and then get the positions that this line represents
          polyFit <- lm(protNow ~ poly(framesNow, 6,raw = TRUE))
          errorFitLine <- predict(polyFit)

          #plot(protNow~framesNow)
          #lines(errorFitLine~framesNow)

          #Save the polynomial fit data to protDF
          for(b in 1:length(errorFitLine)){
            protDF[framesNow[b],1,i] <- errorFitLine[b]
          }

          timeMaxprot <- which(protDF[,1,i]==max(protDF[,1,i], na.rm = TRUE))

          #Save the max best fit protrusion
          finalData$maxProt[i] <- max(protDF[,1,i], na.rm=TRUE)
          finalData$timeMaxProt[i] <- timeMaxprot


          #Now find the protrusion velocity
          #Only do this if the maximum protrusion falls outside the first 10 frames in the video
          if (timeMaxprot>9){
            #Make an expression to represent the best fit polynomial from lm()
            protBestFit <- makeLMintoFunction(pc = polyFit$coefficients, TRUE)

            # Find the derivative of the best fit polynomial
            protDeriv <- D(protBestFit, name = "x")

            #Find max protrusion velocity from frame 5 until the max protrusion is reached
            protVelocityNow <- eval(protDeriv, envir= list(x=2:timeMaxprot))
            #Divide frame rate calibration by 1000 so it's in units of mm/milliseconds
            finalData$protVel[i] <- max(protVelocityNow)*(finalData$frameRate[i]/1000)

            ##plot the protrusion velocity
            #plot(protVelocityNow)


          }
        } else {
          finalData$maxProt[i] <- NA
          finalData$timeMaxProt[i] <- NA
          finalData$protVel[i] <- NA
          skippedVids <- c(skippedVids,i)
          print(paste('Skipped video (protrusion): ', i))
        }








        #########
        # Hyoid #
        #Run through each frame and find the hyoid depression (distance from the eye to the hyoid), save the max value and time of max
        tempVect <- rep(NA, length(dataNow[,1,1,1]))
        for (k in 1:length(dataNow[,1,1,1])) {
          xEyenow <- dataNow[k,5,1,1]
          yEyenow <- dataNow[k,5,2,1]

          xHyoidnow <- dataNow[k,6,1,1]
          yHyoidnow <- dataNow[k,6,2,1]

          hyoidNowpix <- sqrt((xEyenow-xHyoidnow)^2+(yEyenow-yHyoidnow)^2)

          hyoidNow <- hyoidNowpix/finalData$cal_pix.mm[i]

          if(is.na(hyoidNow)==FALSE) {
            if(hyoidNow>=0) {
              tempVect[k] <- hyoidNow
            }
          }
        }


        #This runs a polynomial fit on hyoid depression, then saves out the max hyoid depression from the fitted values

        if (sum(!is.na(tempVect))>10) {

          hyoidNow <- na.omit(tempVect)
          framesNow <- seq(from=1, to=longestVid, by=1)[which(!is.na(tempVect))]


          #Create the first polynomial fit, and then get the positions that this line represents
          polyFit <- lm(hyoidNow ~ poly(framesNow, 6,raw = TRUE))
          errorFitLine <- predict(polyFit)

          #plot(hyoidNow~framesNow)
          #lines(errorFitLine~framesNow)

          #Save the polynomial fit data to hyoidDF
          for(b in 1:length(errorFitLine)){
            hyoidDF[framesNow[b],1,i] <- errorFitLine[b]
          }

          timeMaxHyoid <- which(hyoidDF[,1,i]==max(hyoidDF[,1,i], na.rm = TRUE))

          #Save the max best fit hyoid depression
          finalData$maxHyoidDep[i] <- max(hyoidDF[,1,i], na.rm=TRUE)
          finalData$timeMaxHyoidDep[i] <- timeMaxHyoid


          #Now find the hyoid depression velocity
          #Only do this if the maximum hyoid depression falls outside the first 10 frames in the video
          if (timeMaxHyoid>9){
            #Make an expression to represent the best fit polynomial from lm()
            hyoidBestFit <- makeLMintoFunction(pc = polyFit$coefficients, TRUE)

            # Find the derivative of the best fit polynomial
            hyoidDeriv <- D(hyoidBestFit, name = "x")

            #Find max hyoid depression velocity from frame 5 until the max hyoid depression is reached
            hyoidVelocityNow <- eval(hyoidDeriv, envir= list(x=2:timeMaxHyoid))
            #Divide frame rate calibration by 1000 so it's in units of mm/milliseconds
            finalData$hyoidVel[i] <- max(hyoidVelocityNow)*(finalData$frameRate[i]/1000)

            ##plot the hyoid depression velocity
            #plot(hyoidVelocityNow)

          }
        } else {
          finalData$maxHyoidDep[i] <- NA
          finalData$timeMaxHyoidDep[i] <- NA
          finalData$hyoidVel[i] <- NA
          skippedVids <- c(skippedVids,i)
          print(paste('Skipped video (hyoid): ', i))
        }







        #####################
        # Cranial elevation #

        # Not sure how to do this


        #maxElevation
        #timeMaxElevation
        #elevationVel









        #####################
        # Maxilla rotation #


        xNasalInit <- dataNow[1,3,1,1]
        yNasalInit <- dataNow[1,3,2,1]

        xmaxillaInit <- dataNow[1,4,1,1]
        ymaxillaInit <- dataNow[1,4,2,1]

        maxillaVectInit <- c(xNasalInit,yNasalInit)-c(xmaxillaInit,ymaxillaInit)

        lengthMaxillaVectInit <- sqrt(maxillaVectInit[1]^2+maxillaVectInit[2]^2)

        if (!is.na(lengthMaxillaVectInit)) {

          #Run through each frame and find the maxilla rotation (difference from first frame), save the max value and time of max
          tempVect <- rep(NA, length(dataNow[,1,1,1]))

          for (k in 1:length(dataNow[,1,1,1])) {
            xNasalnow <- dataNow[k,3,1,1]
            yNasalnow <- dataNow[k,3,2,1]


            xmaxillanow <- dataNow[k,4,1,1]
            ymaxillanow <- dataNow[k,4,2,1]

            maxillaVectNow <- c(xNasalnow,yNasalnow)-c(xmaxillanow,ymaxillanow)

            if (!is.na(maxillaVectNow[1]) & !is.na(maxillaVectNow[1])) {
              lengthMaxillaVectNow <- sqrt(maxillaVectNow[1]^2+maxillaVectNow[2]^2)

              angNow <- acos((maxillaVectInit %*% maxillaVectNow)/(abs(lengthMaxillaVectInit)*abs(lengthMaxillaVectNow)))

              maxAngNowDegree <- angNow*180/pi

              if(is.na(angNow)==FALSE) {
                if(angNow>=0) {
                  tempVect[k] <- maxAngNowDegree
                  maxRotDF[k,1,i] <- maxAngNowDegree
                  finalData$maxMaxillaRot [i] <- max(tempVect, na.rm = TRUE)
                  finalData$timeMaxMaxillaRot[i] <- match(max(tempVect, na.rm = TRUE),tempVect)
                }
              }
            }
          }


          if (sum(!is.na(tempVect))>10) {

            maxRotNow <- na.omit(tempVect)
            framesNow <- seq(from=1, to=longestVid, by=1)[which(!is.na(tempVect))]


            #Create the first polynomial fit, and then get the positions that this line represents
            polyFit <- lm(maxRotNow ~ poly(framesNow, 6,raw = TRUE))
            errorFitLine <- predict(polyFit)

            #plot(maxRotNow~framesNow)
            #lines(errorFitLine~framesNow)

            #Save the polynomial fit data to maxRotDF
            for(b in 1:length(errorFitLine)){
              maxRotDF[framesNow[b],1,i] <- errorFitLine[b]
            }

            timeMaxMaxillaRot <- which(maxRotDF[,1,i]==max(maxRotDF[,1,i], na.rm = TRUE))

            #Save the max best fit hyoid depression
            finalData$maxMaxillaRot [i] <- max(maxRotDF[,1,i], na.rm=TRUE)
            finalData$timeMaxMaxillaRot[i] <- timeMaxMaxillaRot


            #Now find the maxilla rotation velocity
            #Only do this if the maximum maxilla rotation falls outside the first 10 frames in the video
            if (timeMaxMaxillaRot>9){
              #Make an expression to represent the best fit polynomial from lm()
              maxillaBestFit <- makeLMintoFunction(pc = polyFit$coefficients, TRUE)

              # Find the derivative of the best fit polynomial
              maxillaDeriv <- D(maxillaBestFit, name = "x")

              #Find max hyoid depression velocity from frame 5 until the max hyoid depression is reached
              maxillaVelocityNow <- eval(maxillaDeriv, envir= list(x=2:timeMaxMaxillaRot))
              #Divide frame rate calibration by 1000 so it's in units of degrees/milliseconds
              finalData$maxillaVel[i] <- max(maxillaVelocityNow)*(finalData$frameRate[i]/1000)

              ##plot the hyoid depression velocity
              #plot(maxillaVelocityNow)


            }
          } else {
            finalData$maxMaxillaRot[i] <- NA
            finalData$timeMaxMaxillaRot[i] <- NA
            finalData$maxillaVel[i] <- NA
            skippedVids <- c(skippedVids,i)
            print(paste('Skipped video (maxilla rotation): ', i))
          }

        } else {
          finalData$maxMaxillaRot[i] <- NA
          finalData$timeMaxMaxillaRot[i] <- NA
          finalData$maxillaVel[i] <- NA
          skippedVids <- c(skippedVids,i)
          print(paste('Skipped video (maxilla rotation): ', i))
        }






        #####################
        # Mandible rotation #


        xLJInit <- dataNow[1,2,1,1]
        yLJInit <- dataNow[1,2,2,1]

        xmaxillaInit <- dataNow[1,4,1,1]
        ymaxillaInit <- dataNow[1,4,2,1]

        ljVectInit <- c(xLJInit,yLJInit)-c(xmaxillaInit,ymaxillaInit)

        lengthljVectInit <- sqrt(ljVectInit[1]^2+ljVectInit[2]^2)

        if (!is.na(lengthljVectInit)) {

          #Run through each frame and find the lj rotation (difference from first frame), save the max value and time of max
          tempVect <- rep(NA, length(dataNow[,1,1,1]))

          for (k in 1:length(dataNow[,1,1,1])) {
            xNasalnow <- dataNow[k,3,1,1]
            yNasalnow <- dataNow[k,3,2,1]


            xljnow <- dataNow[k,4,1,1]
            yljnow <- dataNow[k,4,2,1]

            ljVectNow <- c(xNasalnow,yNasalnow)-c(xljnow,yljnow)

            if (!is.na(ljVectNow[1]) & !is.na(ljVectNow[1])) {
              lengthljVectNow <- sqrt(ljVectNow[1]^2+ljVectNow[2]^2)

              angNow <- acos((ljVectInit %*% ljVectNow)/(abs(lengthljVectInit)*abs(lengthljVectNow)))

              maxAngNowDegree <- angNow*180/pi

              if(is.na(angNow)==FALSE) {
                if(angNow>=0) {
                  tempVect[k] <- maxAngNowDegree
                  mandRotDF[k,1,i] <- maxAngNowDegree
                  finalData$maxMandRot [i] <- max(tempVect, na.rm = TRUE)
                  finalData$timeMaxMandRot[i] <- match(max(tempVect, na.rm = TRUE),tempVect)
                }
              }
            }
          }


          if (sum(!is.na(tempVect))>10) {

            maxRotNow <- na.omit(tempVect)
            framesNow <- seq(from=1, to=longestVid, by=1)[which(!is.na(tempVect))]


            #Create the first polynomial fit, and then get the positions that this line represents
            polyFit <- lm(maxRotNow ~ poly(framesNow, 6,raw = TRUE))
            errorFitLine <- predict(polyFit)

            #plot(maxRotNow~framesNow)
            #lines(errorFitLine~framesNow)

            #Save the polynomial fit data to mandRotDF
            for(b in 1:length(errorFitLine)){
              mandRotDF[framesNow[b],1,i] <- errorFitLine[b]
            }

            timeMaxljRot <- which(mandRotDF[,1,i]==max(mandRotDF[,1,i], na.rm = TRUE))

            #Save the max best fit hyoid depression
            finalData$maxMandRot [i] <- max(mandRotDF[,1,i], na.rm=TRUE)
            finalData$timeMaxMandRot[i] <- timeMaxljRot


            #Now find the lj rotation velocity
            #Only do this if the maximum lj rotation falls outside the first 10 frames in the video
            if (timeMaxljRot>9){
              #Make an expression to represent the best fit polynomial from lm()
              ljBestFit <- makeLMintoFunction(pc = polyFit$coefficients, TRUE)

              # Find the derivative of the best fit polynomial
              ljDeriv <- D(ljBestFit, name = "x")

              #Find max hyoid depression velocity from frame 5 until the max hyoid depression is reached
              ljVelocityNow <- eval(ljDeriv, envir= list(x=2:timeMaxljRot))
              #Divide frame rate calibration by 1000 so it's in units of degrees/milliseconds
              finalData$mandVel[i] <- max(ljVelocityNow)*(finalData$frameRate[i]/1000)

              ##plot the hyoid depression velocity
              #plot(ljVelocityNow)




            }
          } else {
            finalData$maxMandRot[i] <- NA
            finalData$timeMaxMandRot[i] <- NA
            finalData$mandVel[i] <- NA
            skippedVids <- c(skippedVids,i)
            print(paste('Skipped video (LJ rotation): ', i))
          }

        } else {
          finalData$maxMandRot[i] <- NA
          finalData$timeMaxMandRot[i] <- NA
          finalData$mandVel[i] <- NA
          skippedVids <- c(skippedVids,i)
          print(paste('Skipped video (LJ rotation): ', i))
        }









        ################
        # Ram velocity #

        #First find the earliest frame where the eye is visible
        f <- 1
        firstEyeFrame <- 1
        while (f < length(dataNow[,1,1,1])) {

          if (is.na(dataNow[f,5,1,])) {
            f <- f+1
            firstEyeFrame <- f
          } else {

            xEyeInit <- dataNow[f,5,1,1]
            yEyeInit <- dataNow[f,5,2,1]

            f <- length(dataNow[,1,1,1]+1)
          }



        }

        #Now calculate the average velocity between the first frame and the time when they reach max gape
        #Only do this if the eye is visible before maax gape is reached
        if (!is.na(timeMaxGape)) {
          if (timeMaxGape>firstEyeFrame){

            xEyeFinal <- dataNow[timeMaxGape,5,1,1]
            yEyeFinal <- dataNow[timeMaxGape,5,2,1]

            eyeDistpix <- sqrt((xEyeFinal-xEyeInit)^2+(yEyeFinal-yEyeInit)^2)

            eyeDist <- eyeDistpix/finalData$cal_pix.mm[i]
            eyeTime <- timeMaxGape-firstEyeFrame

            ramVelNow <- eyeDist/eyeTime*(finalData$frameRate[i]/1000)

            finalData$ramVel[i] <- ramVelNow

            print(paste('Done with video: ', i))
          } else {
            finalData$ramVel[i] <- NA
            skippedVids <- c(skippedVids,i)
            print(paste('Skipped video (ram velocity): ', i))
          }
        } else {
          finalData$ramVel[i] <- NA
          skippedVids <- c(skippedVids,i)
          print(paste('Skipped video (ram velocity): ', i))
        }






        ##########################
        # Kinematic simultaneity #

        #Measure how well aligned the different kinematic timings are
        #Using the coefficient of variation

        if(!is.na(timeMaxGape) & !is.na(timeMaxprot) & !is.na(timeMaxHyoid)) {

          timingVect <- c(timeMaxGape,timeMaxprot, timeMaxHyoid)
          timeVariance <- sd(timingVect)/mean(timingVect)
          #This is technically called the coefficient of variation

          finalData$simultKinematics[i] <- timeVariance
          finalData$simultKinematics2[i] <- timeMaxGape-timeMaxHyoid

        } else {

          finalData$simultKinematics[i] <- NA

        }

        h5closeAll()


        protPlotNowCali <- protDF[,,i]/finalData$cal_pix.mm[i]
        protTimeCali <- framesNow[which(!is.na(protPlotNowCali))]/finalData$frameRate[i]
        gapePlotNowCali <- gapeDF[,,i]/finalData$cal_pix.mm[i]
        gapeTimeCali <- framesNow[which(!is.na(gapePlotNowCali))]/finalData$frameRate[i]
        hyoidPlotNowCali <- hyoidDF[,,i]/finalData$cal_pix.mm[i]
        hyoidTimeCali <- framesNow[which(!is.na(hyoidPlotNowCali))]/finalData$frameRate[i]
        maxRotPlotNowCali <- maxRotDF[,,i]/finalData$cal_pix.mm[i]
        maxRotTimeCali <- framesNow[which(!is.na(maxRotPlotNowCali))]/finalData$frameRate[i]
        maxTimePlotNowCali <- max(framesNow)/ finalData$frameRate[i]

        #This plots the video, showing the relative timing of each kinematic variable and how good the polynomial fit is

        #png() and dev.off() are the starting and ending calls to save it as a png on the disk.

        fullFilenameNow <- paste0('C:/Users/matth/Documents/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data/Sleap data/Plots/',fileNameNow,'.png')
        png(fullFilenameNow,
            width     = 5,
            height    = 7,
            units     = "in",
            res       = 300)


        par(mfrow = c(4, 1))

        if(sum(protDF[,,i], na.rm=TRUE)>0){
          plot(na.omit(protPlotNowCali)~protTimeCali, main='Protrusion', xlim=c(0,maxTimePlotNowCali), col='grey')
          lines(y=na.omit(protPlotNowCali), x=protTimeCali, lwd=2, col='red')
          abline(v=finalData$timeMaxProt[i]/finalData$frameRate[i])
        } else{plot.new()}

        if(sum(gapeDF[,,i], na.rm=TRUE)>0){
          plot(na.omit(gapePlotNowCali)~gapeTimeCali, main='Gape', xlim=c(0,maxTimePlotNowCali), col='grey')
          lines(y=na.omit(gapePlotNowCali), x=gapeTimeCali, lwd=2, col='red')
          abline(v=finalData$timeMaxGape[i]/finalData$frameRate[i])
        } else{plot.new()}

        if(sum(hyoidDF[,,i], na.rm=TRUE)>0){
          plot(na.omit(hyoidPlotNowCali)~hyoidTimeCali, main='Hyoid Depression', xlim=c(0,maxTimePlotNowCali), col='grey')
        #   lines(y=na.omit(hyoidPlotNowCali), x=hyoidTimeCali, lwd=2, col='red')
        #   abline(v=finalData$timeMaxHyoidDep[i]/finalData$frameRate[i])
        # } else{plot.new()}
        # 
        # if(sum(maxRotDF[,,i], na.rm=TRUE)>0){
        #   plot(na.omit(maxRotPlotNowCali)~maxRotTimeCali, main='Maxilla Rotation', xlim=c(0,maxTimePlotNowCali), col='grey')
        #   lines(y=na.omit(maxRotPlotNowCali), x=maxRotTimeCali, lwd=2, col='red')
        #   abline(v=finalData$timeMaxMaxillaRot[i]/finalData$frameRate[i])
        # } else{plot.new()}
        # 
        # 
        # mtext(fileNameNow, side = 1, line = -2, adj=.9,outer = TRUE)
        # 
        # dev.off()
        
        }
}
  
  
  
  
  
  
  ############################     CALCULATED VARIABLES      ############################
  
  #Relative_ttpg, Relative_tmax, Kinetic_Synchronization
  
  
  
  

finalData$species <- vector(mode="character", length=nrow(finalData))
for (i in 1:nrow(finalData)) {
  if (substring(finalData$Animal.ID[i],1,2)=="AC") {
    finalData$species[i]="AC"
  } else if (substring(finalData$Animal.ID[i],1,2)=="TR") {
    finalData$species[i]="TRC"
  } else {
    finalData$species[i]=substring(finalData$Generation[i],1,2)
  }
}

finalData$family[which(finalData$family=="TRC")] <- "F0"

for (i in 1:nrow(finalData)) {
  if (finalData$Generation[i] == "F1.2" | finalData$Generation[i] == "F1.3" ) {
    finalData$Generation[i] = "F1"
  } 
}

# Save the data as a CSV, analyze in separate script

write.csv(finalData, file='C:/Users/matth/Documents/Science/Postdoc/Albertson lab/Projects/Hybrid Feeding/ACxTRC/Data/Sleap data/kinematicsAll_sleapOutput.csv')


