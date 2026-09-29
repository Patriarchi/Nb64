#Open Field analysis code for the manuscript
#An engineered nanobody inhibitor for molecular-to-circuit control of opioid receptor function
#Figure 5 k-m

#For those unfamiliar with DeepLab Cut: The output is a file in which each keypoint (detected body part) has three assigned columns, where two columns are the predicted
#x and y coordinates and the third is a likelihood, which is an estimate of how sure the model is that the keypoint is at this location.

rm(list = ls())
library(tidyverse)
library(plotly)
library(lattice)
library(htmlwidgets)
library(cowplot)
library(plotrix)
library(viridis)

#set working directory with data (will also be output folder)
setwd("K:\\...")


#Definition of experimental parameters

# Number of detected keypoints (bodyparts) in DeepLab Cut
keypoints <- 11
# The bodypart keypoints used in this analysis are
bodylist <- list("nose", "headcenter", "bodycenter", "tailbase", "tailtip", "fiberbody", "fiberconnector")
fps <- 25
# The corner names must correspond to their name in the DLC file
corners <- c("TopLeftCorner", "TopRightCorner", "BottomLeftCorner", "BottomRightCorner")
# Dimensions of the behavioural arena in cm
arenax <- 50
arenay <- 50
# Dimensions of what is considered the center of the arena in cm
centerx <- 25
centery <- 25
# Threshold of acceptance of likelihood estimated by DLC. Any keypoint that has a lower likelihood will be flagged and interpolated.
likelihood_thresh <- 0.85
# Number of Frames to be binned together
movav <- 4
# Start of baseline measurement
baseline_start <- 0.5 #mins
# Total length of experiment
length_experiment <- 35 #including baseline, in mins

# Injection window in minutes for exclusion
inj_time_start <- 4
inj_time_end <- 7

# Create two 3x4 matrices with the name of the corner (upperleft, upperright, lowerright, lowerleft) and the x-y coordinates of each (0,0) (0, arenay) (arenax, 0) (arenax, arenay)
# The first matrix (arenasize) defines the actual size of the arena and centerarenasize defines what is considered the center of the arena, in case one wants to 
# quantify number of center entries/time spent in the center of the arena
arenasize <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 3))
centerarenasize <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 3))

for(i in 1:length(arenasize[[1]])){
  if(i == 1){
    arenasize[i, 1] <- 0
    arenasize[i, 2] <- 0
    arenasize[i, 3] <- "upperleft"
    
    centerarenasize[i, 1] <- 0+centerx/2
    centerarenasize[i, 2] <- 0+centery/2
    centerarenasize[i, 3] <- "upperleft"
  }else if(i == 2){
    arenasize[i, 1] <- arenax
    arenasize[i, 2] <- 0
    arenasize[i, 3] <- "upperright"
    
    centerarenasize[i, 1] <- arenax-centerx/2
    centerarenasize[i, 2] <- 0-centery/2
    centerarenasize[i, 3] <- "upperright"
  }else if(i == 3){
    arenasize[i, 1] <- 0 
    arenasize[i, 2] <- arenay
    arenasize[i, 3] <- "lowerleft"
    
    centerarenasize[i, 1] <- 0+centerx/2
    centerarenasize[i, 2] <- arenay-centery/2
    centerarenasize[i, 3] <- "lowerleft"
  }else{
    arenasize[i, 1] <- arenax
    arenasize[i, 2] <- arenay
    arenasize[i, 3] <- "lowerright"
    
    centerarenasize[i, 1] <- arenax - centerx/2
    centerarenasize[i, 2] <- arenay - centery/2
    centerarenasize[i, 3] <- "lowerright"
  }
}

#naming convention: date(yyyymmdd)_mousecode(X_Ex/Cn_1 etc)_otherexperimentalinfo

# Data structure: nested lists, where the umbrella list is df.list, containing sublists for each file that is read in (each recording). 

#The first position of the sublist will always be a data frame that includes all the output from DLC and the experimental info such as
#compund that was given, dose, animal ID etc. (df.list[[1]][[1]] = the df from the first file of the directory)
#The second position is a list of Histograms of likelihoods per detected keypoint
#the third position is a data frame with the coordinates of the four corners
#the fourth position is the homography transformation matrix
#the fifth position will be a data frame (df) but with the transformed locations
#the sixth position will be a df but only with a single point representing the mouse (a specific body part is chosen, generally the one that is most consistently detected)
#the seventh position will be a 2D trace plot
#the eighth position will be a df containing binned data for baseline
#the ninth will be a small df with summary statistics about the baseline
#the tenth will be the trajectory plot (tornado) where each timepoint is a dot
#the 11th will be the trajectory plot (tornado) but with lines
#the 12th-16th the same as 7-11 but for post injection
#the 17th position is a df that contains each data point before binning (for the bin df, position 8 and 13)


# Read in files
file.list <- list.files(pattern='*.csv')

df.list <- list()

for(i in 1:length(file.list)){
  sublist <- list()
  df <- read.csv(file.list[i])
  sublist[[1]] <- df
  df.list[[i]] <- sublist
}
rm(sublist)
rm(df)


#=============================== Housekeeping ==================================

# This section exists to clean up the column names, add experimental information and to convert the time vector to minutes

# Fix column names
colnamels <- list()
listcounter <- 0

for(i in 1:length(df.list[[1]][[1]])){
  listcounter <- listcounter + 1
  if(listcounter == 1){
    colnamels[listcounter] <- "Time"
  }else{
    name <- paste(df.list[[1]][[1]][1, i], df.list[[1]][[1]][2, i], sep = "_")
    colnamels[listcounter] <- name
  }
}

ind <- keypoints*3 +2
colnamels[ind:(ind+4)] <- c("ID", "Group", "Cpd", "Day", "other")

for(i in 1:length(df.list)){
  df.list[[i]][[1]][, ind:(ind+4)] <- NA
  colnames(df.list[[i]][[1]]) <- colnamels
  df.list[[i]][[1]] <- tail(df.list[[i]][[1]], -2)
  df.list[[i]][[1]] <- df.list[[i]][[1]] %>% mutate_if(is.character, as.numeric)
}

# Add the mouse code, experimental group, injection received to the df
for(i in 1:length(file.list)){
  info <- str_split(file.list[[i]], "_")
  #mouse code
  df.list[[i]][[1]][[keypoints*3 + 2]] <- paste(info[[1]][2], info[[1]][3], info[[1]][4], sep = "_")
  #con or ex
  df.list[[i]][[1]][[keypoints*3 + 3]] <- info[[1]][3]
  df.list[[i]][[1]][[keypoints*3 + 4]] <- info[[1]][5]
  df.list[[i]][[1]][[keypoints*3 + 5]] <- info[[1]][6]
}

# Creating a time vector (the output from DLC is a count of frames)
for(i in 1:length(df.list)){
  df.list[[i]][[1]] <- mutate(df.list[[i]][[1]], Time = Time/fps, mins = Time*0.0167)
}

colnamels[[40]] <- "mins"

# Checkpoint: Print length of the recordings
for(i in 1:length(df.list)){
  print(paste("recording length", df.list[[i]][[1]]$ID[[1]], df.list[[i]][[1]]$Cpd[[1]], df.list[[i]][[1]]$Dose[[1]],
              df.list[[i]][[1]]$other[[1]], ":", max(df.list[[i]][[1]]$mins)))
}

#========================== Likelihood histograms ==============================

for(i in 1:length(df.list)){
  histograms <- list()
  timeseries <- list()
  index <- 1
  for(j in 1:length(df.list[[i]][[1]])){
    if(grepl("likelihood", colnamels[j], fixed = T)){
      
      histog <- ggplot(df.list[[i]][[1]], aes(x = as.numeric(df.list[[i]][[1]][, j]))) +
        geom_histogram()+
        theme_minimal()+
        xlab("likelihood")+
        ggtitle(paste("likelihood of coordinates of keypoint: ", colnamels[[j]]))+
        labs(subtitle = paste(df.list[[i]][[1]]$ID, df.list[[i]][[1]]$Cpd, df.list[[i]][[1]]$Dose, "mg/kg"))
      
      histograms[[index]] <- histog
      index <- index + 1
      
    }
  }
  df.list[[i]][[2]] <- histograms
}

rm(histog)
rm(histograms)

#========================== Checking the corners ===========================

# The corners are also detected for every frame but assuming they are relatively stable across time, the median of the detected coordinates will be used as the true corner location

for(i in 1:length(df.list)){
  corner_coords <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 3))
  for(k in 1:length(corners)){
    for(j in 1:length(df.list[[i]][[1]])){
      if(colnamels[[j]] == paste(corners[[k]], "x", sep = "_")){
        corner_coords[k, 1] <- median(df.list[[i]][[1]][, j])
        corner_coords[k, 2] <- median(df.list[[i]][[1]][, j+1])
        corner_coords[1, 3] <- "upperleft"
        corner_coords[2, 3] <- "upperright"
        corner_coords[3, 3] <- "lowerleft"
        corner_coords[4, 3] <- "lowerright"
      }
    }
  }
  df.list[[i]][[3]] <- corner_coords
}


for(i in 1:length(df.list)){
  if(i == 1){
    cornerdf <- df.list[[i]][[3]]
    cornerdf[[4]] <- paste(df.list[[i]][[1]]$ID, df.list[[i]][[1]]$Group, df.list[[i]][[1]]$Cpd, df.list[[i]][[1]]$Dose, sep = "_")[[1]]
  }else{
    placeholder <- df.list[[i]][[3]]
    placeholder[[4]] <- paste(df.list[[i]][[1]]$ID, df.list[[i]][[1]]$Group, df.list[[i]][[1]]$Cpd, df.list[[i]][[1]]$Dose, sep = "_")[[1]]
    cornerdf <- rbind(cornerdf, placeholder)
  }
}

colnames(cornerdf) <- c("x_coord", "y_coord", "corner", "recording")

# Checkpoint: plot corner points - they need not be exactly squared since they should still capture the distortion coming from the camera

ggplot(cornerdf, aes(x=x_coord, y = y_coord, colour = recording))+
  geom_point() +
  ggtitle("Median cornerpoints")+
  xlab("x") +
  ylab("y") +
  theme_minimal()


#===================== Homography matrix calculation ===========================

# Based on the detected coordinates of the corner points and the "real" coordinates (size of the behavioural arena), the
# distortion of the camera can be inferred. This is caputred by a homography matrix, which can then be applied to all other 
# detected keypoints to have an undistorted movement trace 

for(i in 1:length(df.list)){
  A <- matrix(data = c(df.list[[i]][[3]][1,1], df.list[[i]][[3]][1,2], 1, 0, 0, 0, -arenasize[1, 1]*df.list[[i]][[3]][1,1], -arenasize[1, 1]*df.list[[i]][[3]][1,2], -arenasize[1, 1],
                       0, 0, 0, df.list[[i]][[3]][1,1], df.list[[i]][[3]][1,2], 1, -arenasize[1, 2]*df.list[[i]][[3]][1,1], -arenasize[1, 2]*df.list[[i]][[3]][1,2], -arenasize[1, 2],
                       
                       df.list[[i]][[3]][2,1], df.list[[i]][[3]][2,2], 1, 0, 0, 0, -arenasize[2, 1]*df.list[[i]][[3]][2,1], -arenasize[2, 1]*df.list[[i]][[3]][2,2], -arenasize[2, 1],
                       0, 0, 0, df.list[[i]][[3]][2,1], df.list[[i]][[3]][2,2], 1, -arenasize[2, 2]*df.list[[i]][[3]][2,1], -arenasize[2, 2]*df.list[[i]][[3]][2,2], -arenasize[2, 2],
                       
                       df.list[[i]][[3]][3,1], df.list[[i]][[3]][3,2], 1, 0, 0, 0, -arenasize[3, 1]*df.list[[i]][[3]][3,1], -arenasize[3, 1]*df.list[[i]][[3]][3,2], -arenasize[3, 1],
                       0, 0, 0, df.list[[i]][[3]][3,1], df.list[[i]][[3]][3,2], 1, -arenasize[2, 2]*df.list[[i]][[3]][3,1], -arenasize[3, 2]*df.list[[i]][[3]][3,2], -arenasize[3, 2],
                       
                       df.list[[i]][[3]][4,1], df.list[[i]][[3]][4,2], 1, 0, 0, 0, -arenasize[4, 1]*df.list[[i]][[3]][4,1], -arenasize[4, 1]*df.list[[i]][[3]][4,2], -arenasize[4, 1],
                       0, 0, 0, df.list[[i]][[3]][4,1], df.list[[i]][[3]][4,2], 1, -arenasize[4, 2]*df.list[[i]][[3]][4,1], -arenasize[4, 2]*df.list[[i]][[3]][4,2], -arenasize[4, 2]),
              nrow = 8, ncol = 9, byrow = T)
  
  lambda <- min(eigen(t(A)%*%A)[[1]])
  ind <- which(eigen(t(A)%*%A)[[1]] == lambda)
  
  homography <- matrix(data = eigen(t(A)%*%A)[[2]][, ind], nrow = 3, ncol = 3, byrow = T)
  
  df.list[[i]][[4]] <- homography
}

#================ Application of Homography to other keypoints =================

colnums <- seq(from=2, to=1 + (keypoints-4)*3, by=3)

for(i in 1:length(df.list)){
  df.list[[i]][[5]]<- df.list[[i]][[1]]
  transformation_points <- keypoints - length(corners)
  for(k in colnums){
    srcmatrix <- matrix(data = 1, nrow = length(df.list[[i]][[5]][[k]]), ncol = 3)
    srcmatrix[, 1] <- df.list[[i]][[5]][[k]]
    srcmatrix[, 2] <- df.list[[i]][[5]][[k+1]]
    
    srcmatrix <- t(srcmatrix)
    
    tilde <- df.list[[i]][[4]] %*% srcmatrix
    tilde[1, ] <- tilde[1,]/tilde[3, ]
    tilde[2, ] <- tilde[2,]/tilde[3, ]
    
    tilde <- t(tilde)
    
    df.list[[i]][[5]][[k]] <- tilde[, 1]
    df.list[[i]][[5]][[k+1]] <- tilde[, 2]
  }
}


#================= Interpolation of low likelihood points ======================

# Steps:
# A column is added to specify experimental phase (baseline, injection, post-injection)
# A keypoint is picked as the representative "mouse point", in this case, BodyCenter
# Loop through each frame and check likelihood of keypoint, flag if likelihood < threshold you set. Keep track of the frames (index) that have low likelihood (indexlist = list containing the first index of a block (may also have length = 1) of low likelihood predictions and countlist = size of low likelihood block starting on index [i]) 
# Interpolate; If the first values at the beginning of the video have low likelihood, they will be replaced by the first high likelihood frame in the df,
#              If the block of low likelihood values is somewhere in the middle of the df, the euclidian distance between the last and next high likelihood frame will be calculated and divided by the number of low likelihood frames in the block. The mouse is assumed to move in equal steps from the last to next high likelihood frame,
#              If it is at the end of the df, then the last high likelihood frame will be used
# If there are any remaining points lying outside of the behavioural arena, they are constrained to the border
# There are some instances of teleportation, which nevertheless have high likelihood. These cases are flagged and interpolated in the same way as described above

for(i in 1:length(df.list)){
  if(max(df.list[[i]][[5]]$mins) > length_experiment){
    baseline <- filter(df.list[[i]][[1]], mins < inj_time_start)
    injection <- filter(df.list[[i]][[1]], mins >= inj_time_start & mins < inj_time_end)
    #you can either introduce an end cutoff or not here.
    post <- filter(df.list[[i]][[1]], mins >= inj_time_end & mins < length_experiment)
    #this serves as control whether all the datapoints were assigned to a phase of the experiment; with a cutoff at the end, this will not be the same.
    print(paste("Total len", length(df.list[[i]][[5]][[1]]), "sum", sum(length(baseline$Time),length(injection$Time), length(post$Time))))
    
    df.list[[i]][[5]][1:length(baseline[[1]]),length(colnamels) + 1] <- "baseline"
    df.list[[i]][[5]][(length(baseline[[1]])+1):(length(baseline[[1]])+length(injection[[1]])), length(colnamels)+1] <- "injection"
    df.list[[i]][[5]][(length(baseline[[1]])+length(injection[[1]])+1):(length(baseline[[1]])+length(injection[[1]])+length(post[[1]])), length(colnamels)+1] <- "post"
    
  }else{
    print(paste("The recording", df.list[[i]][[5]]$ID[[1]], df.list[[i]][[5]]$Cpd[[1]], df.list[[i]][[5]]$Dose[[1]], df.list[[i]][[5]]$other[[1]], "is not long enough"))
    label <- readline(prompt = "how do you want to label this data?")
    df.list[[i]][[5]][, length(colnamels) + 1] <- label
    
  }
}

colnamels[[length(colnamels)+1]] <- "phase"

# Generally speaking, it makes sense to follow any of the center points; implant/shoulders/center body
# not the tail because it's likely occluded quite often

for(i in 1:length(df.list)){
  colnames(df.list[[i]][[5]]) <- colnamels
  
  mouse <- as.data.frame(matrix(data = NA, nrow = length(filter(df.list[[i]][[5]], !is.na(phase))[[1]]), ncol = 5))
  
  mouse[, 1] <- filter(df.list[[i]][[5]], !is.na(phase))$BodyCenter_x
  mouse[, 2]<- filter(df.list[[i]][[5]], !is.na(phase))$BodyCenter_y
  mouse[, 3]<- filter(df.list[[i]][[5]], !is.na(phase))$BodyCenter_likelihood
  mouse[, 4]<- filter(df.list[[i]][[5]], !is.na(phase))$mins
  mouse[, 5]<- filter(df.list[[i]][[5]], !is.na(phase))$phase
  
  df.list[[i]][[6]] <- mouse
  
}


for(a in 1:length(df.list)){
  # list that keeps track of the first position of an NA block, needs to be added when a new block starts
  indexlist <- list()
  # keeps track of block size; should always be added to at the end of a block
  countlist <- list()
  # index of the indexlist & countlist; should increase by 1 every time the counter is set to 1
  list_pos <- 0
  # counts number of successive NAs; should become 1 every time a new block starts, the actual number that goes into the countlist
  counter <- 0
  
  indexlist[[1]] <- 5 #this is a bit of an awkward workaround, please ignore
  
  for(j in 1:length(df.list[[a]][[6]][[1]])){
    #if the phase is not injection and the likelihood of the detected point is lower than the set threshold
    if(df.list[[a]][[6]][j, 5] != "injection" & df.list[[a]][[6]][j, 3] < likelihood_thresh){
      
      #then this should be tracked by the counter and the indexlist, such that it can be interpolated later.
      #low likelihood points can be at different positions of the df (first, last or any in between) and depending on this, the procedure of calculation for interpolation differs a bit
      #secondly, low likelihood points could be followed by another low likelihood point and/or preceded by one, or be the only low likelihood point between higher likelihood points 
      
      if(j == 1){ #if the very first point has low likelihood
        counter <- 1
        list_pos <- j
        # print("first")
        
        if(df.list[[a]][[6]][j+1, 3] >= likelihood_thresh){ #if the next one has high likelihood
          df.list[[a]][[6]][j, 6] <- df.list[[a]][[6]][j+1, 1]
          df.list[[a]][[6]][j, 7] <- df.list[[a]][[6]][j+1, 2]
          df.list[[a]][[6]][j, 8] <- "uncertain_start_single"
          
          countlist[[list_pos]] <- counter
          indexlist[[list_pos]] <- j
          
          # print("last")
          
          
        }else{
          #we do nothing because we want to keep counting the subsequent NAs
          #but we want to save the starting index
          indexlist[[list_pos]] <- j
        }
        
        
      }else if(j == length(df.list[[a]][[6]][[1]])){ #if the very last point has low likelihood
        
        if(df.list[[a]][[6]][j-1, 3] >= likelihood_thresh){
          df.list[[a]][[6]][j, 6] <- df.list[[a]][[6]][j-1, 1]
          df.list[[a]][[6]][j, 7] <- df.list[[a]][[6]][j-1, 2]
          df.list[[a]][[6]][j, 8] <- "uncertain_end_single"
          
          counter <- 1
          list_pos <- list_pos + 1
          countlist[[list_pos]] <- counter
          indexlist[[list_pos]] <- j
          
          # print("first")
          
        }else{
          counter <- counter + 1
          countlist[[list_pos]] <- counter
          # print("last")
        }
        
      }else{ #if j is neither 1 nor -1
        
        
        if(list_pos == 1 & indexlist[[1]] == 1){ #if this is a continuation of an ongoing counter increase starting at the beginning of df
          if(length(countlist) == 0 & df.list[[a]][[6]][j+1, 3] >= likelihood_thresh){ #if this is the last row with low likelihood in a block that started from the beginning
            counter <- counter + 1
            countlist[[list_pos]] <- counter
            list_pos <- list_pos + 1
            
            # print("last")
          }else if(length(countlist) == 0 & df.list[[a]][[6]][j+1, 3] < likelihood_thresh){
            counter <- counter + 1
          }
          
        }else if(indexlist[[1]] != 1){ #if this is a case in which there is not a low likelihood at the start of df
          
          if(df.list[[a]][[6]][j-1, 3] >= likelihood_thresh | !is.na(df.list[[a]][[6]][j-1, 6])){ #if it's a new block
            
            
            counter <- 1
            list_pos <- list_pos + 1
            indexlist[[list_pos]] <- j
            
            # print("first")
            
            if(df.list[[a]][[6]][j+1, 3] >= likelihood_thresh){
              df.list[[a]][[6]][j, 6] <- (df.list[[a]][[6]][j-1, 1]+ df.list[[a]][[6]][j+1, 1])/2
              df.list[[a]][[6]][j, 7] <- (df.list[[a]][[6]][j-1, 2]+ df.list[[a]][[6]][j+1, 2])/2
              df.list[[a]][[6]][j, 8] <- "single_interpol"
              
              countlist[[list_pos]] <- counter
              # print("last")
              
            }else{ #if it's the first of a block of values
              
            }
            
            
          }else{ #if it's not a new block
            
            counter <- counter + 1
            
            
            if(df.list[[a]][[6]][j+1, 3] >= likelihood_thresh | df.list[[a]][[6]][j+1, 5] == "injection"){
              countlist[[list_pos]] <- counter
              # print("last")
            }else{ #if it's also not the last of the block
              #just increasing the counter is fine
            }
            
          }
          
        }else if(indexlist[[1]] == 1 & list_pos > 1){ #if the df started with low likelihood but that first counter was already "closed"
          
          if(df.list[[a]][[6]][j-1, 3] >= likelihood_thresh | !is.na(df.list[[a]][[6]][j-1, 6])){ #if it's a new block
            
            
            counter <- 1
            
            if(list_pos == 2 & length(indexlist) < 2){
              list_pos <- 2
              indexlist[[list_pos]] <- j
            }else{
              list_pos <- list_pos + 1
              indexlist[[list_pos]] <- j
            }
            
            
            
            
            # print("first")
            
            
            if(df.list[[a]][[6]][j+1, 3] >= likelihood_thresh){
              df.list[[a]][[6]][j, 6] <- (df.list[[a]][[6]][j-1, 1]+ df.list[[a]][[6]][j+1, 1])/2
              df.list[[a]][[6]][j, 7] <- (df.list[[a]][[6]][j-1, 2]+ df.list[[a]][[6]][j+1, 2])/2
              df.list[[a]][[6]][j, 8] <- "single_interpol"
              
              countlist[[list_pos]] <- counter
              
              # print("last")
              
            }else{ #if it's the first of a block of values
              
            }
            
            
          }else{ #if it's not a new block
            
            
            counter <- counter + 1
            
            if(df.list[[a]][[6]][j+1, 3] >= likelihood_thresh | df.list[[a]][[6]][j+1, 5] == "injection"){
              
              countlist[[list_pos]] <- counter
              # print("last")
            }else{ #if it's also not the last of the block
              #just increasing the counter is fine
            }
            
          }
          
        }
        
      }
      
      
    }else{
      df.list[[a]][[6]][j, 6] <- df.list[[a]][[6]][j, 1]
      df.list[[a]][[6]][j, 7] <- df.list[[a]][[6]][j, 2]
      df.list[[a]][[6]][j, 8] <- "unchanged"
    }
    
    
  }
  
  #this is to check that the two lists have the same length and also to see how far we've progressed
  print(a)
  # print(length(indexlist))
  # print(length(countlist))
  
  if(length(indexlist) == 1 & indexlist[[1]] == 5 & length(countlist) == 0){
    
  }else{
    
    
    
    for(i in 1:length(indexlist)){
      #if you have a block of NAs
      if(countlist[[i]] != 1){
        
        # and if it includes the first value of the data frame, then the x and y values will be equal to the first high likelihood point
        if(i == 1 & indexlist[[i]] == 1){
          df.list[[a]][[6]][i:countlist[[i]], 6] <- df.list[[a]][[6]][i+countlist[[i]], 1]
          df.list[[a]][[6]][i:countlist[[i]], 7] <- df.list[[a]][[6]][i+countlist[[i]], 2]
          df.list[[a]][[6]][i:countlist[[i]], 8] <- "first_highprob"
          # if it is the last couple of frames
        }else if((indexlist[[i]] + countlist[[i]] - 1) == length(df.list[[a]][[6]][[1]])){
          df.list[[a]][[6]][indexlist[[i]]:(countlist[[i]] + indexlist[[i]] - 1), 6] <- df.list[[a]][[6]][indexlist[[i]]-1, 1]
          df.list[[a]][[6]][indexlist[[i]]:(countlist[[i]] + indexlist[[i]] - 1), 7] <- df.list[[a]][[6]][indexlist[[i]]-1, 2]
          df.list[[a]][[6]][indexlist[[i]]:(countlist[[i]] + indexlist[[i]] - 1), 8] <- "last_highprob"
          # if it is not the first couple of frames
        }else{
          
          
          y_diff <- as.numeric(df.list[[a]][[6]][indexlist[[i]]+countlist[[i]], 7]-as.numeric(df.list[[a]][[6]][indexlist[[i]]-1, 7]))/(countlist[[i]]+1)
          x_diff <- as.numeric(df.list[[a]][[6]][indexlist[[i]]+countlist[[i]], 6]-as.numeric(df.list[[a]][[6]][indexlist[[i]]-1, 6]))/(countlist[[i]]+1)
          
          
          for(k in 1:(countlist[[i]])){
            # print(k)
            df.list[[a]][[6]][indexlist[[i]]-1+k, 6] <- as.numeric(df.list[[a]][[6]][indexlist[[i]]+k-2, 6]) + x_diff
            df.list[[a]][[6]][indexlist[[i]]-1+k, 7] <- as.numeric(df.list[[a]][[6]][indexlist[[i]]+k-2, 7]) + y_diff
            df.list[[a]][[6]][indexlist[[i]]-1+k, 8] <- "interpol"
            
          }
          
        }
        
      }
    }
    
  }  
  
  
  colnames(df.list[[a]][[6]]) <- c("x_og", "y_og", "likelihood", "mins", "phase", "x_new", "y_new", "origin")
  
  # constrain to arena borders
  
  for(i in 1:length(df.list[[a]][[6]][[1]])){
    if(df.list[[a]][[6]][i, 6] < 0){
      df.list[[a]][[6]][i, 6] <- 0
    }else if(df.list[[a]][[6]][i, 6] > arenax){
      df.list[[a]][[6]][i, 6] <- arenax
    }else if(df.list[[a]][[6]][i, 7] > arenay){
      df.list[[a]][[6]][i, 7] <- arenay
    }else if(df.list[[a]][[6]][i, 7] < 0){
      df.list[[a]][[6]][i, 7] <- 0
    }
  }
  
  print(paste("There are", sum(is.na(df.list[[a]][[6]][[6]])), "NAs left for the recording",
              df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]], df.list[[a]][[5]]$Dose[[1]], df.list[[a]][[5]]$other[[1]]))
  
  testdf <- df.list[[a]][[6]]
  for(i in 1:(length(df.list[[a]][[6]][[1]])-1)){
    # d = euclidian distance
    d <- sqrt((df.list[[a]][[6]][i+1, 6] - df.list[[a]][[6]][i, 6])^2 + (df.list[[a]][[6]][i+1, 7] - df.list[[a]][[6]][i, 7])^2)
    speed <- d/0.04008
    testdf[i, 9] <- d
    testdf[i, 10] <- speed
  }

  # ggplot(filter(testdf, phase != "injection"), aes(x = mins, y = V9)) +
  #   geom_line() +
  #   geom_hline(yintercept = mean(testdf$V9, na.rm = T) + sd(testdf$V9, na.rm = T)*2, col = "red") +
  #   labs(subtitle = paste("Recording:", df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]],
  #                         ifelse(!is.na(df.list[[a]][[5]]$Dose[[1]]),df.list[[a]][[5]]$Dose[[1]], "" ),
  #                         ifelse(!is.na(df.list[[a]][[5]]$other[[1]]),df.list[[a]][[5]]$other[[1]], "" ),
  #                         ifelse(k == 1, "baseline", "post injection")))
  
  testdf <- mutate(testdf, teleport = ifelse(testdf$V9 > mean(testdf$V9, na.rm = T) + sd(testdf$V9, na.rm = T)*2, "tele", "normal"))

  zone<- "norm"
  for(i in 1:(length(testdf[[1]])-1)){
    if(testdf[i, 5] != "injection"){
      if(testdf[i, 11] == "tele" & zone == "norm"){
        if(testdf[i+1, 6] == 0 | testdf[i+1, 6] == 50 | testdf[i+1, 7] == 0 | testdf[i+1, 7] == 50){
          zone <- "telezone"
          
        }
      }else if(zone == "telezone"){
        if( testdf[i+1, 11] == "tele"){
          zone <- "norm"
          testdf[i, 11] <- "tele"
        }else if(testdf[i+1, 6] == 0 | testdf[i+1, 6] == 50 | testdf[i+1, 7] == 0 | testdf[i+1, 7] == 50){
          zone <- "telezone"
          testdf[i, 11] <- "tele"
        }
      }
    }
    
  }
  
  
  counter <- 0
  lastx <- 0
  lasty <- 0
  lastxpos <- 0
  
  for(i in 1:(length(testdf[[1]])-1)){
    if(counter == 0 & testdf[i, 11] == "tele"){
      lastx <- testdf[i, 6]
      lasty <- testdf[i, 7]
      lastxpos <- i
      if(testdf[i+1, 11] == "tele"){
        counter <- 1
      }else if(testdf[i+1, 11] != "tele"){
        testdf[i+1, 12] <- (testdf[i+2, 6] - lastx)/2 + lastx
        testdf[i+1, 13] <- (testdf[i+2, 7] - lasty)/2 + lasty
        counter <- 0
      }
    }else if(counter != 0 & testdf[i, 11] == "tele"){
      if(testdf[i+1, 11] == "tele"){
        counter <- counter + 1
      }else{
        if(testdf[i+1, 6] != 0 | testdf[i+1, 6] != 50 | testdf[i+1, 7] != 0 | testdf[i+1, 7] != 50){
          counter <- counter + 1
          xdiff <- (testdf[i+1, 6] - lastx)/counter
          ydiff <- (testdf[i+1, 7] - lasty)/counter
          for(j in 1:(counter-1)){
            testdf[lastxpos+j, 12] <- lastx + xdiff*j
            testdf[lastxpos+j, 13] <- lasty + ydiff*j
          }
          counter <- 0
        }else{
          counter <- counter + 1
        }
        
      }
      
    }else if(counter != 0 & testdf[i, 11] != "tele"){
      if(testdf[i+1, 6] != 0 | testdf[i+1, 6] != 50 | testdf[i+1, 7] != 0 | testdf[i+1, 7] != 50){
        counter <- counter + 1
      }
    }
  }
  
  testdf <- mutate(testdf, V12 = ifelse(is.na(V12), x_new, V12), V13 = ifelse(is.na(V13), y_new, V13))
  
  
  df.list[[a]][[6]] <- testdf
  
}

#=================== Calculation of Outcome Variables ==========================


for(a in 1:length(df.list)){
  if(is.na(df.list[[a]][[5]]$other[[1]])){
    subdflist <- list()
    post <- filter(df.list[[a]][[6]], mins > inj_time_end)
    baseline <- filter(df.list[[a]][[6]], mins < inj_time_start)
    subdflist[[1]] <- baseline
    subdflist[[2]] <- post
    
    for(k in 1:length(subdflist)){
      speed_df <- as.data.frame(matrix(data = NA, nrow = length(subdflist[[k]][[1]])-movav-1, ncol = 6))
      
      for(i in 1: (length(subdflist[[k]][[1]])-1)){
        if(i <= length(subdflist[[k]][[1]]) - movav-1){ #any of the "middle" values
          speed_df[i, 1] <- subdflist[[k]][i, 12]
          speed_df[i, 2] <- subdflist[[k]][i+movav, 12]
          speed_df[i, 3] <- subdflist[[k]][i, 13]
          speed_df[i, 4] <- subdflist[[k]][i+movav, 13]
          speed_df[i, 5] <- subdflist[[k]][i, 4]*60
          speed_df[i, 6] <- subdflist[[k]][i+movav, 4]*60
        }else{ #end bit
        }
        
      }
  
      
      colnames(speed_df) <- c("xval", "nextx", "yval", "ynext", "sec", "nexttime")
      speed_df <- mutate(speed_df, dist = sqrt((nextx - xval)^2 + (ynext - yval)^2), timediff = nexttime - sec,
                         velocity = dist/timediff)
      

      if(k == 1){
        ind_seq <- seq(1, length(filter(speed_df, sec/60 >= baseline_start)[[1]]), movav)
      }else{
        ind_seq <- seq(1, length(filter(speed_df, sec/60 <= length_experiment)[[1]]), movav)
      }
      
      
      bin <- as.data.frame(matrix(data = NA, nrow = length(ind_seq), ncol = length(speed_df)+1))
      count <- 1
      for(i in ind_seq){
        bin[count, ] <- speed_df[i, ]
        bin[count, length(speed_df)+1] <- (count-1)*speed_df[i, 8]
        count <- count + 1
        
      }
      
      lb <- length(bin)
      
      for(i in 1:length(bin[[1]])){
        if(0 + centerx/2 < bin[i, 1] & arenax - centerx/2 > bin[i, 1] &
           0 + centery/2 < bin[i, 3] & arenay - centery/2 > bin[i, 3]){
          bin[i, lb+1] <- "center"
        }else{
          bin[i, lb+1] <- "edge"
        }
        
        if(i != 1){
          if(bin[i-1, lb+1] == "edge" & bin[i, lb+1] == "center"){
            bin[i, lb+2] <- 1
          }else{
            bin[i, lb+2] <- 0
          }
        }else{
          bin[i, lb+2] <- ifelse(bin[i, lb+1] == "center", 1, 0)
        }
        
      }
      
      
      
      colnames(bin) <- c("xval", "nextx", "yval", "ynext", "sec", "nexttime", "dist", "time_diff", "velocity", "reltime", "center","centerentry")
      
      
      twoDtrace <- ggplot(speed_df, aes(x = xval, y = yval))+
        annotate("rect", xmin = 0, xmax = arenax, ymin = 0, ymax = arenay, fill = "seashell")+
        annotate("rect", xmin = 0 + 0.5*centerx, xmax = arenax - 0.5*centerx, ymin = 0 + 0.5*centery, ymax = arenay - 0.5*centery, fill = "seashell3")+
        geom_path() +
        theme_minimal()+
        ggtitle("Trajectory in OF")+
        xlab("x") +
        xlim(0, arenax) +
        ylab("y")+
        ylim(0, arenay) +
        coord_fixed(ratio = 1) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_blank()) +
        labs(subtitle = paste("Recording:", df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]],
                              ifelse(!is.na(df.list[[a]][[5]]$Dose[[1]]),df.list[[a]][[5]]$Dose[[1]], "" ), 
                              ifelse(!is.na(df.list[[a]][[5]]$other[[1]]),df.list[[a]][[5]]$other[[1]], "" ),
                              ifelse(k == 1, "baseline", "post injection")))
      
      print("2D trace plot done")
      
      
      summary_stats <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 14))
      colnames(summary_stats) <- c("distance_moved", "med_speed", "perc_highspeed", "max_speed", 
                                   "centertime", "percent_center", "center_entries", "nrow", "df", "phase", "ID", "Compound", "Day",
                                   "Group")
      
      
      tornado_line <- plot_ly(bin, x = ~xval, y = ~yval, z = ~reltime/60, type = "scatter3d", mode = "lines",
                              opacity = 0.7,
                              line = list(width = 10, color = ~velocity, colorscale = list(c(0,'#4A7BB7'),
                                                                                           c(0.7, "#FEDA8B"), # 0 corresponds to min(z)
                                                                                           c(1, '#DD3D2D')), showscale = TRUE)) %>%
        layout(title = paste("Trajectory across time<br><sup>", df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]],
                             ifelse(!is.na(df.list[[a]][[5]]$Dose[[1]]),df.list[[a]][[5]]$Dose[[1]], "" ),
                             ifelse(!is.na(df.list[[a]][[5]]$other[[1]]),df.list[[a]][[5]]$other[[1]], "" ),
                             ifelse(k == 1, "baseline", "post injection"),"</sup>"),
               scene = list(xaxis = list(title = "X Coord"),
                            yaxis = list(title = "Y Coord"),
                            zaxis = list(title = "Time (mins)"),
                            aspectmode = "manual", aspectratio = list(x=3, y=3, z=20)),
               xaxis = list(range = c(0, 50)),
               yaxis = list(range = c(0, 50)))
      
      tornado_point <- plot_ly(bin, x = ~xval, y = ~yval, z = ~reltime/60,type = "scatter3d", mode = "marker",
                               opacity = 0.7,
                               marker = list(size = 4, color = ~velocity, colorscale = list(c(0,'#4A7BB7'),
                                                                                            c(0.7, "#FEDA8B"), # 0 corresponds to min(z)
                                                                                            c(1, '#DD3D2D')), showscale = TRUE)) %>%
        layout(title = paste("Trajectory across time<br><sup>", df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]],
                             ifelse(!is.na(df.list[[a]][[5]]$Dose[[1]]),df.list[[a]][[5]]$Dose[[1]], "" ),
                             ifelse(!is.na(df.list[[a]][[5]]$other[[1]]),df.list[[a]][[5]]$other[[1]], "" ),
                             ifelse(k == 1, "baseline", "post injection"),"</sup>"),
               scene = list(xaxis = list(title = "X Coord"),
                            yaxis = list(title = "Y Coord"),
                            zaxis = list(title = "Time (mins)"),
                            aspectmode = "manual", aspectratio = list(x=3, y=3, z=20)),
               xaxis = list(range = c(0, 50)),
               yaxis = list(range = c(0, 50)))
      
      
      print("plotly done")
      
      summary_stats[1, 1] <- sum(bin$dist)/100
      summary_stats[1, 2] <- median(bin$velocity)
      summary_stats[1, 3] <- length(filter(bin, velocity > 3))[[1]]/length(bin[[1]])
      summary_stats[1, 4] <- max(bin$velocity)
      summary_stats[1, 5] <- sum(filter(bin, center == "center")[[8]])/60
      summary_stats[1, 6] <- sum(filter(bin, center == "center")[[8]])/sum(bin[[8]])*100
      summary_stats[1, 7] <- sum(bin$centerentry)
      summary_stats[1, 8] <- dim(bin)[[1]]
      summary_stats[1, 9] <- "bin"
      summary_stats[1, 10] <- ifelse(k == 1, "baseline", "post")
      summary_stats[1, 11] <- df.list[[a]][[5]]$ID[[1]]
      summary_stats[1, 12] <- df.list[[a]][[5]]$Cpd[[1]]
      summary_stats[1, 13] <- df.list[[a]][[5]]$Day[[1]]
      summary_stats[1, 14] <- paste(df.list[[a]][[5]]$Group[[1]])
      
      if(k == 1){
        df.list[[a]][[7]] <- twoDtrace
        df.list[[a]][[8]] <- bin
        df.list[[a]][[9]] <- summary_stats
        df.list[[a]][[10]] <- tornado_point
        df.list[[a]][[11]] <- tornado_line

        # 
      }else{
        df.list[[a]][[12]] <- twoDtrace
        df.list[[a]][[13]] <- bin
        df.list[[a]][[14]] <- summary_stats
        df.list[[a]][[15]] <- tornado_point
        df.list[[a]][[16]] <- tornado_line

        
      }

      
      df.list[[a]][[17]] <- speed_df
    }
    
  }
  
  
  print(paste("Everything was saved for the recording",
              df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]], df.list[[a]][[5]]$Run[[1]], df.list[[a]][[5]]$other[[1]]))
  next
}

#============== Combining the data for all mice in a single df =================

combodf <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 14))
colnames(combodf) <- c("distance_moved", "med_speed", "%highspeed", "max_speed", 
                       "centertime", "centertime%", "center_entries", "nrow", "df", "phase", "ID", "Compound", "Run",
                       "Group")

for(a in 1:length(df.list)){
  for(b in 1:2){
    if(b == 1 & length(df.list[[a]]) > 9 ){
      if(is.na(combodf[1, 1]) & length(df.list[[a]][[9]]) != 0){
        combodf <- df.list[[a]][[9]]
      }else if(is.na(combodf[1, 1]) & length(df.list[[a]][[9]]) == 0){
        
      }else if(!is.na(combodf[1, 1]) & length(df.list[[a]][[9]]) != 0){
        combodf <- rbind(combodf, df.list[[a]][[9]])
      }
    }else{
      if(length(df.list[[a]]) > 14){
        if(is.na(combodf[1, 1]) & length(df.list[[a]][[14]]) != 0){
          combodf <- df.list[[a]][[14]]
        }else if(is.na(combodf[1, 1]) & length(df.list[[a]][[14]]) == 0){
          
        }else if(!is.na(combodf[1, 1]) & length(df.list[[a]][[14]]) != 0){
          combodf <- rbind(combodf, df.list[[a]][[14]])
        }
      }
    }
  }
}

#Fen22 is a second Fentanyl injection that was done on this particular mouse
combodf_use <- filter(combodf, Compound != "Fen22")

#Change the levels of Run for clarity when plotting and add extra column for the Nanobody that was injected
for(i in 1:length(combodf_use[[1]])){
  if(combodf_use[i, 14] == "Cn"){
    combodf_use[i, 14] <- "NbAlpha"
  }else{
    combodf_use[i, 14] <- "Nb64"
  }
  
  if(combodf_use[i, 13] == "R1"){
    combodf_use[i, 13] <- "Run1"
  }else{
    combodf_use[i, 13] <- "Run2"
  }
}


colnames(combodf_use)[14] <- c("Nanobody")


#Calculate the mean fentanyl response between the two injections
fen <- filter(combodf_use, Compound != "NaCl", df == "bin", phase == "post")
mean_fen <- fen %>%
  group_by(ID, Nanobody, Day) %>%
  summarise(distance = mean(distance_moved), speed = mean(med_speed))

#Extend the combined data by the mean fentanyl response
len <- length(combodf_use[[1]])
for(i in 1:length(mean_fen[[1]])){
  combodf_use[len+i, 1] <- mean_fen[i, 4]
  combodf_use[len+i, 2] <- mean_fen[i, 5]
  combodf_use[len+i, 9] <- "bin"
  combodf_use[len+i, 10] <- "post"
  combodf_use[len+i, 11] <- mean_fen[i, 1]
  combodf_use[len+i, 13] <- mean_fen[i, 3]
  combodf_use[len+i, 14] <- mean_fen[i, 2]
}


combodf_use <- mutate(combodf_use, Cpd = ifelse(is.na(Compound), "mean_Fen", Compound))
combodf_use <- mutate(combodf_use, Compound = as.character(Compound))
for(i in 1:length(combodf_use[[1]])){
  if(is.na(combodf_use[i, 12])){
    combodf_use[i, 12] <- combodf_use[i, 15]
  }
}
combodf_use <- mutate(combodf_use, Compound = as.factor(Compound))

#Order the factor levels, this is useful to set the order of the factors on plots
combodf_use$Nanobody <- factor(combodf_use$Nanobody, levels = c("NbAlpha", "Nb64")) 
combodf_use$Compound <- factor(combodf_use$Compound, levels = c("NaCl", "Fen1", "Fen2", "mean_Fen")) 

colnames(combodf_use)[13] <- c("Run")

combodf_use <- select(combodf_use, !Cpd)

#This df can now be used to plot, either in R or also using softwares like GraphPad Prism


#=============================== Speed DF ======================================

for(i in 1:length(df.list)){
  mouseinfo_long <- cbind(rep(df.list[[i]][[5]][1, 35], length(df.list[[i]][[13]][[1]])), rep(df.list[[i]][[5]][1, 36], length(df.list[[i]][[13]][[1]])), rep(df.list[[i]][[5]][1, 37], length(df.list[[1]][[13]][[1]])),
                          rep(df.list[[i]][[5]][1, 38], length(df.list[[i]][[13]][[1]])))
  mouseinfo_short <- cbind(rep(df.list[[i]][[5]][1, 35], length(df.list[[i]][[8]][[1]])), rep(df.list[[i]][[5]][1, 36], length(df.list[[i]][[8]][[1]])), rep(df.list[[i]][[5]][1, 37], length(df.list[[1]][[8]][[1]])),
                           rep(df.list[[i]][[5]][1, 38], length(df.list[[i]][[8]][[1]])))
  if(i == 1){
    speed <- cbind(df.list[[i]][[13]]$reltime, df.list[[i]][[13]]$velocity, mouseinfo_long, rep("Post", length(df.list[[i]][[13]][[1]])))
    colnames(speed)[1:7] <- c("Time", "Velocity", "ID", "Group", "Compound", "Run", "Phase")
    
    
    speed_base <- cbind(df.list[[i]][[8]]$reltime, df.list[[i]][[8]]$velocity, mouseinfo_short, rep("Baseline", length(df.list[[i]][[8]][[1]])))
    colnames(speed_base)[1:7] <- c("Time", "Velocity", "ID", "Group", "Compound", "Run", "Phase")
    
  }else{
    
    midstep1 <- cbind(df.list[[i]][[13]]$reltime, df.list[[i]][[13]]$velocity, mouseinfo_long, rep("Post", length(df.list[[i]][[13]][[1]])))
    colnames(midstep1)[1:7] <- c("Time", "Velocity", "ID", "Group", "Compound", "Run", "Phase")
    
    speed <- rbind(speed, midstep1)
    
    midstep2 <- cbind(df.list[[i]][[8]]$reltime, df.list[[i]][[8]]$velocity, mouseinfo_short, rep("Baseline", length(df.list[[i]][[8]][[1]])))
    colnames(midstep2)[1:7] <- c("Time", "Velocity", "ID", "Group", "Compound", "Run", "Phase")
    
    speed_base <- rbind(speed_base, midstep2)
    
  }
  
}

speed <- rbind(speed, speed_base)
speed <- as.data.frame(speed)
speed <- mutate(speed, Time = as.numeric(Time)/60, Velocity = as.numeric(Velocity), ID = as.character(ID), Group = as.character(Group), Compound = as.character(Compound),
                Run = as.character(Run), Phase = as.character(Phase))

#for plotting purposes
baseline_length <- 5
gap <- 0.5

speed <- speed %>%
  mutate(Mins = ifelse(Phase == "Baseline", Time, Time + baseline_length + gap))

speed <- speed %>%
  mutate(Group = ifelse(Group == "Cn", "NbAlpha", "Nb64"))

speed$Compound <- factor(speed$Compound, levels = c("NaCl", "Fen1", "Fen2")) 
speed$Nanobody <- factor(speed$Group, levels = c("NbAlpha", "Nb64")) 

speed <- speed %>%
  arrange(Group, ID) %>%
  mutate(ID = factor(ID, levels = unique(ID)))

# Heatmap code
ggplot(speed, aes(x = Mins, y = ID, fill = Velocity))+
  geom_tile() +
  facet_grid(Nanobody ~ Compound, scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_viridis_c(option = "magma", name = "Velocity (cm/s)", trans = "log1p")+
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_blank()) +
  xlab("Time (mins)")

#================ Tornado plots with comparable colour scheme ==================

# In the above code, the speed colour code is not comparable across recordings.
# If a common colour code is wanted, one has to extract the bin data from each mouse 
# and create a single data frame with speed measurements for all timepoints (= speed in line 989).
# Find code for this in the previous section
# With this kind of data frame, one can also show the speed as heat maps (like for photometry measurements) etc.


for(a in 1:length(df.list)){
  for(i in 1:length(df.list[[a]])){
    if(i == 8 | i == 13){
      tornado_line <- plot_ly(df.list[[a]][[i]], x = ~xval, y = ~yval, z = ~reltime/60, type = "scatter3d", mode = "lines",
                              opacity = 0.7, 
                              line = list(width = 10, color = ~velocity, colorscale = list(c(0,'#4A7BB7'),
                                                                                           c(0.25, "#FEDA8B"), # 0 corresponds to min(z)
                                                                                           c(0.75, '#DD3D2D'),
                                                                                           c(1,'#434343')), cmin = 0, cmax = quantile(speed$Velocity, 0.98, na.rm = T), showscale = TRUE)) %>%
        layout(title = paste("Trajectory across time<br><sup>", df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]],
                             ifelse(!is.na(df.list[[a]][[5]]$Dose[[1]]),df.list[[a]][[5]]$Dose[[1]], "" ), 
                             ifelse(!is.na(df.list[[a]][[5]]$other[[1]]),df.list[[a]][[5]]$other[[1]], "" ),
                             ifelse(i == 8, "baseline", "post injection"),"</sup>"),
               scene = list(xaxis = list(title = "X Coord"),
                            yaxis = list(title = "Y Coord"),
                            zaxis = list(title = "Time (mins)"),
                            aspectmode = "manual", aspectratio = list(x=3, y=3, z=20)),
               xaxis = list(range = c(0, 50)),
               yaxis = list(range = c(0, 50)))
      

      if(i == 8){
        fig <- tornado_line %>%
          plotly_build
        name <- paste(df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]], df.list[[a]][[5]]$Dose[[1]], df.list[[a]][[5]]$other[[1]],
                      "binned_baseline", sep = "_")
        saveWidget(fig, file = paste(name, "html", sep = "."))
      }else{
        fig <- tornado_line %>%
          plotly_build
        name <- paste(df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]], df.list[[a]][[5]]$Dose[[1]], df.list[[a]][[5]]$other[[1]],
                      "binned_post", sep = "_")
        saveWidget(fig, file = paste(name, "html", sep = "."))
      }
    }
  }
  
}


#============================= Saving Figures ==================================

for(a in 1:length(df.list)){
  for(i in 1:length(df.list[[a]])){
    if(i == 7 & length(df.list[[a]][[i]]) != 0){
      
      ggsave(paste(df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]], df.list[[a]][[5]]$Dose[[1]],
                   "twoDtrace_baseline.jpeg", sep = "_"), plot = df.list[[a]][[i]])
      
    }else if(i == 12 & length(df.list[[a]][[i]]) != 0){
       
      ggsave(paste(df.list[[a]][[5]]$ID[[1]], df.list[[a]][[5]]$Cpd[[1]], df.list[[a]][[5]]$Dose[[1]],
                   "twoDtrace_post.jpeg", sep = "_"), plot = df.list[[a]][[i]])
      
    }
    
  }
}



