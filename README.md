# Nb64
This repository contains the code written to analyse the fiber photometry and behavioural data in the manuscript: "An engineered nanobody inhibitor for molecular-to-circuit control of opioid receptor function" (BioRxiv). doi: https://doi.org/10.64898/2026.01.12.698943


### Notes on the File "Behavioural_analysis_Manuscript.R"

#### This was written on RStudio (ver. 1.4.1106, R ver. 4.0.4)

##### Data structure
This is a code to analyse Open Field video data processed in DeepLabCut [(DLC)](https://deeplabcut.github.io/DeepLabCut/README.html).  
The data processed in DLC has the following structure:  
Column 1: Frame number  
Columns 2-4: Predicted x coordinate of a keypoint k, predicted y coordinate of k, likelihood of the prediction   
Columns 5-7: Predicted x coordinate of a keypoint w, predicted y coordinate of w, likelihood of the prediction  
... for all detected keypoints.  
It is necessary to have at least one keypoint tracking the mouse (the best suited bodypart will depend on the behavioural arena and camera angle; typically it will be the center of the mouse body).  
Suitability depends on factors such as occlusion, and tracking confidence (likelihood of the prediction given by DLC). In addition, all corners of the behavioural arena have to be tracked as keypoints.  

The naming convention of the output files is: date(yyyymmdd)_mousecode(X_Ex/Cn_1 etc)_otherexperimentalinfo

##### Libraries on R
The code assumes that the necessary packages are already installed (if not, use install.packages("package_of_interest"). If the code is only used to attain the outcome variables, only the tidyverse package is needed. If plotting is desired, the packages plotly, lattice, htmlwidgets, cowplot, plotrix and viridis may be used.  

##### Output
- Main readouts (in numbers) are (separately for baseline and post injection): Distance moved, Median velocity, percent of time the mouse is running at high speed (>3 cm/s), maximum speed, time spent in center, percent of time spent in center, number of center entries,
- A file with the velocity across time for each recording
- A 2D trajectory plot for each baseline and post injection
- A line tornado plot with comparable speed colour coding across recordings and mice for each baseline and post injection (there is also code for a tornado plot with each timepoint represented as a dot (in the section "Calculation of Outcome Variables), so this could easily be used and saved also)
- There is also a code chunk plotting speed across the recording as heatmaps

##### Outline by Code Section
More details can be found as comments in the code.  
- First, some variables concerning experimental information are defined (adjust to your needs).
- All files (.csv) in the working directory are read in.
- Data is collected in a large data frame (df), where each position is a sublist for one recording, where the data frames and graphs created during analysis will be stored.
- **Housekeeping:** Renaming of columns in of the raw data, adding mouse information to the df and creating a time variable based on the frame number and frame rate.
- **Likelihood histograms:** The likelihood of each keypoint is plotted across time (could be used to judge which keypoint should be used to follow the mouse, or to see how stable the detection of certain keypoints is).
- **Checking the corners:** Extraction of the corner points (median coordinates of all frames) of the arena and plotting them to check for detection problems (cannot be expected to be completely squared due to image distortion by the camera).
- **Homography matrix calculation:** Corner coordinates and real arena bounds are used to estimate distortion of the image.
- **Application of Homography to other keypoints:** Using the homography matrix calculated in the previous step and apply it to all other detected keypoints to account for image distortion when calculating distance moved, speed, etc.
- **Interpolation of low likelihood points:** Based on the likelihood given by DLC, low likelihood detections are flagged and interpolated, assuming the mouse moves at a constant speed between the last and next high likelihood detection. In case there are either high likelihood detections outside of arena bounds and/or teleportation events (euclidian distance > mean distance moved + standard deviation of distance moved*2) are flagged and interpolated in the same day.
- **Calculation of Outcome Variables:** Speed is binned in steps of 4 (the movav variable defined at the beginning of the code), the 2D trace is plotted, as well as the tornado plots (but with non-comparable colour scales). Variables such as center entries are created based on mouse location and the defined center of the arena (at the beginning of the code). Summary statistics (outcome variables listed in the first point under the Output section) are collected and saved in the large data frame.
- **Combining the data for all mice in a single df:** Collecting the outcome variables from all runs in a single df that can be saved and used to plot in any software of preference. In addition, since there were two fentanyl injections, there is also a code chunk calculating the mean for distance moved and median speed for the two injection days.
- **Speed DF:** This is a data frame where the speed data across time is extracted for each recording. This can be used to plot speed across time in multiple ways and is necessary if a comparable speed colour code is necessary for the tornado plots (next section). Contains a code chunk for a heatmap of speed, which can be saved.
- **Tornado plots with comparable colour scheme:** The tornado plots as of now do not have comparable colour coding of speeds across recordings, which is necessary when the data should be visually comparable. For this, it is necessary to have all speed data from all recordings (previous section). The tornado plots can be plotted as points instead of lines (code in the "Calculation of Outcome Variables" section). The plots are automatically saved in the working directory. 
- **Saving Figures:** Currently only saves the 2D trajectory plots, also other plots and dfs could be saved with small adjustments. 



### Notes on File "nLightG2_HPC_Photometry_v20250213pl.m"

#### This was written on Matlab for analysis of the Fiber Photometry data
