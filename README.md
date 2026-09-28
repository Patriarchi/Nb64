# Nb64
This repository contains the code written to analyse the fiber photometry and behavioural data in the manuscript: "An engineered nanobody inhibitor for molecular-to-circuit control of opioid receptor function" (BioRxiv). doi: https://doi.org/10.64898/2026.01.12.698943

====================================================================
Notes on the File "Behavioural_analysis_Manuscript.R"
        *        *        *        *        *        *
This is a code to analyse Open Field video data processed in DeepLabCut (https://deeplabcut.github.io/DeepLabCut/README.html, DLC).
The data processed in DLC has the following structure:
Column 1: Frame number
Columns 2-4: Predicted x coordinate of a keypoint k, predicted y coordinate of k, likelihood of the prediction 
Columns 5-7: Predicted x coordinate of a keypoint w, predicted y coordinate of w, likelihood of the prediction
... for all detected keypoints.
It is necessary to have at least one keypoint tracking the mouse (the best suited bodypart will depend on the behavioural arena and camera angle; typically it will be the center of the mouse body).
Suitability depends on factors such as occlusion, and tracking confidence (likelihood of the prediction given by DLC). In addition, all corners of the behavioural arena have to be tracked as keypoints.

The naming convention of the output files is: date(yyyymmdd)_mousecode(X_Ex/Cn_1 etc)_otherexperimentalinfo
