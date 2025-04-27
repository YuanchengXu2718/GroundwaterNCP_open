# GroundwaterNCP_open
Here you can find the data, code, and source used in our paper about rapid groundwater recovery in the North China Plain (NCP)
  
## GroundwaterDepth.xlsx
In-situ groundwater depth data in the NCP from 2005 to 2024
### Type of monitoring well
If 10 < type <= 19, the well monitors unconfined aquifers.  
Type 11 represents the first unconfined aquifers.  
Type 12 represents the second unconfined aquifers, etc.  
  
If 20 < type <= 29, the well monitors confined aquifers.  
Type 21 represents the first confined aquifers.  
Type 22 represents the second confined aquifers, etc.


## GW_Dph_NN.m
Spatial interpolation based on Thiessen polygons (Nearest neighbour)  
  
## dphAnalysis.m
Code for spatial interpolation and uncertainty analysis  
GoundwaterDepth.xlsx and GW_Dph_NN.m are used in this code 
  
Check the "CoreNum = 4" in LINE 19 and ensure the number of cores used is smaller than you have!  
  
Some users have reported encountering unexpected errors when using "parfor" in MATLAB, which may be related to MATLAB's registration. If you experience the same issue and are unable to resolve it, please consider replacing "parfor" with "for".  
  
Tips: Monte Carlo simulations for uncertainty analysis are time-consuming. For example, running 1,000 simulations for average groundwater depth in the NCP takes ~20 min with Intel i9-14900K and parallel computing. For older laptops with Intel i7-7700HQ, the same process can last over 3 hours.


## dphAnalysis_preprocessingPara.m
Code for examining the uncertainty introduced by preprocessing  
GoundwaterDepth.xlsx and GW_Dph_NN.m are used in this code   
  
## Folder: mask
Find the mask for different regions:
### North China Plain
ncp.tif
### City: Beijing, Baoding, Hengshui, and Tangshan
city_beijing.tif, city_baoding.tif, city_hengshui.tif, and city_tangshan.tif
  
  
## sourceDataFigure.xlsx
Source data for the figures
  
  
## WaterSupplyUse_HRB_2005_2023.xlsx
Annual water supply and usage in the Hai River basin (HRB), including precipitation and renewable water resources
