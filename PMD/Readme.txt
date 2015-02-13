PMD Codes
++++++++++++++++++++++++++
Author: Satyajit Jena
Tue Dec 16 17:12:47 CET 2014 
+++++++++++++++++++++++++

CMakeLists.txt : List of modules

PMDbase: PMD base classes        
PMDsim : Simulation Module        
PMDrec : Reconstructoin Module         
DA     : Detector Algorithm Module            

anal
    : Analysis related codes, calibration, cleanup
    : template reader
           
info/
    : To keep brief infos

data/ 
    : Ascii mapping info + some default data           
macro/ 
    : Contains all Example Files

ocdbmacro:/
It contains all macros that creates the default CDB 
objects we have in the repository's OCDB


MakePMDHotCDB.C             : OCDB/PMD/Calib/Hot     
MakePMDPedCDB.C             : OCDB/PMD/Calib/Ped
MakePMDGainCDB.C            : OCDB/PMD/Calib/Gain
MakePMDRecoParam.C          : OCDB/PMD/Calib/RecoParam 
MakePMDMappingCDB.C         : OCDB/PMD/Calib/Mapping
MakePMDDDLinfoCDB.C         : OCDB/PMD/Calib/Ddlinfo
MakePMDNoiseCutCDB.C        : OCDB/PMD/Calib/NoiseCut  
MakePMDResMisAlignment.C    : OCDB/PMD/Align/Data
MakePMDZeroMisAlignment.C   : OCDB/PMD/Align/Data
MakePMDFullMisAlignment.C   : OCDB/PMD/Align/Data



Readme.txt : This file    

