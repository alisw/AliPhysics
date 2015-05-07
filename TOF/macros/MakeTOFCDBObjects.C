#include "TOF/MakeTOFZeroMisAlignment.C"
#include "TOF/CreateConfigMap.C"
#include "TOF/CreateConfigMapNoise.C"
#include "TOF/MakeCDBEntryCTPLatency.C"
#include "TOF/MakeCDBEntryDeltaBCOffset.C"
#include "TOF/CreateIdealOnlineCalibPars.C"
#include "TOF/CreateCalibPars_Ideal.C"
#include "TOF/MakeCDBEntryProblematic.C"
#include "TOF/MakeCDBEntryReadoutEfficiency.C"
#include "TOF/MakeTOFRecoParam.C"
#include "TOF/MakeCDBEntryRunParams.C"
#include "TOF/MakeCDBEntryT0Fill.C"
#include "TOF/MakeCDBEntryT0FillOnlineCalib.C"

void MakeTOFCDBObjects()
{
  MakeTOFZeroMisAlignment();         // TOF/Align/Data                            
  CreateConfigMap();                 // TOF/Calib/Config                          
  CreateConfigMapNoise();            // TOF/Calib/ConfigNoise                     
  MakeCDBEntryCTPLatency();          // TOF/Calib/CTPLatency                      
  MakeCDBEntryDeltaBCOffset();       // TOF/Calib/DeltaBCOffset                   
  CreateIdealOnlineCalibPars();      // TOF/Calib/HW                              
  CreateIdealOnlineCalibPars();      // TOF/Calib/Noise                           
  CreateCalibPars_Ideal();           // TOF/Calib/ParOffline                      
  CreateIdealOnlineCalibPars();      // TOF/Calib/ParOnline                       
  CreateIdealOnlineCalibPars();      // TOF/Calib/ParOnlineDelay                  
  MakeCDBEntryProblematic();         // TOF/Calib/Problematic                     
  CreateIdealOnlineCalibPars();      // TOF/Calib/Pulser                          
  MakeCDBEntryReadoutEfficiency();   // TOF/Calib/ReadoutEfficiency               
  MakeTOFRecoParam();                // TOF/Calib/RecoParam                       
  MakeCDBEntryRunParams(0.,0.);      // TOF/Calib/RunParams                       
  CreateCalibPars_Ideal();           // TOF/Calib/SimHisto                        
  CreateIdealOnlineCalibPars();      // TOF/Calib/Status                          
  MakeCDBEntryT0Fill();              // TOF/Calib/T0Fill                          
  MakeCDBEntryT0FillOnlineCalib();   // TOF/Calib/T0FillOnlineCalib               
}
