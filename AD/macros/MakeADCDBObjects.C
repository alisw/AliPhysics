#include "AD/macros/MakeADZeroMisAlignment.C"
#include "AD/macros/MakeADRecoParamEntry.C"
#include "AD/macros/MakeADTimeDelaysEntry.C"
#include "AD/macros/MakeADLightYieldsEntry.C"
#include "AD/macros/MakeADPMGainsEntry.C"
#include "AD/macros/MakeADCalibEntry.C"
#include "AD/macros/MakeADTimeSlewingEntry.C"
#include "AD/macros/MakeADQAParamEntry.C"

void MakeADCDBObjects()
{
  MakeADZeroMisAlignment();   // AD/Align/Data                   
  MakeADRecoParamEntry();     // AD/Calib/RecoParam              
  MakeADTimeDelaysEntry();    // AD/Calib/TimeDelays             
  MakeADLightYieldsEntry();   // AD/Calib/LightYields            
  MakeADPMGainsEntry();       // AD/Calib/PMGains                
  MakeADCalibEntry();         // AD/Calib/Data                   
  MakeADTimeSlewingEntry();   // AD/Calib/TimeSlewing            
  MakeADQAParamEntry();       // AD/Calib/QAParam
}
