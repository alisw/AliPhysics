#include "VZERO/MakeVZEROZeroMisAlignment.C" 
//#include "VZERO/MakeVZEROCalibEntry.C" gets the entry from AliEn
#include "VZERO/MakeVZEROEqualizationFactorsEntry.C"
#include "VZERO/MakeVZEROLightYieldsEntry.C"
#include "VZERO/MakeVZEROPMGainsEntry.C"
#include "VZERO/MakeVZERORecoParam.C"
#include "VZERO/MakeVZEROSaturationEntry.C"
#include "VZERO/MakeVZEROTimeDelaysEntry.C"
#include "VZERO/MakeVZEROTimeSlewingEntry.C"
// #include "VZERO/MakeVZEROTriggerEntry.C" gets the entry from AliEn

void MakeVZEROCDBObjects()
{
  MakeVZEROZeroMisAlignment();         // VZERO/Align/Data                   
  //MakeVZEROCalibEntry(0);               // VZERO/Calib/Data       HVfix  PbPb  PbPbNew  Scaled
  MakeVZEROEqualizationFactorsEntry(); // VZERO/Calib/EqualizationFactors    
  // MakeVZEROEqualizationFactorsFor2010pp.C // VZERO/Calib/EqualizationFactors    
  MakeVZEROLightYieldsEntry();         // VZERO/Calib/LightYields            
  // MakeVZEROLightYieldsEntryAging.C        // VZERO/Calib/LightYields            
  MakeVZEROPMGainsEntry();             // VZERO/Calib/PMGains                
  MakeVZERORecoParam();                // VZERO/Calib/RecoParam              
  MakeVZEROSaturationEntry();          // VZERO/Calib/Saturation             
  MakeVZEROTimeDelaysEntry();          // VZERO/Calib/TimeDelays             
  MakeVZEROTimeSlewingEntry();         // VZERO/Calib/TimeSlewing            
  //MakeVZEROTriggerEntry();             // VZERO/Trigger/Data                 
}

/*  Macro not yet defined for:
VZERO/PbPb/VZERO/Calib/Data
*/
