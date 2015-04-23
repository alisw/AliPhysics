#include "PHOS/MakePHOSZeroMisAlignment.C"
#include "PHOS/MakePHOSEmcCalib.C"
#include "PHOS/MakePHOSEmcCalib.C"
#include "PHOS/MakePHOSEmcCalib.C"
#include "PHOS/MakePHOSAltroMapping.C"
#include "PHOS/MakePHOSRecoParam.C"
#include "PHOS/MakePHOSTrigger.C"

void MakePHOSCDBObjects()
{
   MakePHOSZeroMisAlignment();  // PHOS/Align/Data
   MakePHOSEmcCalib();          // PHOS/Calib/CpvGainPedestals PHOS/Calib/EmcBadChannels PHOS/Calib/EmcGainPedestals
   MakePHOSAltroMapping();      // PHOS/Calib/Mapping
   MakePHOSRecoParam();         // PHOS/Calib/RecoParam
   MakePHOSTrigger();           // PHOS/Trigger/Parameters
}
