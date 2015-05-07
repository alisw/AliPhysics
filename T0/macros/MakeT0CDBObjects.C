#include "T0/MakeT0ZeroMisAlignment.C"
#include "T0/setlookUpTable.C"
#include "T0/MakeT0RecoParam.C"

void MakeT0CDBObjects()
{
  MakeT0ZeroMisAlignment(); // T0/Align/Data                             
  setlookUpTable();         // T0/Calib/LookUp_Table                     
  MakeT0RecoParam();        // T0/Calib/RecoParam                        
}

/*  Macro not yet defined for:
T0/Calib/Latency
T0/Calib/Slewing_Walk
T0/Calib/TimeAdjust
T0/Calib/TimeDelay
T0/Calib/TriggerParam
*/
