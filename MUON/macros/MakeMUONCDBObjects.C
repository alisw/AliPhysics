#include "AliMpCDB.h"
#include "AliMUONCDB.h"
#include "MUON/macros/MakeMUONRecoParamArray.C"
void MakeMUONCDBObjects()
{
  AliMpCDB::LoadAll2();
  AliMUONCDB::WriteTrigger();  // MUON/Calib/LocalTriggerBoardMasks MUON/Calib/TriggerDCS MUON/Calib/TriggerEfficiency MUON/Calib/TriggerLut
  AliMUONCDB::WriteTracker();  // MUON/Calib/MappingData MUON/Calib/MappingRunData MUON/Calib/HV MUON/Calib/Pedestals MUON/Calib/OccupancyMap MUON/Calib/RejectList MUON/Calib/Config MUON/Calib/LV
  AliMUONCDB::WriteRegionalTriggerConfig();  //MUON/Calib/RegionalTriggerConfig
  AliMUONCDB::WriteGlobalTriggerConfig();    //MUON/Calib/GlobalTriggerCrateConfig
  MakeMUONRecoParamArray();                  //MUON/Calib/RecoParam
}

/* to be done
MUON/Align/Baseline
MUON/Align/Data
MUON/Calib/GlobalTriggerBoardMasks
MUON/Calib/Neighbours                                  on the TODO list
MUON/Calib/RegionalTriggerBoardMasks
*/
