#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSTriggerConditions.h"
#include "AliCDBRunRange.h"
#endif
void defaultObjSPD(Int_t runNrStart=0, Int_t runNrEnd=AliCDBRunRange::Infinity(), const Char_t* storage = "local:///tmp/defaultSPD/")
{
  AliITSOnlineCalibrationSPDhandler *h = new AliITSOnlineCalibrationSPDhandler();
  h->WriteDeadToDB(runNrStart,runNrEnd,storage);
  h->WriteNoisyToDB(runNrStart,runNrEnd,storage);
  h->WriteSparseDeadToDB(runNrStart,runNrEnd,storage);
  //h->ReadPITConditionsFromDB(0,"local://$ALICE_ROOT/OCDB/");
  h->ReadPITConditionsFromDB(0,"alien://folder=/alice/data/2015/OCDB");
  AliITSTriggerConditions* tr = h->GetTriggerConditions();
  tr->ResetAll();
  h->WritePITConditionsToDB(runNrStart,runNrEnd,storage);
}
