void MakeLHCClockPhaseEntry(const char *cdbStorage = "local://$ALICE_ROOT/OCDB")
{
  // Example macro to put in OCDB the default (=0) LHC-clock phase
  // It is valid fro runs from 0 to inf
  // The timestamp range is also inf (we store the first and last value for
  // each beam)
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbStorage);

  AliLHCClockPhase phaseObj;

  phaseObj.AddPhaseB1DP(0,0.);
  phaseObj.AddPhaseB2DP(0,0.);

  phaseObj.AddPhaseB1DP(2147483647,0.);
  phaseObj.AddPhaseB2DP(2147483647,0.);

  AliCDBMetaData* metadata = new AliCDBMetaData();
  metadata->SetResponsible("Cvetan Cheshkov");
  metadata->SetComment("Default LHC-clock phase object");
  AliCDBId id("GRP/Calib/LHCClockPhase",0,AliCDBRunRange::Infinity());

  man->Put(&phaseObj,id,metadata);

  return;
}
