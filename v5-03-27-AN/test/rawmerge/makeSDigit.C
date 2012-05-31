void makeSDigit(Int_t run=117048,
                const char *fname="/alice/data/2010/LHC10b/000117048/raw/10000117048041.30.root")
{
  AliSimulation sim;
  //AliLog::SetGlobalDebugLevel(1);
  AliLog::SetModuleDebugLevel("STEER",1);
  AliCDBManager *cdbm = AliCDBManager::Instance();
  cdbm->SetRun(run);
  cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");     
  //cdbm->SetDefaultStorage("raw://");     
  cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s","$ALICE_ROOT/test/rawmerge"));  
  sim.SetMakeSDigits("ITS TPC");
  sim.SetNumberOfEvents(10000);
  sim.ConvertRaw2SDigits(fname);
  return;
}
