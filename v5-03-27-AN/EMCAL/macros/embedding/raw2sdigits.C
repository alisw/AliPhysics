
// Transform raw data into digits
// Input: run number, raw data file and ESD file from this data

// After executing this, generate simulation and merge with it using
// the macro simmerge.C

// Full exercise can be done automatically with other alice systems
// in the examples in $ALICE_ROOT/test/embedding

// Author: GCB (from example in $ALICE_ROOT/test/embedding)

void raw2sdigits(Int_t runNumber) 
{
  AliSimulation sim;
  AliCDBManager *cdbm = AliCDBManager::Instance();
  cdbm->SetRun(runNumber);
  cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");     
  cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  sim.SetMakeSDigits("EMCAL");  
  sim.ConvertRaw2SDigits("raw/raw.root","raw/AliESDs.root");
}
