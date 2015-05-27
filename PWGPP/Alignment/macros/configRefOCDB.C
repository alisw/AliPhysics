// this is a macro to setup the OCDB whose objects will be used
// as a reference for the alignment/calibration, i.e. coorections
// will be evaluated wrt these objects

void configRefOCDB(int run = 188503) 
{
  //
  AliAlgAux::CleanOCDB();
  AliCDBManager* man = AliCDBManager::Instance();
  //
  if (gSystem->AccessPathName("data/OCDB.root", kFileExists)==0) {        
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");    
    man->SetRun(run);
    man->SetSnapshotMode("data/OCDB.root");
  }
  else {
    man->SetRaw(1);
    man->SetRun(run);    
  }
  //
}
