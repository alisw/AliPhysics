// this is a macro to setup the OCDB whose objects will be used
// to undo the alignment/calibration applied during reconstruction

void configRecoOCDB(int run = 188503) 
{
  //
  AliAlgAux::CleanOCDB();
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->SetRaw(1);
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    man->SetSnapshotMode("OCDB.root");
  }
  //
  man->SetRun(run);

}
