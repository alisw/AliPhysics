// this is a macro to setup the OCDB whose objects will be used
// as a reference for the alignment/calibration, i.e. coorections
// will be evaluated wrt these objects

TString specOCDB = "";//"/alice/cern.ch/user/s/shahoian/tstOCDB/outOCDB_LHC12tst0";
//TString specOCDB = "outOCDB";

void configRefOCDB(int run = 225105) 
{
  //
  int ver = -1;
  int sver = -1;

  AliAlgAux::CleanOCDB();
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->SetRaw(1);
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    man->SetSnapshotMode("OCDB.root");
  }
  //
  if (!specOCDB.IsNull()) {
    if (specOCDB.BeginsWith("/alice/cern.ch")) {
      man->SetSpecificStorage("ITS/Align/Data",Form("alien://folder=%s",specOCDB.Data()),ver<0 ? -1:ver, sver<0 ? -1:sver);
      man->SetSpecificStorage("TRD/Align/Data",Form("alien://folder=%s",specOCDB.Data()),ver<0 ? -1:ver, sver<0 ? -1:sver);
      man->SetSpecificStorage("TOF/Align/Data",Form("alien://folder=%s",specOCDB.Data()),ver<0 ? -1:ver, sver<0 ? -1:sver);
    }
    else {
      man->SetSpecificStorage("ITS/Align/Data",Form("local://%s",specOCDB.Data()),ver<0 ? -1:ver, sver<0 ? -1:sver);
      man->SetSpecificStorage("TRD/Align/Data",Form("local://%s",specOCDB.Data()),ver<0 ? -1:ver, sver<0 ? -1:sver);
      man->SetSpecificStorage("TOF/Align/Data",Form("local://%s",specOCDB.Data()),ver<0 ? -1:ver, sver<0 ? -1:sver);
    }
  }
  //
  man->SetRun(run);
  //
}
