/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * OCDB configuration script
 *
 */

enum EOCDB_t {
  kOCDBDefault,
  kOCDBCustom,
  kNOCDBs
};

const Char_t *OCDBName[kNOCDBs] = {
  "Default",
  "Custom"
};

/*****************************************************************/

OCDBConfig(Int_t tag, Int_t run, Int_t type)
{

  switch (tag) {
    
    // kOCDBDefault
  case kOCDBDefault:
    OCDBDefault(run, type);
    break;
    
    // kOCDBCustom
  case kOCDBCustom:
    if ((gROOT->LoadMacro("OCDBCustom.C")) != 0) {
      printf("ERROR: cannot find OCDBCustom.C\n");
      abort();
      return;
    }
    OCDBCustom(run);
    break;
    
  }

}

/*****************************************************************/

OCDBDefault(Int_t run, Int_t mode)
{

  gROOT->LoadMacro("Sim/Utils.C");
  Int_t year = RunToYear(run);
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(Form("alien://Folder=/alice/data/%d/OCDB", year));
  man->SetRun(run);
  
  //
  // set detector specific paths
  DefaultSpecificStorage(man, run, mode);

}

/*****************************************************************/

DefaultSpecificStorage(AliCDBManager *man, Int_t run, Int_t mode)
{
  
  gROOT->LoadMacro("Sim/Utils.C");
  Int_t year = RunToYear(run);

  const Char_t *Raw      = Form("alien://Folder=/alice/data/%d/OCDB", year);
  const Char_t *Ideal    = "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/";
  const Char_t *Residual = "alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/";
  const Char_t *Full     = "alien://Folder=/alice/simulation/2008/v4-15-Release/Full/";
  
  // DEFAULT
  const Int_t nSpecificStorages = 12;
  const Char_t *SpecificStorageList[nSpecificStorages][3] = {
    // path                    sim       rec
    //
    // ITS
    "ITS/Align/Data",          Ideal,    Residual,
    "ITS/Calib/SPDSparseDead", Full,     Residual,
    // TPC
    "TPC/Calib/Parameters",    Residual, Residual,
    "TPC/Calib/ClusterParam",  Residual, Residual,
    "TPC/Calib/RecoParam",     Full,     Full,
    "TPC/Calib/TimeGain",      Ideal,    Residual,
    "TPC/Calib/Correction",    Ideal,    Residual,
    "TPC/Align/Data",          Ideal,    Residual,
    "TPC/Calib/TimeDrift",     Ideal,    Residual,
    // MUON
    "MUON/Align/Data",         Ideal,    Residual,
    // ZDC
    "ZDC/Align/Data",          Ideal,    Ideal,
    "ZDC/Calib/Pedestals",     Ideal,    Ideal
  };

  for (Int_t isto = 0; isto < nSpecificStorages; isto++) {
    if (SpecificStorageList[mode][isto]) {
      printf("Setting specific storage: %s -> %s\n", SpecificStorageList[isto][0], SpecificStorageList[isto][mode+1]);
      man->SetSpecificStorage(SpecificStorageList[isto][0], SpecificStorageList[isto][mode+1]);      
    }
  }

  // temporary HACK 
  if (year < 2015) {
    man->SetSpecificStorage("AD/Calib/QAParam", "alien://Folder=/alice/data/2015/OCDB");
  }

}

