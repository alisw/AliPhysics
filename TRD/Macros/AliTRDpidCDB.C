// Functions to create/store/inspect TRD PID DB.
// 
// Author:
// Alex Bercuci (A.Bercuci@gsi.de)
// 


//___________________________________________________________________
void makePIDRefs(const char *dir = ".", const char *file="Refs.root")
{
// Build the reference data for PID. The simulations have to fulfill
// the directory structure defined inside AliTRDCalPIDRefMaker.
// Parameters:
// 1. "dir" - the root directory of the production
// 2. "file" - output file containing reference data saved in directory
//             "dir" 

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);

  AliTRDCalPIDRefMaker maker;
  maker.BuildLQReferences(file, dir);
}

//___________________________________________________________________
void generatePIDDB()
{
// Write TRD PID DB using the reference data from file "file"


  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("TRD/Calib/PIDLQ", "local://.");
  man->SetRun(0);

  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) return;

  AliTRDpidRefMakerLQ pidMaker;
  TObject *o = pidMaker.GetOCDBEntry("20091101");
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetObjectClassName("AliTRDCalPIDLQ");
  md->SetResponsible("Alexandru Bercuci");
  md->SetBeamPeriod(1);
  md->SetAliRootVersion("v4-16-Release"); //root version
  md->SetComment("2D PID for TRD");
  gStorLoc->Put(o, AliCDBId("TRD/Calib/PIDLQ", 0, 999999999, 0, 1), md, AliCDBManager::kReference);
}

//___________________________________________________________________
void generatePIDDBNN(const char *fileNN = "NN.root")
{
// Write TRD PID DB using the reference data from file "file"


  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("TRD/Calib/PIDNN", "local://.");
  man->SetRun(0);

  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) return;

  AliTRDCalPID *pidNN = new AliTRDCalPIDNN("pidNN", "NN TRD PID object"); 	 
  pidNN->LoadReferences(fileNN); 	 
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetObjectClassName("AliTRDCalPIDNN");
  md->SetResponsible("Alexander Wilk");
  md->SetBeamPeriod(1);
  md->SetAliRootVersion("v4-16-Release"); //root version
  md->SetComment("NN PID for TRD");
  gStorLoc->Put(pidNN, AliCDBId("TRD/Calib/PIDNN", 0, 999999999, 0, 1), md, AliCDBManager::kReference);
}

//___________________________________________________________________
AliTRDCalPID* getPIDObject(const char *method="NN")
{
// Returns PIDLQ object.
// In order to browse histos do:
//   > AliTRDCalPID *pid = getPIDObject();
//   > pid->GetHistogram(0, 3);

  gStyle->SetOptStat(0);
  
  AliCDBManager *CDBManager = AliCDBManager::Instance();
  CDBManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  CDBManager->SetRun(0);

  AliCDBEntry *wrap = CDBManager->Get(Form("TRD/Calib/PID%s", method), 0);
  AliTRDCalPID *pid = dynamic_cast<const AliTRDCalPID *>wrap->GetObject();
  AliCDBMetaData *meta = wrap->GetMetaData();
  if(!pid){
    printf("Error while retriving pid object from DB.\n");
    return 0x0;
  }
  meta->PrintMetaData();
  return pid;
}

