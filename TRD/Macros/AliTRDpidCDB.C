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
	man->SetDefaultStorage("local://$ALICE_ROOT");
	man->SetRun(0);

	AliTRDCalPIDRefMaker maker;
	maker.BuildLQReferences(file, dir);
}

//___________________________________________________________________
void generatePIDDB(const char *file = "Refs.root")
{
// Write TRD PID DB using the reference data from file "file"


	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT");
	man->SetRun(0);

	AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT");
  if (!gStorLoc) return;

  
	AliTRDCalPIDLQ *pid = new AliTRDCalPIDLQ("pid", "TRD PID object");
	pid->LoadLQReferences(file);

	AliCDBMetaData *md= new AliCDBMetaData();
  md->SetObjectClassName("AliTRDCalPIDLQ");
  md->SetResponsible("Alex Bercuci");
  md->SetBeamPeriod(1);
  md->SetAliRootVersion("v4-05-Release"); //root version
  md->SetComment("2D PID for TRD");

	gStorLoc->Put(pid, AliCDBId("TRD/Calib/PIDLQ", 0, 0), md);
}

//___________________________________________________________________
AliTRDCalPIDLQ* getPIDObject()
{
// Returns PIDLQ object.
// In order to browse histos do:
//   > AliTRDCalPIDLQ *pid = getPIDObject();
//   > pid->GetHistogram(0, 3);

	gStyle->SetOptStat(0);
	
  AliCDBManager *CDBManager = AliCDBManager::Instance();
	CDBManager->SetDefaultStorage("local://$ALICE_ROOT");
	CDBManager->SetRun(0);

	AliCDBEntry *wrap = CDBManager->Get("TRD/Calib/PIDLQ", 0);
	AliTRDCalPIDLQ *pid = dynamic_cast<const AliTRDCalPIDLQ *>wrap->GetObject();
	AliCDBMetaData *meta = wrap->GetMetaData();
	if(!pid){
		printf("Error while retriving pid object from DB.\n");
		return 0x0;
	}
	meta->PrintMetaData();
	return pid;
}

