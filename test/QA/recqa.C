void recqa()
{
	const char * kYear = "08" ; 
 	gEnv->SetValue("Root.Stacktrace","no");
	gEnv->SetValue("Root.Stacktrace","no");
	AliCDBManager * man = AliCDBManager::Instance();
	//man->SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08d/OCDB/") ; 
	man->SetDefaultStorage("local://$ALICE_ROOT");
  man.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
	TString detectors("ITS TPC TRD TOF PHOS HMPID EMCAL/*MUON*/ FMD ZDC PMD T0 VZERO");
	
	//AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
	AliQA::SetQARefStorage("local://$ALICE_ROOT") ;
	//AliQA::SetQARefDataDirName(AliQA::kMONTECARLO) ; //RUN_TYPE
   
  AliQADataMakerSteer qas("rec") ; 
  qas.Run(detectors.Data(), AliQA::kRECPOINTS);
}
