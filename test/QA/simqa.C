void simqa()
{
	const char * kYear = "08" ; 
	gEnv->SetValue("Root.Stacktrace","no");
	gEnv->SetValue("Root.Stacktrace","no");
	AliCDBManager * man = AliCDBManager::Instance();
	//man->SetDefaultStorage("alien://Folder=/alice/simulation/2007/PDC07_v4-09-Rev-00/Ideal/CDB/");
	man->SetDefaultStorage("local://$ALICE_ROOT");
	TString detectors("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD PMD ZDC T0 VZERO"); 
	
	AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
	AliQA::SetQARefDataDirName("Sim") ; //Data, Pedestals, BlackEvent, .....
  

	AliQADataMakerSteer qas ; 
	qas.Run(detectors.Data(), AliQA::kHITS);
	qas.Run(detectors.Data(), AliQA::kSDIGITS);
	qas.Run(detectors.Data(), AliQA::kDIGITS);
}
