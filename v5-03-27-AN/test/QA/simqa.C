void simqa()
{
	const char * kYear = "08" ; 
	gEnv->SetValue("Root.Stacktrace","no");
	gEnv->SetValue("Root.Stacktrace","no");
	AliCDBManager * man = AliCDBManager::Instance();
	//man->SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08d/OCDB/");
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	TString detectors("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD PMD ZDC T0 VZERO"); 
	
	//AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
	AliQA::SetQARefStorage("local://$ALICE_ROOT/QAref") ;
  

	AliQAManager *  qas = AliQAManager::QAManager("sim") ; 
	qas->SetDefaultStorage("local://$ALICE_ROOT/QAref");
        qas->SetEventSpecie(AliRecoParam::kLowMult); 
	qas->Run(detectors.Data(), AliQA::kHITS);
	qas->Run(detectors.Data(), AliQA::kSDIGITS);
	qas->Run(detectors.Data(), AliQA::kDIGITS);
}
