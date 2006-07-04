void Shuttle(const char* param = "listen") {

	gSystem->Load("libSHUTTLE");
        gSystem->Load("$ROOTSYS/lib/libRLDAP");
	gSystem->Load("$ROOTSYS/lib/libThread");
	gSystem->Load("test/libTest.so");

//	AliLog::SetGlobalDebugLevel(1);

	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://MainCDB");

	AliShuttleConfig config("pcalice290.cern.ch", 389,
			"o=alice,dc=cern,dc=ch");
	config.SetProcessAll(kTRUE);
        config.Print();

	AliShuttleTrigger trigger(&config);

	AliShuttle* shuttle = trigger.GetShuttle();

	// Add here detectors preprocessor ...
	TestTPCPreprocessor *tpcPrep = new TestTPCPreprocessor("TPC",shuttle);
	TestITSPreprocessor *itsPrep = new TestITSPreprocessor("ITS",shuttle);

	TString paramStr(param);
	
	if (paramStr.IsDigit()) {
		Int_t run = paramStr.Atoi();
		trigger.Collect(run);
	} else if (paramStr == "new") {
		trigger.CollectNew();
	} else if (paramStr == "all") {
		trigger.CollectAll();
	} else if (paramStr == "listen") {
		trigger.Run();
	} else {
		cout<<"Bad parameter: "<<param<<endl;
		cout<<"Parameter options: "<<endl;
		cout<<"<run> - collect data for the given run"<<endl;
		cout<<"new - collect data only for the new runs"<<endl;
		cout<<"all - collect data for all runs"<<endl;
		cout<<"listen - start listening for DAQ notification"<<endl;
		cout<<"<empty parameter> - the same as 'listen'"<<endl;
	}

	AliCDBManager::Destroy();
}


