Bool_t Shuttle(const char* param = "listen") {

	// WARNING: if ldap is built with ssl support it may cause confilcts with the 
	// AliEn interface. If this happens, grid storage activation must be done BEFORE 
	// loading LDAP libraries!!!

	gSystem->Load("libRLDAP.so");
	gSystem->Load("libSHUTTLE");
	gSystem->Load("$ROOTSYS/lib/libThread");
	gSystem->Load("$ALICE_ROOT/SHUTTLE/test/libTest.so");

//	AliLog::SetGlobalDebugLevel(1);

	// Setting local CDB and reference storage locations
	AliShuttle::SetMainCDB("alien://user=aliprod?folder=colla/GridShuttleCDB");
	AliShuttle::SetMainRefStorage("alien://user=aliprod?folder=colla/GridShuttleRefStorage");

//	AliShuttle::SetMainCDB("local://testLeakCDB");
//	AliShuttle::SetMainRefStorage("local://testLeakRef");

	AliShuttle::SetLocalCDB("local://LocalShuttleCDB");
	AliShuttle::SetLocalRefStorage("local://LocalShuttleRefStorage");

	AliShuttle::SetProcessDCS(kFALSE);


//	AliCDBManager *man = AliCDBManager::Instance();
//	man->SetDefaultStorage("local://MainCDB");
//	man->SetDefaultStorage("alien://DBFolder=ShuttleMainCDB");


	AliShuttleConfig config("pcalice290.cern.ch", 389, "o=alice,dc=cern,dc=ch");
	config.SetProcessAll(kTRUE);
        config.Print();

	AliShuttleTrigger trigger(&config);

	AliShuttle* shuttle = trigger.GetShuttle();

	// Add here detectors preprocessor ...
	TestTPCPreprocessor *tpcPrep = new TestTPCPreprocessor(shuttle);
	TestITSPreprocessorSPD *spdPrep = new TestITSPreprocessorSPD("SPD",shuttle);
	TestRICHPreprocessor *richPrep = new TestRICHPreprocessor("HMP",shuttle);
	TestZDCPreprocessor *zdcPrep = new TestZDCPreprocessor("ZDC",shuttle);

	TString paramStr(param);
	
	if (paramStr.IsDigit()) {
		Int_t run = paramStr.Atoi();
		trigger.Collect(run);
	} else if (paramStr == "new") {
		Bool_t result = trigger.Collect();
	} else if (paramStr == "listen") {
		trigger.Run();
	} else {
		cout<<"Bad parameter: "<<param<<endl;
		cout<<"Parameter options: "<<endl;
		cout<<"<run> - collect data for the given run"<<endl;
		cout<<"new - collect data only for the new runs"<<endl;
		cout<<"listen - start listening for DAQ notification"<<endl;
		cout<<"<empty parameter> - the same as 'listen'"<<endl;
	}

	AliCDBManager::Destroy();
}


