void Shuttle(const char* param = "listen") {

	gSystem->Load("libSHUTTLE");
        gSystem->Load("$ROOTSYS/lib/libRLDAP");
	gSystem->Load("$ROOTSYS/lib/libThread");

//	AliLog::SetGlobalDebugLevel(1);

        AliShuttleConfig config("pcepalice60.cern.ch", 389,
			"cn=Shuttle,dc=alice,dc=cern,dc=ch", "passhuttle");
        config.Print();

	AliCDBStorage* cdbStorage = AliCDBManager::Instance()->
			GetStorage("local://~/temp/DCS");

	AliShuttleTrigger trigger(&config, cdbStorage);

	AliShuttle* shuttle = trigger.GetShuttle();
	// Add here detectors preprocessor ...
	//shuttle->RegisterCDBPreProcessor(new TestITSPreProcessor());

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


