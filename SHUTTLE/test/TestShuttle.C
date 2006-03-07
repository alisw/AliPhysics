void TestShuttle() {

	gSystem->Load("libSHUTTLE");
        gSystem->Load("libRLDAP");

	gSystem->Load("libTest");

	AliLog::SetGlobalDebugLevel(1);

        AliShuttleConfig config("pcepalice60.cern.ch", 389,
			"cn=Shuttle,dc=alice,dc=cern,dc=ch", "passhuttle");
        config.Print();

	AliCDBStorage* cdbStorage = AliCDBManager::Instance()->
			GetStorage("local://~/temp/DCS");

	AliShuttleTrigger trigger(&config, cdbStorage);

	AliShuttle* shuttle = trigger.GetShuttle();
	shuttle->RegisterCDBPreProcessor(new TestITSPreProcessor());

	trigger.CollectNew();
	
/*	TTimeStamp currentTime;
	shuttle.Process(1, currentTime.GetSec() - 18 * 3600, 
		currentTime.GetSec());*/
}


