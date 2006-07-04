void TestShuttle() {

	gSystem->Load("libSHUTTLE");
        gSystem->Load("libRLDAP");
	gSystem->Load("libTest");

	AliLog::SetGlobalDebugLevel(1);

 	AliShuttleConfig config("pcalice290.cern.ch", 389, "o=alice,dc=cern,dc=ch");
        config.Print();

	AliShuttleTrigger trigger(&config);

	AliShuttle* shuttle = trigger.GetShuttle();
	TestITSPreprocessor *itsPrep = new AliITSPreprocessor("ITS",shuttle);

	trigger.CollectNew();
	
/*	TTimeStamp currentTime;
	shuttle.Process(1, currentTime.GetSec() - 18 * 3600, 
		currentTime.GetSec());*/
}


