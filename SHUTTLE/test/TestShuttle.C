void TestShuttle() {

	gSystem->Load("libSHUTTLE");
        gSystem->Load("libRLDAP");

	gSystem->Load("libTest");

	AliLog::SetGlobalDebugLevel(1);

        AliShuttleConfig config("localhost", 5000);
        config.Print();

	AliShuttle shuttle(&config, "local://~/temp/DCS");
	shuttle.RegisterCDBPreProcessor(new TestITSPreProcessor());
	
	TTimeStamp currentTime;
	shuttle.Process(1, currentTime.GetSec() - 18 * 3600, 
		currentTime.GetSec());
}


