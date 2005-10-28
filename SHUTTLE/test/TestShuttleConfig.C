void TestShuttleConfig() {

	gSystem->Load("libSHUTTLE");
	gSystem->Load("libRLDAP");
	
	AliShuttleConfig config("localhost", 5000);
	config.Print();
}
