void TestShuttleConfig() {

	gSystem->Load("libSHUTTLE");
	gSystem->Load("libRLDAP");

	AliShuttleConfig config("pcalice290.cern.ch", 389, "o=alice,dc=cern,dc=ch");
	config.Print();
}
