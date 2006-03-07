void TestShuttleConfig() {

	gSystem->Load("libSHUTTLE");
	gSystem->Load("libRLDAP");
	
	AliShuttleConfig config("pcepalice60.cern.ch", 389,
			"cn=Shuttle,dc=alice,dc=cern,dc=ch", "passhuttle");
	config.Print();
}
