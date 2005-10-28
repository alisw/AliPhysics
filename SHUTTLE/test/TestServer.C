void TestServer(Int_t port) {

	gSystem->Load("libSHUTTLE");
	
	gSystem->Load("libTest");

	AliLog::EnableDebug(kTRUE);
	AliLog::SetGlobalDebugLevel(1);

	TestServer server(port);

	server.Run(1000, 200);	
}
