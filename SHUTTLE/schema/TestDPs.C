void TestDPs(const char* detector, int tsFrom, int tsTo)
{
	// Query amanda for DPs for detector from timestamp tsFrom to tsTo
	//
	gSystem->Load("libRLDAP.so");
	gSystem->Load("libMonaLisa");
	gSystem->Load("libSHUTTLE");

	AliDCSClient client("alidcsamanda.cern.ch",1337,1000,500,4000);

	AliShuttleConfig config("pcalishuttle02.cern.ch", 389, "", "", "o=shuttle,dc=cern,dc=ch");

	TObjArray* list = config.GetDCSAliases(detector, 0);
	for (Int_t i=0; i<list->GetEntries(); i++)
	{
	  //Printf("%s", list->At(i)->GetName());
	  client.GetAliasValues(list, tsFrom, tsTo, i, i+1);
	}
}
