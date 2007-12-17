void TestDPs(const char* detector)
{
	gSystem->Load("libRLDAP.so");
	gSystem->Load("libMonaLisa");
	gSystem->Load("libSHUTTLE");

	AliDCSClient client("alidcsamanda.cern.ch",1337,1000,500,100);

	AliShuttleConfig config("pcalishuttle01.cern.ch", 389, "", "", "o=shuttle_prod,dc=cern,dc=ch");

	TObjArray* list = config.GetDCSAliases(detector, 0);
	for (Int_t i=0; i<list->GetEntries(); i++)
	{
	  //Printf("%s", list->At(i)->GetName());
	  client.GetAliasValues(list, 1197825708, 1197825808, i, i+1);
	}
}
