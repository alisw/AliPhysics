void RetrieveDPs(const char* detector, int tsFrom, int tsTo)
{
	//
	// Query amanda for DPs for detector from timestamp tsFrom to tsTo
	//
	// date --> unix date
	// date --date='2008-02-28' '+%s'
	
	gSystem->Load("libRLDAP.so");
	gSystem->Load("libMonaLisa");
	gSystem->Load("libSHUTTLE");

	AliDCSClient client("alidcsamanda.cern.ch",1337,1000,500,100);
	//AliLog::SetClassDebugLevel("AliDCSClient",5);

	AliShuttleConfig config("pcalishuttle01.cern.ch", 389, "", "", "o=shuttle_prod,dc=cern,dc=ch");
	
	TObjArray* list = config.GetDCSAliases(detector, 0);
	TMap* map = client.GetAliasValues(list, tsFrom, tsTo);
	
	TFile* file = TFile::Open("DCSMap.root", "RECREATE");
	map->Write("DCSMap", TObject::kSingleKey);
	file->Close();
	
}
