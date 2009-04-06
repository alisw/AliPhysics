void CreateConfigMapNoise(const char* noiseThr="1"){

	// Create Configuration Map entry for OCDB 
	// to configure AliTOFPreprocessor to run over Noise runs
	// USAGE: 
	// - "NoiseThr" is the threshold to declare a channel as noisy 
	
	TMap *mapTOF = new TMap();
	mapTOF->Add(new TObjString("NoiseThr"),new TObjString(noiseThr));
	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	AliCDBId id("TOF/Calib/ConfigNoise",0,AliCDBRunRange::Infinity());
	AliCDBMetaData* md = new AliCDBMetaData();
	md->SetResponsible("Chiara Zampolli");
	md->SetComment("Configuration parameter for nose runs");
	man->Put(mapTOF,id,md);
}
