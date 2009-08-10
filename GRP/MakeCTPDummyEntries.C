void MakeCTPDummyEntries(){

	// Example macro to put in OCDB the dummy entries for CTP configuration and scalers
	// The entries are at present taken from $ALICE_ROOT 
	// Should be used to test the GRP preprocessor 

	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	Char_t * filenameConfig = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.cfg");
	Char_t * filenameScalers = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.cnt");

	AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfiguration(filenameConfig);
	AliTriggerRunScalers *scalers = AliTriggerRunScalers::ReadScalers(filenameScalers);

	AliCDBMetaData* metaconfig = new AliCDBMetaData();
	metaconfig->SetResponsible("Roman Lietava");
	metaconfig->SetComment("Dummy CTP configuration for standalone runs");
	AliCDBId idconfig("GRP/CTP/DummyConfig",0,AliCDBRunRange::Infinity());
	man->Put(runcfg,idconfig, metaconfig);

	AliCDBMetaData* metascalers = new AliCDBMetaData();
	metascalers->SetResponsible("Roman Lietava");
	metascalers->SetComment("Dummy CTP scalers for standalone runs");
	AliCDBId idscalers("GRP/CTP/DummyScalers",0,AliCDBRunRange::Infinity());
	man->Put(scalers,idscalers, metascalers);
	return;
}
