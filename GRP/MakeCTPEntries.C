void MakeCTPEntries(){

	// Macro to put in OCDB the entries for CTP timing and scalers

	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	Char_t * filenameScalers = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.cnt");
        Char_t * filenameCTPtime = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.tim");

	AliTriggerRunScalers *scalers = AliTriggerRunScalers::ReadScalers(filenameScalers);
	AliCTPTimeParams *ctptime = AliCTPTimeParams::LoadCTPTimeParams(filenameCTPtime);

	// CTP scalers
	AliCDBMetaData* metascalers = new AliCDBMetaData();
	metascalers->SetResponsible("Roman Lietava");
	metascalers->SetComment("Dummy CTP scalers for local reconstruction");
	AliCDBId idscalers("GRP/CTP/Scalers",0,AliCDBRunRange::Infinity());
	man->Put(scalers,idscalers, metascalers);

	// CTP time parameters
	AliCDBMetaData* metactptime = new AliCDBMetaData();
	metactptime->SetResponsible("Roman Lietava");
	metactptime->SetComment("Dummy CTP time params for standalone runs");
	AliCDBId idctptime("GRP/CTP/CTPtiming",0,AliCDBRunRange::Infinity());
	man->Put(ctptime,idctptime, metactptime);

	return;
}
