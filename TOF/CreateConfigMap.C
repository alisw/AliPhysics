void CreateConfigMap(const char* delayFlag="kFALSE", const char* startRun="0", const char* binRangeAve="13", const char* integralThr="100", const char* thrPar="0.013"){

	// Create Configuration Map entry for OCDB
	// to configure AliTOFPreprocessor to compute delays/write online calibration obj on CDB
	// USAGE: 
	// - "delayFlag" should be set to kTRUE in case the delays have to be calculated, to 
	//   kFALSE otherwise.
	// - "startRun" indicates the starting run for the online calibration object validity
	//   for delays
	
	AliTOFcalib *tofcalib = new AliTOFcalib();
	TMap *mapTOF = (TMap*)tofcalib->GetConfigMap();
	mapTOF->Add(new TObjString("ComputingDelays"),new TObjString(delayFlag));
	mapTOF->Add(new TObjString("StartingRun"),new TObjString(startRun));
	mapTOF->Add(new TObjString("BinRangeAve"),new TObjString(binRangeAve));
	mapTOF->Add(new TObjString("IntegralThr"),new TObjString(integralThr));
	mapTOF->Add(new TObjString("ThrPar"),new TObjString(thrPar));
	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	tofcalib->WriteConfigMapOnCDB("TOF/Calib");
}
