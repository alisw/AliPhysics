Bool_t Shuttle(const char* param = "listen", const char* dets=0) {

	// WARNING: if ldap is built with ssl support it may cause confilcts with the 
	// AliEn interface. If this happens, grid storage activation must be done BEFORE 
	// loading LDAP libraries!!!

	gSystem->Load("libRAliEn.so");
	gSystem->Load("libRLDAP.so");
	gSystem->Load("libMonaLisa");
	gSystem->Load("libSHUTTLE");
	gSystem->Load("libThread");
//	gSystem->Load("$ALICE_ROOT/SHUTTLE/test/libTest.so");

	AliLog::SetGlobalDebugLevel(2);
	
	const char* pwd = gSystem->pwd();
	
	if( !gSystem->Getenv("SHUTTLE_DIR")) 
	{
		printf("Setting base SHUTTLE folder to %s\n", pwd);
		gSystem->Setenv("SHUTTLE_DIR", pwd);
	} else {
		printf("Shuttle base folder is %s\n", gSystem->Getenv("SHUTTLE_DIR"));
	}

// 	Setting local CDB and reference storage locations
	  AliShuttle::SetMainCDB("alien://user=aliprod?folder=testShuttle/OCDB");
	  AliShuttle::SetMainRefStorage("alien://user=aliprod?folder=testShuttle/Reference");

//        AliShuttle::SetMainCDB("local://$SHUTTLE_DIR/LocalShuttleCDB");
//        AliShuttle::SetMainRefStorage("local://$SHUTTLE_DIR/LocalShuttleRefStorage");

	AliShuttle::SetLocalCDB("local://$SHUTTLE_DIR/LocalShuttleCDB");
	AliShuttle::SetLocalRefStorage("local://$SHUTTLE_DIR/LocalShuttleRefStorage");

// 	Setting Shuttle log and temp folders
	AliShuttle::SetShuttleLogDir("$SHUTTLE_DIR/log");
	AliShuttle::SetShuttleTempDir("$SHUTTLE_DIR/temp");
	
	
	
	AliShuttle::SetProcessDCS(kTRUE);

//	AliShuttleConfig config("pcalice290.cern.ch", 389, "o=alice,dc=cern,dc=ch");
	AliShuttleConfig config("pcalishuttle01.cern.ch", 389, "", "", "o=alice,dc=cern,dc=ch");
	config.SetProcessAll(kTRUE);
        config.Print();

	AliShuttleTrigger trigger(&config, 100000);

	AliShuttle* shuttle = trigger.GetShuttle();

	// Add here detectors preprocessor ...

	TString detector = dets;
	
	printf ("Processing detectors: %s \n", detector.Data());

	if (detector.Contains("SPD")) 
		new AliITSPreprocessorSPD(shuttle);
	if (detector.Contains("SDD")) 
		new AliITSPreprocessorSDD(shuttle);
	if (detector.Contains("SSD")) 
		new AliITSPreprocessorSSD(shuttle);
	if (detector.Contains("TPC")) 
		new AliTPCPreprocessor(shuttle);
	if (detector.Contains("TRD")) 
		new AliTRDPreprocessor(shuttle);
	if (detector.Contains("TOF")) 
		new AliTOFPreprocessor(shuttle);
	if (detector.Contains("PHS")) 
	{
		gSystem->Load("libPHOSshuttle");
		new AliPHOSPreprocessor(shuttle);
	}
	if (detector.Contains("CPV")) 
		new AliCPVPreprocessor(shuttle);
	if (detector.Contains("HMP")) 
		new AliHMPIDPreprocessor(shuttle);
	if (detector.Contains("EMC")) 
		new AliEMCALPreprocessor(shuttle);
	if (detector.Contains("MCH")) 
		new AliMUONPreprocessor(shuttle);
	if (detector.Contains("MTR")) 
		new AliMTRPreprocessor(shuttle);
	if (detector.Contains("FMD")) 
		new AliFMDPreprocessor(shuttle);
	if (detector.Contains("ZDC")) 
		new AliZDCPreprocessor(shuttle);
	if (detector.Contains("PMD")) 
		new AliPMDPreprocessor("PMD", shuttle);
	if (detector.Contains("T00")) 
	{
		gSystem->Load("libT0shuttle");
		new AliT0Preprocessor("T00", shuttle);
	}
	if (detector.Contains("V00")) 
		new AliVZEROPreprocessor(shuttle);
	if (detector.Contains("GRP")) 
		new AliGRPPreprocessor(shuttle);

//	AliTOFPreprocessor *tofPrep = new AliTOFPreprocessor(shuttle);
//	AliTRDPreprocessor *trdPrep = new AliTRDPreprocessor(shuttle);
//	AliGRPPreprocessor *grpPrep = new AliGRPPreprocessor(shuttle);
	
	TString paramStr(param);
	
	TStopwatch stopwatch;
	stopwatch.Start();

      if (paramStr.IsDigit() || paramStr == "-1") {
	      Int_t run = paramStr.Atoi();
	      trigger.Collect(run);
      } else if (paramStr == "new") {
	      Bool_t result = trigger.Collect();
      } else if (paramStr == "listen") {
	      trigger.Run();
      } else {
	      cout<<"Bad parameter: "<<param<<endl;
	      cout<<"Parameter options: "<<endl;
	      cout<<"<run> - collect data for the given run"<<endl;
	      cout<<"new - collect data only for the new runs"<<endl;
	      cout<<"listen - start listening for DAQ notification"<<endl;
	      cout<<"<empty parameter> - the same as 'listen'"<<endl;
      }

	printf("Execution time: R:%.2fs C:%.2fs \n",
	       stopwatch.RealTime(),stopwatch.CpuTime());

	AliCDBManager::Destroy();
}


