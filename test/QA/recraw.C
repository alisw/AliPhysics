void recraw() {
  const char * kYear = "08" ; 
 
	AliLog::SetGlobalLogLevel(AliLog::kError);
	
	gSystem->Load("libRAliEn.so");
  gSystem->Load("libNet.so");
  gSystem->Load("libMonaLisa.so");
  new TMonaLisaWriter(0, "GridAliRoot-rec.C", 0, 0, "global");
  gSystem->Setenv("APMON_INTERVAL", "120");
	
	// Set the CDB storage location
	AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  //man->SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08b/OCDB/");

	// Example in case a specific CDB storage is needed
	man->SetSpecificStorage("EMCAL/*","local://DB");
  
  // ITS settings
  AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetFactorSAWindowSizes(20);
  itsRecoParam->SetClusterErrorsParam(2);
  itsRecoParam->SetFindV0s(kFALSE);
  itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
  itsRecoParam->SetUseAmplitudeInfo(kFALSE);
  // In case we want to switch off a layer
  //  itsRecoParam->SetLayerToSkip(<N>);
  //  itsRecoParam->SetLayerToSkip(4);
  //  itsRecoParam->SetLayerToSkip(5);
  itsRecoParam->SetLayerToSkip(2);
  itsRecoParam->SetLayerToSkip(3);
  //itsRecoParam->SetSAOnePointTracks();
  itsRecoParam->SetClusterMisalError(1.0); // [cm]
  itsRecoParam->SetSAUseAllClusters();
  AliITSReconstructor::SetRecoParam(itsRecoParam);
	
  // TPC settings
  //AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kFALSE);
  tpcRecoParam->SetTimeInterval(60,940);
  Double_t sysError[5]={0.3,1, 0.3/150., 1./150.,0.3/(150*150.)};
  tpcRecoParam->SetSystematicError(sysError);
  tpcRecoParam->SetMinMaxCutAbs(4.);
  tpcRecoParam->SetMinLeftRightCutAbs(6.);
  tpcRecoParam->SetMinUpDownCutAbs(6.);
  //  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
	
  // TRD setting
  // Settings for the TRD Raw Reader
  AliTRDrawStreamBase::SetRawStreamVersion("TB");
  AliTRDrawStreamTB::SetNoErrorWarning();
  AliTRDrawStreamTB::AllowCorruptedData();
  AliTRDrawStreamTB::DisableStackNumberChecker();
  AliTRDrawStreamTB::DisableStackLinkNumberChecker();
  AliTRDrawStreamTB::SetSubtractBaseline(10);
  
  // TRD reconstruction params
  AliTRDrecoParam *fTRDrecoParam = AliTRDrecoParam::GetCosmicTestParam();
  AliTRDReconstructor::SetRecoParam(fTRDrecoParam);
  AliTRDtrackerV1::SetNTimeBins(30);
	
  // PHOS settings
  AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  recEmc->SetSubtractPedestals(kTRUE);
  recEmc->SetMinE(0.05);
  recEmc->SetClusteringThreshold(0.10);
  AliPHOSReconstructor::SetRecoParamEmc(recEmc);
	
  // T0 settings
  //AliLog::SetModuleDebugLevel("T0", 10);
	
  // MUON settings
  //AliLog::SetClassDebugLevel("AliMUONRawStreamTracker",3);
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  muonRecoParam->CombineClusterTrackReco(kTRUE);
  muonRecoParam->SetCalibrationMode("NOGAIN");
  //muonRecoParam->SetClusteringMode("PEAKFIT");
  //muonRecoParam->SetClusteringMode("PEAKCOG");
  muonRecoParam->Print("FULL");
 // AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);

	// Tracking settings
	// **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps", "Maps", 2, 1., 10., AliMagFMaps::k5kG);
  Bool_t uniform = kFALSE;
  AliTracker::SetFieldMap(field, uniform);
	Double_t mostProbPt=0.35;
  AliExternalTrackParam::SetMostProbablePt(mostProbPt);
	
	// AliReconstruction settings
	AliReconstruction reco;
	reco.SetUniformFieldTracking(uniform);
	reco.SetWriteESDfriend(kTRUE);
  reco.SetWriteAlignmentData();
	reco.SetInput("raw.root");
	reco.SetUseTrackingErrorsForAlignment("ITS");
	
	// In case some detectors have to be switched off...
	reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS MUON FMD PMD T0 VZERO ZDC ACORDE");

	reco.SetRunVertexFinder(kTRUE);
	
	// all events in one single file
  reco.SetNumberOfEventsPerFile(-1);

	// switch off cleanESD
  reco.SetCleanESD(kFALSE);
	
 	reco.SetRunQA("ALL:ALL") ;
	//AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefStorage("local://$ALICE_ROOT") ;
  
	AliLog::Flush();
	
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}

