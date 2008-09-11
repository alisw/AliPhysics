void recraw() {
  const char * kYear = "08" ; 
 
  AliLog::SetGlobalLogLevel(AliLog::kError);
	
//	gSystem->Load("libRAliEn.so");
//  gSystem->Load("libNet.so");
// gSystem->Load("libMonaLisa.so");
//  new TMonaLisaWriter(0, "GridAliRoot-rec.C", 0, 0, "global");
// gSystem->Setenv("APMON_INTERVAL", "120");
	
	// Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  //man->SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08d/OCDB/");

	// Example in case a specific CDB storage is needed
	//man->SetSpecificStorage("EMCAL/*","local://DB");
  
	// Tracking settings
	// **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps", "Maps", 2, 1., 10., AliMagFMaps::k5kG);
  Bool_t uniform = kFALSE;
  AliTracker::SetFieldMap(field, uniform);
  Double_t mostProbPt=0.35;
  AliExternalTrackParam::SetMostProbablePt(mostProbPt);
	
	// AliReconstruction settings
  AliReconstruction reco;
  reco.SetRecoParam("TPC",AliTPCRecoParam::GetLowFluxParam());
  reco.SetRecoParam("TRD",AliTRDrecoParam::GetLowFluxParam());
  reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
  reco.SetRecoParam("MUON",AliMUONRecoParam::GetLowFluxParam());
  reco.SetUniformFieldTracking(uniform);
  reco.SetWriteESDfriend(kTRUE);
  reco.SetWriteAlignmentData();
  reco.SetInput("raw.root");
  reco.SetUseTrackingErrorsForAlignment("ITS");
	
	// In case some detectors have to be switched off.aliextr..
  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS MUON FMD PMD T0 VZERO ZDC ACORDE");
	//reco.SetRunReconstruction("ITS TRD TOF HMPID PHOS MUON FMD PMD VZERO ZDC ACORDE");

  reco.SetRunVertexFinder(kTRUE);
	
	// all events in one single file
  reco.SetNumberOfEventsPerFile(-1);

	// switch off cleanESD
  reco.SetCleanESD(kFALSE);
	
 	reco.SetRunQA("ALL:ALL") ;
	//AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefStorage("local://$ALICE_ROOT") ;
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    reco.SetQACycles(det, 9999) ;
    reco.SetQAWriteExpert(det) ; 
  }
  
  AliLog::Flush();
	
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}

