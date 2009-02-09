void recraw() {

	// AliReconstruction settings
  AliReconstruction reco;

  reco.SetWriteESDfriend(kTRUE);
  reco.SetWriteAlignmentData();
  reco.SetInput("raw.root");

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/..",gSystem->pwd()));


  reco.SetUseTrackingErrorsForAlignment("ITS");
	
  // all events in one single file
  reco.SetNumberOfEventsPerFile(-1);

  // switch off cleanESD
  reco.SetCleanESD(kFALSE);
	
  reco.SetRunQA("ALL:ALL") ;

  AliQA::SetQARefStorage("local://$ALICE_ROOT/OCDB") ;
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

