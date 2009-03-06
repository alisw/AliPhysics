void rec() {

  AliReconstruction reco;

  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetInput("raw.root");

  reco.SetNumberOfEventsPerFile(-1); // all events in one single file

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/..",gSystem->pwd()));

  reco.SetRunQA("ALL:ALL") ;
  
  reco.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    reco.SetQACycles(det, 999) ;
    reco.SetQAWriteExpert(det) ; 
  }

  TStopwatch timer;
  timer.Start();
  reco.Run();
  reco.MergeQA() ; 
  timer.Stop();
  timer.Print();
}
