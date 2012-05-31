void rec() {

  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/..",gSystem->pwd()));
  reco.SetSpecificStorage("VZERO/Calib/Data",
			  "local://$ALICE_ROOT/OCDB/VZERO/PbPb");

  reco.SetRunQA("ALL:ALL") ;
  
  reco.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  reco.SetRunPlaneEff(kTRUE);
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    reco.SetQACycles(det, 999) ;
    reco.SetQAWriteExpert(det) ; 
  }

  reco.SetInput("raw.root");

  reco.SetNumberOfEventsPerFile(-1); // all events in one single file

  // High multiplicity settings
  Double_t cuts[]={33, 0.1, 0.1, 0.05, 0.99, 0.9, 100}; 
  AliV0vertexer::SetDefaultCuts(cuts); 
  Double_t cts[]={33., 0.05, 0.008, 0.035, 0.1, 0.9985, 0.9,100}; 
  AliCascadeVertexer::SetDefaultCuts(cts); 

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
