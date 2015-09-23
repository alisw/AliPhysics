// Reconstruction steering macro

// mikael.mieskolainen@cern.ch, 22.9.2015

void rec() {

  gSystem->Load("libEVGEN"); // Needs to be!

  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

//  reco.SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  reco.SetRunPlaneEff(kTRUE);

  // Quality Assurance
  reco.SetRunQA("ALL:ALL") ;
  reco.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    reco.SetQACycles((AliQAv1::DETECTORINDEX_t)det, 999);
    reco.SetQAWriteExpert((AliQAv1::DETECTORINDEX_t)det);
  }

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
