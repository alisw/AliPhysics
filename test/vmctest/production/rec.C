// $Id$

void rec() {

  AliReconstruction reco;
// switch off cleanESD
  reco.SetCleanESD(kFALSE);

  
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetRunQA("kFALSE:kFALSE");
  
  reco.SetRunPlaneEff(kTRUE);
  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  reco.SetSpecificStorage("GRP/GRP/Data",
                          Form("local://%s",gSystem->pwd()));


  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
