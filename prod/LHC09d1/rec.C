void rec() {

  AliReconstruction reco;
  
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
//   reco.SetSpecificStorage("GRP/GRP/Data",
// 			  Form("local://%s",gSystem->pwd()));
  // We store the object in AliEn during the simulation
  reco.SetSpecificStorage("GRP/GRP/Data",
			  "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

  // ITS Plane efficiency
  reco.SetRunPlaneEff(kTRUE);
  AliITSRecoParam *itspar = AliITSRecoParam::GetPlaneEffParam(-1);
  itspar->SetOptTrackletsPlaneEff(kTRUE);
  reco.SetRecoParam("ITS",itspar);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
