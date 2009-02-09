void rec() {
  if (!strcmp(gSystem->GetBuildArch(),"win32gcc")) {
    gSystem->Load("libProof");
    gSystem->Load("libGui");
    gROOT->Macro("loadlibsrec.C");
    new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  }
  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));

  AliGRPRecoParam *grpRecoParam = AliGRPRecoParam::GetLowFluxParam();
  grpRecoParam->SetVertexerTracksConstraintITS(kFALSE);
  grpRecoParam->SetVertexerTracksConstraintTPC(kFALSE);
  reco.SetRecoParam("GRP",grpRecoParam);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
