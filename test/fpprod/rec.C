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
  reco.SetRecoParam("ITS",AliITSRecoParam::GetLowFluxParam());
  reco.SetRecoParam("TPC",AliTPCRecoParam::GetLowFluxParam());
  reco.SetRecoParam("TRD",AliTRDrecoParam::GetLowFluxParam());
  reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
  reco.SetRecoParam("MUON",AliMUONRecoParam::GetLowFluxParam());
  AliGRPRecoParam *grpRecoParam = AliGRPRecoParam::GetLowFluxParam();
  grpRecoParam->SetVertexerTracksConstraintITS(kFALSE);
  grpRecoParam->SetVertexerTracksConstraintTPC(kFALSE);
  reco.SetRecoParam("GRP",grpRecoParam);

  //   reco.SetInput("raw.root");
  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON VZERO T0 FMD PMD ZDC");
  reco.SetRunQA(":");
  reco.SetRunGlobalQA(kTRUE);

// **** The field map settings must be the same as in Config.C !
  AliMagWrapCheb* field = 0x0;
  field = new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);  // tracking with the real map

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
