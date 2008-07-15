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
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
  //   reco.SetInput("raw.root");
  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON VZERO T0 FMD PMD ZDC");
// **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps","Maps",2,1.,10.,AliMagFMaps::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
