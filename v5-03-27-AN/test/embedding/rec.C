void rec(Int_t embrun=0) {
//  new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  AliReconstruction reco;
  reco.SetWriteESDfriend(kTRUE);
  reco.SetWriteAlignmentData(kFALSE);
//    reco.SetRecoParam("ITS",AliITSRecoParam::GetHighFluxParam());
//    reco.SetRecoParam("TPC",AliTPCRecoParam::GetHighFluxParam());
//    reco.SetRecoParam("TRD",AliTRDrecoParam::GetHighFluxParam());
//    reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
//    reco.SetRecoParam("MUON",AliMUONRecoParam::GetHighFluxParam());
   
//   AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
//   AliTPCReconstructor::SetRecoParam(tpcRecoParam);
//   AliTPCReconstructor::SetStreamLevel(0);
  reco.SetRunReconstruction("ITS TPC TRD TOF VZERO");
//   reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
//   reco.SetRunQA(kFALSE);
//   reco.SetRunGlobalQA(kFALSE);

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if (embrun == 1) {
    reco.SetSpecificStorage("GRP/GRP/Data",
			    Form("local://%s/../BackgroundFull",gSystem->pwd()));
  }
  else {
    reco.SetSpecificStorage("GRP/GRP/Data",
			    Form("local://%s",gSystem->pwd()));
  }
  reco.SetRunPlaneEff(kTRUE);
  reco.SetRunQA("ALL:ALL") ;
  
  AliQA::SetQARefStorage("local://$ALICE_ROOT/OCDB") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    reco.SetQACycles(det, 999) ;
    reco.SetQAWriteExpert(det) ; 
  }
  
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
