void rec() {

  AliReconstruction MuonRec;

  MuonRec.SetInput("raw.root");
  MuonRec.SetRunLocalReconstruction("MUON ITS VZERO FMD");
  MuonRec.SetRunTracking("MUON ITS VZERO FMD");
  MuonRec.SetRunVertexFinder(kTRUE);
  MuonRec.SetFillESD("MUON VZERO FMD HLT");
  MuonRec.SetRunQA("MUON:ALL");
  
  MuonRec.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

  // QA reference
  MuonRec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // HLT
  MuonRec.SetSpecificStorage("HLT/ConfigMUON/DecisionComponent","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  MuonRec.SetSpecificStorage("HLT/ConfigMUON/HitReconstructor","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  MuonRec.SetSpecificStorage("HLT/ConfigMUON/MansoTrackerFSM","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  MuonRec.SetSpecificStorage("HLT/ConfigMUON/TriggerReconstructor","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  
  // CTP
  MuonRec.SetSpecificStorage("GRP/CTP/Config","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
    
  // reconstruction parameters
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  muonRecoParam->SaveFullClusterInESD(kTRUE,100.);
  //muonRecoParam->Print("FULL");
  MuonRec.SetRecoParam("MUON",muonRecoParam);

  MuonRec.SetNumberOfEventsPerFile(500);

  TStopwatch timer;
  timer.Start();
  MuonRec.Run();
  timer.Stop();
  timer.Print();
  
}
