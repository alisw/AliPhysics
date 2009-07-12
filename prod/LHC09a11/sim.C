void sim(Int_t nev=50) {

  AliSimulation MuonSim;
  MuonSim.SetMakeTrigger("MUON");
  MuonSim.SetMakeSDigits("MUON ITS");
  MuonSim.SetMakeDigits("MUON ITS VZERO");
  MuonSim.SetWriteRawData("MUON ITS VZERO HLT","raw.root",kTRUE);
  MuonSim.SetRunHLT("libAliHLTMUON.so chains=dHLT-sim");
  MuonSim.SetRunQA("MUON:ALL");
  
  MuonSim.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"); 
  
  // QA reference
  MuonSim.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  // HLT
  MuonSim.SetSpecificStorage("HLT/ConfigMUON/DecisionComponent","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  MuonSim.SetSpecificStorage("HLT/ConfigMUON/HitReconstructor","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  MuonSim.SetSpecificStorage("HLT/ConfigMUON/MansoTrackerFSM","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  MuonSim.SetSpecificStorage("HLT/ConfigMUON/TriggerReconstructor","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  
  // CTP
  MuonSim.SetSpecificStorage("GRP/CTP/Config","alien://Folder=/alice/cern.ch/user/b/bogdan/prod2009/cdb");
  
  // alignment
  MuonSim.SetSpecificStorage("MUON/Align/Data","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  TStopwatch timer;
  timer.Start();
  
  MuonSim.Run(nev);
  
  timer.Stop();
  timer.Print();
  
}
