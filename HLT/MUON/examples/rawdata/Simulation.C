void Simulation(char config = "Config_MUON_test.C", Int_t nevents = 100)
{

  AliSimulation MuonSim("Config_MUON_test.C");
  MuonSim.SetMakeSDigits("MUON");
  MuonSim.SetMakeDigits("MUON");
  //  MuonSim.SetWriteRawData ("MUON");
  
  MuonSim.Run (nevents);

}
