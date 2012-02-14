void runSimulation(Int_t seed,
		   Int_t nevents,
		   const Char_t *config,
		   Int_t runNumber) {
  
  AliSimulation *simulator = new AliSimulation(config);

  simulator->SetSeed(seed);
  simulator->SetRunNumber(runNumber);
  simulator->SetTriggerConfig("MUON");
  simulator->SetMakeDigits("MUON MFT");
  simulator->SetMakeSDigits("MUON MFT");
  simulator->SetRunQA("ALL");
  simulator->SetRunHLT("");

  // MUON Tracker
  // simulator->SetSpecificStorage("MUON/Align/Data", "local:///$OCDB/simulation/2008/v4-15-Release/Ideal");

  // The rest
  TStopwatch timer;
  timer.Start();
  simulator->Run(nevents);
  timer.Stop();
  timer.Print();

}
