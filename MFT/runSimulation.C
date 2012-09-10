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
  simulator->SetRunQA("DetectorList:ActionList");
  simulator->SetRunHLT("");

  // MUON Tracker -> local:///$OCDB should reflect the content of alien://folder=/alice
  simulator->SetDefaultStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  simulator->SetSpecificStorage("MUON/Align/Data", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  simulator->SetSpecificStorage("MFT/Align/Data",  "alien://folder=/alice/cern.ch/user/a/auras/OCDB/");

  // The rest
  TStopwatch timer;
  timer.Start();
  simulator->Run(nevents);
  timer.Stop();
  timer.Print();

}
