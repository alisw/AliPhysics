void sim(Int_t nev=200) {
  AliSimulation simulator;
  simulator.SetRunHLT("");
  simulator.SetMakeSDigits("VZERO");
  simulator.SetMakeDigitsFromHits("ITS");
  simulator.SetMakeDigits("ITS VZERO");
  
//  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

//   simulator.SetSpecificStorage("GRP/GRP/Data",
// 			       Form("local://%s",gSystem->pwd()));

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
