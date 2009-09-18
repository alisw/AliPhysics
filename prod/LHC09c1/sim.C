void sim(Int_t nev=200) {

  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");

  // The raw data are not written due to the huge increase of the 
  // virtual memory in HLT
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
