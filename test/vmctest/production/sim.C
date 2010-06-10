void sim(Int_t nev=20) {

  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetRunHLT("");
  
  // The raw data are not written due to the huge increase of the 
  // virtual memory in HLT
  // simulator.SetWriteRawData("ALL","raw.root",kTRUE);
  
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
 
  simulator.SetSpecificStorage("GRP/GRP/Data",
                                Form("local://%s",gSystem->pwd()));
  
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
