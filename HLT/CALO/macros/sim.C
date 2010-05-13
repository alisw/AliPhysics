void sim(Int_t nev=1) {
  
  AliSimulation simulator;
  simulator.SetMakeSDigits("PHOS EMCAL");
  //  simulator.SetMakeDigitsFromHits("ITS TPC");
  //  simulator.SetWriteRawData("ALL","raw.root",kTRUE);
  simulator.SetWriteRawData("ALL");
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  
  simulator.SetRunQA(":") ; 

  
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
