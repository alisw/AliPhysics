void sim(Int_t nev=1, bool dophos = true, bool doemcal = true, bool dotm = true) {
  
  AliSimulation simulator;
  TString sdigits; 
  if(dophos) sdigits += " PHOS";
  if(doemcal) sdigits += " EMCAL";
  simulator.SetMakeSDigits(sdigits);
  
  if(dotm) simulator.SetMakeDigitsFromHits("TPC");
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
