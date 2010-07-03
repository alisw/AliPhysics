void cosmicSim(Int_t nev=10) {
  AliSimulation simulator("Config_AliGenACORDE.C");
  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetRunQA("ALL:ALL") ; 
 simulator.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    simulator.SetQACycles(det, nev+1) ;
  }
simulator.SetEventSpecie(AliRecoParam::kCosmic);
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
