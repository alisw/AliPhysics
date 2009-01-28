void sim(Int_t nev=1) {
  AliSimulation simulator;
  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetRunQA("ALL:ALL") ; 
  AliQA::SetQARefStorage("local://$ALICE_ROOT") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    simulator.SetQACycles(det, nev+1) ;
  }
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
