void sim(Int_t nev=20) {
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libhijing");
  gSystem->Load("libTHijing");
  gSystem->Load("libgeant321");

  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  
  AliSimulation simulator;
  simulator.SetMakeSDigits("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0");
  //simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

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
