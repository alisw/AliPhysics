void sim(Int_t nev=1) {

  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libgeant321");
  
  gSystem->Load("libHIJING");
  gSystem->Load("libTHijing");

  AliSimulation simulator;
  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetSpecificStorage("VZERO/Calib/Data",
			       "local://$ALICE_ROOT/OCDB/VZERO/PbPb");
  simulator.SetRunQA("ALL:ALL") ; 
  simulator.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    simulator.SetQACycles((AliQAv1::DETECTORINDEX_t)det, nev+1) ;
  }

  simulator.SetRunHLT("default"); // In case we do not have ancored production

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
