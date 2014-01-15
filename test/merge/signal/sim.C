void sim(Int_t nev=1) {
   // libraries required by geant321

    gSystem->Load("liblhapdf");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libgeant321");
    
    gSystem->Load("libhijing");
    gSystem->Load("libTHijing");

  AliSimulation simulator;
  simulator.MergeWith("../backgr/galice.root",3);
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
