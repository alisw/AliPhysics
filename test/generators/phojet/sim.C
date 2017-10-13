void sim(Int_t nev=5) {
  // Libraries required by the simulation

  gSystem->Load("liblhapdf"); // Parton density functions
  gSystem->Load("libEGPythia6"); // TGenerator interface
  gSystem->Load("libgeant321");

  gSystem->Load("libDPMJET"); // DPMJET, PhoJet and Pythia6115 library
  gSystem->Load("libAliPythia6"); // ALICE specific implementations
  gSystem->Load("libTDPMjet"); // DPMJET interface


  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
 
  simulator.SetRunHLT("default"); // In case we do not have ancored production

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
