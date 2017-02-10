// $Id$
//
// Macro for running simulation in test/vmctest/gun.
// From test/gun. 

void simG4(Int_t nev, const char * config) {
  // Load Geant3 + Geant3 VMC libraries
  //
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations

  gROOT->Macro("$G4VMCINSTALL/share/examples/macro/g4libs.C");
  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/gun/commonConfig.C");

  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/gun/genExtFileConfig.C");

  AliSimulation simulator(config);
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));

  simulator.SetRunHLT("default"); // In case we do not have ancored production

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
