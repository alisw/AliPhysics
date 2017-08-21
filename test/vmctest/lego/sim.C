// $Id$
//
// Macro for running lego simulation in test/vmctest/lego.

void sim(Int_t nev, const char * config, const char * det)  {

   // Load Pythia related libraries
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
 // Load Geant3 + Geant3 VMC libraries
  //
  gSystem->Load("libgeant321");
  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/lego/commonConfig.C");

  AliSimulation simulator(config);
  simulator.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  TString configFunction("Config(\"");
  configFunction.Append(det);
  configFunction.Append(TString("\");"));
  cout << configFunction.Data() << endl;
  gAlice->SetConfigFunction(configFunction.Data());
  TStopwatch timer;
  timer.Start();
  simulator.RunLego(config);
  timer.Stop();
  timer.Print();
}
