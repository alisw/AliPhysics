// $Id$
//
// Macro for running lego simulation in test/vmctest/lego.

void simG4(Int_t nev, const char * config, const char * det)  {

  // Load Pythia related libraries
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  // Load Geant4 + Geant4 VMC libraries
  //
  gROOT->Macro("$G4VMCINSTALL/share/examples/macro/g4libs.C");
  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/lego/commonConfig.C");

  AliSimulation simulator(config);
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
