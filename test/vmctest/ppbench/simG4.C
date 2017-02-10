// $Id$
//
// Macro for running simulation in test/vmctest/ppbench.
// From test/ppbench. 

void simG4(Int_t nev, const char * config) {
  // Load Pythia related libraries
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libHIJING");
  gSystem->Load("libTHijing");

  // Load Geant4 libraries
  gROOT->Macro("$G4VMCINSTALL/share/examples/macro/g4libs.C");

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/commonConfig.C");

  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/genPPbenchConfig.C");
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/genExtFileConfig.C");


  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  
  AliSimulation simulator(config);
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetWriteRawData("ALL","raw.root",kTRUE); 

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  
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
