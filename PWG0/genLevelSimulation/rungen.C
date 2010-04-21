typedef enum {kPhojet = -1, kPyTuneCDFA=100,kPyTuneAtlasCSC=306, kPyTuneCMS6D6T=109, kPyTunePerugia0=320 } Tune_t;

void rungen(Tune_t tune = kPyTuneCDFA, Float_t energy, Int_t nev=1, TString process, Int_t random_index = 0){
  // Simulation and reconstruction
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");
  gSystem->Load("liblhapdf.so");      // Parton density functions
  if (tune == kPhojet) {
    cout << "Loading phojet" << endl;
    
    // => phojet
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libdpmjet.so");     // 
    gSystem->Load("libTDPMjet.so");     // 
  } 
  else if (tune == kPyTunePerugia0) {
    gSystem->Load("libEGPythia6.so");   // TGenerator interface 
    gSystem->Load("libpythia6.4.21.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
  }

  else {
    gSystem->Load("libEGPythia6.so");   // TGenerator interface 
    gSystem->Load("libqpythia.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
   }
  gROOT->LoadMacro("fastGen.C");
  fastGen(tune, energy, nev, process);


  timer.Stop();
  timer.Print();
}
