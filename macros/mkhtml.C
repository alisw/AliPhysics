void mkhtml (char *macro=0, Int_t force=0) {
// to run this macro, you must have the correct .rootrc file
// in your galice directory.
// The gAlice classes summary documentation go to directory html
// The gAlice classes source  documentation go to directory html/src
// The example macros documentation go to directory html/examples
   
  // gROOT->LoadMacro("loadlibs.C");
  // loadlibs();
  THtml html;
  TStopwatch timer;
  timer.Start();
  if(macro) {
    gROOT->LoadMacro(macro);
    html.Convert(macro,"Example Macro");
  } else {
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
    gSystem->Load("libRALICE.so");
    html.MakeAll(force,"[A-Z]*");
  }
  timer.Stop();
  timer.Print();
}
