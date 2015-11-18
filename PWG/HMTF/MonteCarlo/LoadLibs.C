void LoadLibs(Int_t tune) {
  std::cout << "Loading Libs" << std::endl;
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");// -I$ALICE_ROOT/include/pythia");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libEVGEN.so");
  gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
  gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
  gSystem->Load("liblhapdf.so");      // Parton density functions
  if (tune == -3111) {
    std::cout << "Loading epos" << std::endl;
    gSystem->Load("libpythia6.so");     // Pythia
  } 
  else if (tune < 0) {
    std::cout << "Loading phojet" << std::endl;
    
    // => phojet
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libdpmjet.so");     // 
    gSystem->Load("libTDPMjet.so");     // 
  } 
  else if (tune == 0) {
//    gSystem->Load("libpythia6.4.25.so");     // Pythia
    gSystem->Load("libpythia6_4_25.so");     // Pythia
    gSystem->Load("libEGPythia6.so");   // TGenerator interface 
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
  }
  else { //if tune>0 
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8205/xmldoc"));
    gSystem->Load("libEGPythia6");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8205");
    gSystem->Load("libAliPythia8");
  }
  gSystem->ListLibraries();
}
