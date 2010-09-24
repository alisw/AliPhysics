{
  if (0) {
    cout << "Loading " << gSystem->Getenv("PWD") 
         << "/rootlogon.C" << endl;
    cout << "Using ROOT version " 
         << gROOT->GetVersion() << endl;
  }

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  gSystem->Load("libTree.so");
  if (gSystem->Getenv("TMPDIR")) 
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  if (1) {
      delete gRandom;
      gRandom = new TRandom3(0);
   }
  gROOT->LoadMacro("TreeClasses.C+g");
  gROOT->LoadMacro("EventPool.C+g");
  gROOT->LoadMacro("AutoCorr.C+g");
  if (0) {
    gROOT->LoadMacro("EventPoolManager.C+g");
    gROOT->LoadMacro("anaCorr+g");
  }
}



