{
  cout << "Loading " << gSystem->Getenv("PWD") 
       << "/rootlogon.C" << endl;
  cout << "Using ROOT version " 
       << gROOT->GetVersion() << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  gSystem->Load("libTree.so");
  if (gSystem->Getenv("TMPDIR")) 
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  if (1) {
      delete gRandom;
      gRandom = new TRandom3(0);
   }

  // Include directories
  // gSystem->AddIncludePath("-I\"common\" ");
  // gROOT->LoadMacro("common/Utils.C+g");
}
