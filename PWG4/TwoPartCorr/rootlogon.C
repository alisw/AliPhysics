{
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
}



