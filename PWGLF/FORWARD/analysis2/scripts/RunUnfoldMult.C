void
RunUnfoldMult()
{
  TString rooUnfold = gSystem->Getenv("ROOUNFOLD");
  if (!rooUnfold.IsNull()) {
    gSystem->AddIncludePath(Form("-I%s/src", rooUnfold.Data()));
    gSystem->AddDynamicPath(rooUnfold);
  }
  gSystem->Load("libRooUnfold.so");
  gROOT->Macro("UnfoldMult.C++");
}

