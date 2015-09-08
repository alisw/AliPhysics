// macro to merge outputs of different QA tasks

void mergeQAgroups(const char* lst, const char* out="QAresults.root")
{
  TString lstS = lst;
  TString outS = out;
  if (lstS.IsNull()) exit(1);
  if (outS.IsNull()) exit(1);
  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass1/LoadLibraries.C");
  gSystem->Load("libCORRFW.so");
  AliFileMerger fm;
  fm.IterTXT(lstS.Data(),outS.Data());
  //  
}
