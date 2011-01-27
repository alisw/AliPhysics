void FitSpectrum(const char* filename, const char * listName = "clambdak0Histo_00", const char* suffix = "test", Int_t iparticle) {
  
  // load basic libs, needed to 
  gROOT->LoadMacro("run.C");
  InitAndLoadLibs();

  // Load Lee's Macro
  gROOT->LoadMacro("FitControl.h+g");
  gROOT->LoadMacro("PtMassAna2.C");
  gROOT->LoadMacro("MultYields2.C");

  char* histName = 0;
  switch (iparticle) {
  case 1:
    histName = "h2PtVsMassLambda";
    break;
  default:
    cout << "Particle "<< iparticle << " to be implemented" << endl;  
  }

  TFile *file = new TFile(filename);
  TList *list = file->Get(listName); 

  TH2F * h2 = (TH2F*) list->FindObject(histName);

  TString suffixFull = histName;
  if(strlen(suffix)) suffixFull = suffixFull + "_" + suffix;
  MultYields2((TH3F*)h2,iparticle,0,suffixFull); // FIXME: modify MultYields2 to handle 1D histos

  
}
