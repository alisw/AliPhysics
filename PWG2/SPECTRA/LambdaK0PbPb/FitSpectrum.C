void FitSpectrum(const char* filename, const char * listName = "clambdak0Histo_00", const char* suffix = "test", Int_t iparticle) {
  
  // load basic libs, needed to 
  gROOT->LoadMacro("run.C");
  InitAndLoadLibs();

  // Load Lee's Macro
  gROOT->LoadMacro("AliMassFitControl.h+g");
  gROOT->LoadMacro("PtMassAna2.C");
  gROOT->LoadMacro("MultYields2.C");

  char* histName = 0;
  switch (iparticle) {
  case 0:
    histName = "h2PtVsMassK0";
    break;
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
  //MultYields3((TH3F*)h2,iparticle,0,suffixFull); // FIXME: modify MultYields2 to handle 1D histos
  MultYields2(h2,iparticle,0,suffixFull); // FIXME: modify MultYields2 to handle 1D histos

  
}
