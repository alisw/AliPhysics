

AliAnalysisTaskGammaConversion * AddTaskGammaConversion(TString arguments,AliAnalysisDataContainer *cin_esd){

  gROOT->LoadMacro("./ConfigGammaConversion.C"); // load the CreateChain macro

  ConfigGammaConversion(arguments,cin_esd);
  
  return NULL;
}
