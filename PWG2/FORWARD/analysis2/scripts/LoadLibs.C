/** 
 * Load the libraries of PWG2/FORWARD/analsysis2
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
LoadLibs()
{
  gROOT->LoadClass("TVirtualMC",           "libVMC");
  gROOT->LoadClass("AliVEvent",            "libSTEERBase");
  gROOT->LoadClass("AliESDEvent",          "libESD");
  gROOT->LoadClass("AliAnalysisManager",   "libANALYSIS");
  gROOT->LoadClass("AliAnalysisTaskSE",    "libANALYSISalice");
  gROOT->LoadClass("AliAODForwardMult",    "libPWG2forward2");

#if 0
  const char* test = gSystem->GetLibraries("PWG2forward2","D",false);
  if (test && test[0] != '\0') { 
    // TInterpreter* inter = gROOT->GetInterpreter();
    // inter->ClearFileBusy();
    // inter->UnloadFile(inter->GetCurrentMacroName());
    return;
  }
  gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG2forward2");
#endif
}
//
// EOF
//
