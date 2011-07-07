/** 
 * Load the libraries of PWG2/FORWARD/analsysis2
 * 
 * @ingroup pwg2_forward_scripts
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
}
//
// EOF
//
