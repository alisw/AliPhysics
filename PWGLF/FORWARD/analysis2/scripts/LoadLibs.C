/** 
 * Load the libraries of PWGLF/FORWARD/analsysis2
 * 
 * @ingroup pwglf_forward_scripts
 */
void
LoadLibs()
{

  gROOT->LoadClass("TVirtualMC",              "libVMC");
  gROOT->LoadClass("TLorentzVector",          "libPhysics");
  gROOT->LoadClass("TLinearFitter",           "libMinuit");
  gROOT->LoadClass("TTree",                   "libTree");
  gROOT->LoadClass("AliVEvent",               "libSTEERBase");
  gROOT->LoadClass("AliESDEvent",             "libESD");
  gROOT->LoadClass("AliESDEvent",             "libAOD");
  gROOT->LoadClass("AliAnalysisManager",      "libANALYSIS");
  gROOT->LoadClass("AliAnalysisTaskSE",       "libANALYSISalice");
  gROOT->LoadClass("AliOADBPhysicsSelection"  "libOADB");
  gROOT->LoadClass("AliAODForwardMult",       "libPWGLFforward2");
}
//
// EOF
//
