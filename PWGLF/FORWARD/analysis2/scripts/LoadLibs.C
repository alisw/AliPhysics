/** 
 * Load the libraries of PWGLF/FORWARD/analsysis2
 * 
 * @ingroup pwglf_forward_scripts
 */
void
LoadLibs(bool alsoHit=false)
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

  if (!alsoHit) return;
  
  gROOT->LoadClass("TProof",                  "libProof");
  gROOT->LoadClass("TGFrame",                 "libGui");
  gROOT->LoadClass("AliCDBManager",           "libCDB");
  gROOT->LoadClass("AliRawVEvent",            "libRAWDatabase");
  gROOT->LoadClass("AliHit",                  "libSTEER");
  gROOT->LoadClass("AliFMDDigit"              "libFMDbase");
  gROOT->LoadClass("AliFMDHit",               "libFMDsim");
  gROOT->LoadClass("AliFMDMCHitEnergyFitter", "libPWGLFforwardhit");
  
}
//
// EOF
//
