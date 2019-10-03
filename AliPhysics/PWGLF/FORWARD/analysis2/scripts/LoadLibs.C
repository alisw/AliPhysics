void
LoadOne(const char* cls, const char* lib,
	Bool_t verbose=false,
	Bool_t pause=false)
{
  gROOT->LoadClass(cls, lib);
  if (!verbose) return;
  Info("", "Current libraries after loading %s", lib);
  gSystem->ListLibraries();
}

/** 
 * Load the libraries of PWGLF/FORWARD/analsysis2
 * 
 * @ingroup pwglf_forward_scripts
 */
void
LoadLibs(bool alsoBase=false, bool alsoHit=false)
{

  LoadOne("TVirtualMC",              "libVMC");
  LoadOne("TLorentzVector",          "libPhysics");
  LoadOne("TLinearFitter",           "libMinuit");
  LoadOne("TTree",                   "libTree");
  LoadOne("AliVEvent",               "libSTEERBase");
  LoadOne("AliESDEvent",             "libESD");
  LoadOne("AliESDEvent",             "libAOD");
  LoadOne("AliAnalysisManager",      "libANALYSIS");
  LoadOne("AliAnalysisTaskSE",       "libANALYSISalice");
  LoadOne("AliOADBPhysicsSelection", "libOADB");
  LoadOne("AliAODForwardMult",       "libPWGLFforward2");

  if (!alsoBase) return;
  LoadOne("TProof",                  "libProof");
  LoadOne("TGFrame",                 "libGui");
  LoadOne("TSAXParser",              "libXMLParser");
  LoadOne("AliCDBManager",           "libCDB");
  LoadOne("AliRawVEvent",            "libRAWDatabase");
  LoadOne("AliHit",                  "libSTEER");
  LoadOne("AliGenMC",                "libEVGEN");
  LoadOne("AliGenMC",                "libFASTSIM");
  
  // Printf("AliFMDMCTrackELoss=%p", gROOT->GetClass("AliFMDMCTrackELoss"));

  if (!alsoHit) return;
  
  LoadOne("AliFMDDigit"              "libFMDbase");
  LoadOne("AliFMDHit",               "libFMDsim");
  LoadOne("AliFMDMCHitEnergyFitter", "libPWGLFforwardhit");
  
}
//
// EOF
//
