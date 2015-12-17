void
RunEPos(const char* url="lite://${PWD}/index.root?events=-1&run=138190#Particle", 
	const char* opt="")
{
  TString fwd = ""; // gSystem->Getenv("ANA_SRC");
  if (fwd.IsNull()) 
    fwd = gSystem->ExpandPathName("${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2");
  gSystem->AddIncludePath(Form("-I${ALICE_ROOT}/include "
			       "-I${ALICE_PHYSICS}/include "
			       "-I%s/include",
			       fwd.Data()));
  gROOT->SetMacroPath(Form("%s:%s/sim", gROOT->GetMacroPath(), fwd.Data()));

  // Remember to copy changes to FastSim.C(FastSim::ProofLoadLibs)
  TList clsLib;
  clsLib.Add(new TNamed("TVirtualMC",              "libVMC"));
  clsLib.Add(new TNamed("TLorentzVector",          "libPhysics"));
  clsLib.Add(new TNamed("TLinearFitter",           "libMinuit"));
  clsLib.Add(new TNamed("TTree",                   "libTree"));
  clsLib.Add(new TNamed("TProof",                  "libProof"));
  clsLib.Add(new TNamed("TGFrame",                 "libGui"));
  clsLib.Add(new TNamed("TSAXParser",              "libXMLParser"));
  clsLib.Add(new TNamed("AliVEvent",               "libSTEERBase"));
  clsLib.Add(new TNamed("AliESDEvent",             "libESD"));
  clsLib.Add(new TNamed("AliAODEvent",             "libAOD"));
  clsLib.Add(new TNamed("AliAnalysisManager",      "libANALYSIS"));
  clsLib.Add(new TNamed("AliCDBManager",           "libCDB"));
  clsLib.Add(new TNamed("AliRawVEvent",            "libRAWDatabase"));
  clsLib.Add(new TNamed("AliHit",                  "libSTEER"));
  clsLib.Add(new TNamed("AliGenMC",                "libEVGEN"));
  clsLib.Add(new TNamed("AliFastEvent",            "libFASTSIM"));

  TIter next(&clsLib);
  TObject* obj = 0;
  while ((obj = next())) {
    gROOT->LoadClass(obj->GetName(), obj->GetTitle());
  }


  // Uncomment next line to use number of diffractive processes for SD
  // detection.
  gSystem->AddIncludePath("-DNO_DPMJET_TYPE");
  // gDebug = 7;
  gROOT->LoadMacro(Form("%s/sim/FastShortHeader.C", fwd.Data()));
  gROOT->LoadMacro(Form("%s/sim/FastCentEstimators.C+%s",fwd.Data(),opt));
  gROOT->LoadMacro(Form("%s/sim/FastMonitor.C+%s",fwd.Data(),opt));
  gROOT->LoadMacro(Form("%s/sim/FastSim.C+%s",fwd.Data(),opt));

  const char* cleanFiles[] = { "grp.dat",
			       "galice.root",
			       "Kinematics.root",
			       "fort.8",
			       "fort.16",
			       0 };
  const char** pClean = cleanFiles;
  while (*pClean) { 
    gSystem->Unlink(*pClean);
    pClean++;
  }
  EPosSim::Run(url, opt);
}
