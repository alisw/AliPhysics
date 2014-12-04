void
RunFast(Bool_t      proof=false,
	Long64_t    maxEvents=100,
	UInt_t      runNo=138190,
	Double_t    bMin=0,
	Double_t    bMax=20,
	const char* eg="default",
	Int_t       monitor=5)
{
  TString ali = gSystem->ExpandPathName("${ALICE_ROOT}");
  // TString fwd = gSystem->ExpandPathName("$ANA_SRC");
  TString fwd = ali + "/PWGLF/FORWARD/analysis2";
  gSystem->AddIncludePath(Form("-I%s/include", ali.Data()));
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

  const char* opt="";
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

  const char* url = "lite:///?workers=8";
  ::Info("runFast", "Monitor=%d", monitor);
  if (proof) FastSim::Proof(url,maxEvents, runNo, eg, bMin, bMax, monitor);
  else       FastSim::Run(maxEvents, runNo, eg, bMin, bMax, monitor);
}
