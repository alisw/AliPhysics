// for alternative markers, add 0x20000 to the flags 
void
ExtractMCWeights(const char* nFile="dpmjet.root",
		 const char* dFile="hijing.root",
		 const char* oFile=0)
{
  // gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  // gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ANA_SRC/dndeta/tracklets3");
  gROOT->SetMacroPath(Form("%s:scripts:$(ANA_SRC)/dndeta/tracklets3",
			   gROOT->GetMacroPath()));  
  gROOT->LoadMacro("AliTrackletWeights.C+g");  
  gROOT->LoadMacro("AliTrackletAODUtils.C+g");  
  gROOT->LoadMacro("PureMCWeights.C+g");  

  PureMCWeights p;
  p.Run(nFile, dFile, oFile);
}

