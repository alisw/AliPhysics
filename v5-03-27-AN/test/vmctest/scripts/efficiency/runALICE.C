void runALICE(TString mode="local",
	      TString analysisMode="full",
	      Bool_t useMC = kTRUE,
	      Int_t nEvents=1.0*1e6, 
	      Int_t nEventsSkip=0,
	      TString format="esd") {
  
  gSystem->Load("libTree.so"); 
  gSystem->Load("libGeom.so"); 
  gSystem->Load("libVMC.so"); 
  gSystem->Load("libPhysics.so"); 
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");

  
  // Create and configure the alien handler plugin 
  gROOT->LoadMacro("SetupAnalysis.C"); 
  SetupAnalysis(mode.Data(),
		analysisMode.Data(),
		useMC,
		nEvents,
		nEventsSkip,
		format.Data());
  
}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);
                                                                                                                                               
    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      //       printf("*******************************\n");
      //       printf("*** Building PAR archive    ***\n");
      //       printf("*******************************\n");
      
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      //       printf("*******************************\n");
      //       printf("*** Setup PAR archive       ***\n");
      //       printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory("../");
  }
  return 1;
}
