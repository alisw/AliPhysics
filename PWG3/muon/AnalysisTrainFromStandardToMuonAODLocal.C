void AnalysisTrainFromStandardToMuonAODLocal(char* filein= "AliAODs.root", 
                                             char* fileout= "AliMuonAOD.root", 
				             char* dirChain= ".",
				             char* dirData= ".",
					     Int_t nev=123456789){
     
// Macro to produce a MUON-AOD, i.e. a replica of the standard AOD, containing only events
// where at least one muon is present
// 
// - The input files are the standard AOD and the AOD.tag.root files 
// - The AOD.tag file can be:
//      1) the one previously created together with the AOD file (i.e. from
//         AnalysisTrainMuonLocal.C)
//      2) created on the fly with this macro
// - The selection of the muon events is based on the AOD tags
// - The content of the MUON-AOD can be defined by the user with some settings as
//   SetNeedsTracksBranchReplication(), SetNeedsVerticesBranchReplication() 
//   (defined in STEER/AliAODHandler.h)...


    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEER.so");              // for aliroot based analysis
    gSystem->Load("libPWG3muon.so");  // for aliroot based analysis

    // Load par files, if the analysis is par based
    // SetupPar("STEERBase");
    // SetupPar("ESD");
    // SetupPar("AOD");
    // SetupPar("ANALYSIS");
    // SetupPar("ANALYSISalice");
    // SetupPar("PWG3muon");   
     
    // Uncomment the following lines if the AOD tag file has to be created on the fly	 
    // printf("Creating AOD Tags on the fly\n");
    // AliAODTagCreator *t = new AliAODTagCreator();
    // t->SetStorage(0);
    // t->ReadLocalCollection(dirData); 
    
    AliTagAnalysis *TagAna = new AliTagAnalysis("AOD"); 
   
    // Define tag cuts to select events containing at least one muon in the dimuon spectrometer
    printf("Defining Tags cuts to select events containing at least one muon in the dimuon spectrometer\n");
    AliRunTagCuts *runCuts = new AliRunTagCuts();
    AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
    AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
    AliEventTagCuts *evCuts = new AliEventTagCuts();
    evCuts->SetNFWMuonRange(1,10);
   
    // Create the chain of interesting events
    TChain* chain = 0x0;
    TagAna->ChainLocalTags(dirChain);
    TagAna->SetType("AOD");
    chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
    Info("AnalysisTrainFromStandardToMuonAOD",Form("CHAIN HAS %d ENTRIES",(Int_t)chain->GetEntries()));

    // Define aod input handler
    AliAODInputHandler* aodInputHandler = new AliAODInputHandler();
    
    // Define aod output handler
    AliAODHandler* aodOutputHandler = new AliAODHandler();

    // Create non standard AOD
    aodOutputHandler->SetCreateNonStandardAOD(); 
    
    // Select the branches to be replicated in the MUON-AOD
    aodOutputHandler->SetNeedsHeaderReplication();
    aodOutputHandler->SetNeedsTracksBranchReplication();
    aodOutputHandler->SetNeedsVerticesBranchReplication();
    aodOutputHandler->SetNeedsV0sBranchReplication();
    aodOutputHandler->SetNeedsTrackletsBranchReplication(); 
    aodOutputHandler->SetNeedsPMDClustersBranchReplication();
    aodOutputHandler->SetNeedsJetsBranchReplication();
    aodOutputHandler->SetNeedsFMDClustersBranchReplication();
    aodOutputHandler->SetNeedsCaloClustersBranchReplication();
    
    aodOutputHandler->SetOutputFileName(fileout);

     // Define the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("AOD Manager", "AOD Manager");
    mgr->SetInputEventHandler(aodInputHandler);
    mgr->SetOutputEventHandler(aodOutputHandler);
    //mgr->SetDebugLevel(10);
    
    AliAnalysisTaskFromStandardToMuonAOD *aodfilter = new AliAnalysisTaskFromStandardToMuonAOD("AOD Filter");
    mgr->AddTask(aodfilter);
  
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
							       
    mgr->ConnectInput(aodfilter,0,cinput1);
    mgr->ConnectOutput(aodfilter,0,coutput1);
 
    // Run the analysis    
    printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain,nev);
}


Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }

    gSystem->ChangeDirectory("../");
  }

  return 1;
}

//______________________________________________________________________________
void SetupPar(char* pararchivename)
{
    if (pararchivename) {
	char processline[1024];
	sprintf(processline,".! tar xvzf %s.par",pararchivename);
	gROOT->ProcessLine(processline);
	TString ocwd = gSystem->WorkingDirectory();
	gSystem->ChangeDirectory(pararchivename);
	
	// check for BUILD.sh and execute
	if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
	    printf("*******************************\n");
	    printf("*** Building PAR archive    ***\n");
	    printf("*******************************\n");
	    
	    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
		Error("runProcess","Cannot Build the PAR Archive! - Abort!");
		return -1;
	    }
	}
	// check for SETUP.C and execute
	if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
	    printf("*******************************\n");
	    printf("*** Setup PAR archive       ***\n");
	    printf("*******************************\n");
	    gROOT->Macro("PROOF-INF/SETUP.C");
	}
	
	gSystem->ChangeDirectory(ocwd.Data());
   printf("Current dir: %s\n", ocwd.Data());
    } 
}
