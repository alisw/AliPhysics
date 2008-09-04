void AnalysisTrainFromStandardToMuonAODLocal(char* filein= "AliAODs.root", 
                                             char* fileout= "AliMuonAOD.root", 
				             char* dirChain= ".",
				             char* dirData= "."){
     
// Macro to produce a MUON-AOD, i.e. a replica of the standard AOD, containing only events
// where at least one muon is present
// 
// - The input files are the ESD (used only for tag creation) and the standard AOD.
//   Tags files are created from all the ESD/AOD files placed in the directory dirData 
//   and in its subdirectories.
// - The selection of the muon events is based on the AOD tags.
// - The content of the MUON-AOD can be defined by the user with some settings as
//   SetNeedsTracksBranchReplication(), SetNeedsVerticesBranchReplication() 
//   (defined in STEER/AliAODHandler.h)...


    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");

    // Load par files  
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    SetupPar("PWG3base");      
    SetupPar("PWG3muon");   
        
    // Create ESD tag files   
    printf("Creating ESD Tags\n");
    AliESDTagCreator *tesd = new AliESDTagCreator();
    tesd->SetStorage(0);     
    tesd->ReadLocalCollection(dirData);  
 
    // Create AOD tag files   
    printf("Creating AOD Tags\n");
    AliAODTagCreator *t = new AliAODTagCreator();
    t->SetStorage(0);
    t->ReadLocalCollection(dirData);    
    
    AliTagAnalysis *TagAna = new AliTagAnalysis("AOD"); 
    
    // Define tag cuts to select only events with muons
    printf("Defining Tags cuts to select events containing at least one muon\n");
    AliRunTagCuts *runCuts = new AliRunTagCuts();
    AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
    AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
    AliEventTagCuts *evCuts = new AliEventTagCuts();
    evCuts->SetNMuonRange(1,10);
    
    // Create the chain of interesting events
    TChain* chain = 0x0;
    TagAna->ChainLocalTags(dirChain);
    chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);

    // Define aod input handler
    AliAODInputHandler* aodInputHandler = new AliAODInputHandler();
    
    // Define aod output handler
    AliAODHandler* aodOutputHandler = new AliAODHandler();

    // Create non standard AOD
    aodOutputHandler->SetCreateNonStandardAOD(); 
    
    // Select the branches to be written in the MUON-AOD
    aodOutputHandler->SetNeedsHeaderReplication();
    aodOutputHandler->SetNeedsTracksBranchReplication();
    aodOutputHandler->SetNeedsVerticesBranchReplication();
    //aodOutputHandler->SetNeedsV0sBranchReplication();
    //aodOutputHandler->SetNeedsTrackletsBranchReplication(); 
    //aodOutputHandler->SetNeedsPMDClustersBranchReplication();
    //aodOutputHandler->SetNeedsJetsBranchReplication();
    //aodOutputHandler->SetNeedsFMDClustersBranchReplication();
    //aodOutputHandler->SetNeedsCaloClustersBranchReplication();
    
    aodOutputHandler->SetOutputFileName(fileout);

     // Define the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("AOD Manager", "AOD Manager");
    mgr->SetInputEventHandler(aodInputHandler);
    mgr->SetOutputEventHandler(aodOutputHandler);
    //mgr->SetDebugLevel(10);
    
    AliAnalysisTaskFromStandardToMuonAOD *aodfilter = new AliAnalysisTaskFromStandardToMuonAOD("AOD Filter");
    mgr->AddTask(aodfilter);
  
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");

    mgr->ConnectInput(aodfilter,0,cinput1);
    mgr->ConnectOutput(aodfilter,0,coutput1);
    
    // Run the analysis    
    printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
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
