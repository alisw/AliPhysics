void AnalysisTrainMuonCompLocal()
{
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  // Common packages
  SetupPar("STEERBase");
  SetupPar("ESD");
  SetupPar("AOD");
  SetupPar("ANALYSIS");
  SetupPar("ANALYSISalice");
  // Analysis-specific packages
  SetupPar("PWG3muon");      
  
  gROOT->LoadMacro("AliAnalysisTaskAODvsESD.cxx+");  

  AliTagAnalysis *TagAna = new AliTagAnalysis("ESD"); 
  TagAna->ChainLocalTags("./");

  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  
  TChain* chain = 0x0;
  chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);
  
  // Define the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Analysis train");
  
  // ESD input handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdHandler);
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("AOD.root");  
  mgr->SetOutputEventHandler(aodHandler);
  
  // Set of cuts for the ESD filter
  // standard cuts
  AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
  esdTrackCutsL->SetMinNClustersTPC(50);
  esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsL->SetRequireTPCRefit(kTRUE);
  esdTrackCutsL->SetMinNsigmaToVertex(3);
  esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
  // hard cuts
  AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
  esdTrackCutsH->SetMinNClustersTPC(100);
  esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
  esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsH->SetRequireTPCRefit(kTRUE);
  esdTrackCutsH->SetMinNsigmaToVertex(2);
  esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
  
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsL);
  trackFilter->AddCuts(esdTrackCutsH);
  
  // first task - ESD filter task putting standard info in the output generic AOD 
  AliAnalysisTaskESDfilter *esdFilter = new AliAnalysisTaskESDfilter("ESD Filter");
  //esdFilter->SetTrackFilter(trackFilter);
  mgr->AddTask(esdFilter);
  
  // second task - ESD filter task putting muon info in the output generic AOD 
  AliAnalysisTaskESDMuonFilter *esdMuonFilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
  mgr->AddTask(esdMuonFilter);
  
  // third task - compare created AOD and exixting ESDs
  AliAnalysisTaskAODvsESD *AODvsESD = new AliAnalysisTaskAODvsESD("aodVsEsd");
  mgr->AddTask(AODvsESD);

  // Input ESD container
  AliAnalysisDataContainer *esdIn = mgr->GetCommonInputContainer();
  // Output AOD container. 
  AliAnalysisDataContainer *aodOut = mgr->GetCommonOutputContainer();
  // Output comparison
  AliAnalysisDataContainer *listOut = mgr->CreateContainer("output2", TList::Class(), AliAnalysisManager::kOutputContainer, "AODvsESDoutput.root");

  // Connect containers to tasks slots
  mgr->ConnectInput(esdFilter,0,esdIn); 
  mgr->ConnectOutput(esdFilter,0,aodOut);
  
  mgr->ConnectInput(esdMuonFilter,0,esdIn);
  mgr->ConnectOutput(esdMuonFilter,0,aodOut);

  mgr->ConnectInput(AODvsESD,0,esdIn);
  mgr->ConnectOutput(AODvsESD,0,listOut);
  
  // Run the analysis
  if (mgr->InitAnalysis()){
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }   
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
      printf("*** Building PAR archive    ***\n");
      
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	Error("runProcess","Cannot Build the PAR Archive! - Abort!");
	return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*** Setup PAR archive       ***\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory(ocwd.Data());
    printf("Current dir: %s\n", ocwd.Data());
  } 
}
