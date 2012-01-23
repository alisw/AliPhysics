// 19 Nov 2007
// Macro for the running of the AliAnalysisTaskMuonAODfromGeneral
//

void RunMuonAODfromGeneral(char* filein = "../FromESDToGenAOD/AliAOD.root", char* fileout = "AliMuonAOD.root" ){
     
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");			
    gSystem->Load("libPWG3base.so");
             
    // Input AOD file
    TChain* chain = new TChain("aodTree");
    chain->Add(filein);

    // Make aod output handler
    AliAODHandler* aodHandler = new AliAODHandler();
    aodHandler->SetOutputFileName(fileout);
    
    // Make the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("AOD Manager", "AOD Manager");
    mgr->SetOutputEventHandler(aodHandler);
    mgr-> SetDebugLevel(10);
    
    // Task for AOD generation from PWG0base directory
    AliAnalysisTaskMuonAODfromGeneral *aodfilter = new AliAnalysisTaskMuonAODfromGeneral("AOD Filter",7000.);
    mgr->AddTask(aodfilter);
  
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histo", TH1F::Class(),
   							      AliAnalysisManager::kOutputContainer, "default");

    mgr->ConnectInput  (aodfilter,  0, cinput1  );
    mgr->ConnectOutput (aodfilter,  0, coutput1 );
    
    // Run the analysis    
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
