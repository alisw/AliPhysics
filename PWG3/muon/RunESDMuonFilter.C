// Macro to run AliAnalysisTaskESDMuonFilter
//

void RunESDMuonFilter(char* filein = "AliESDs.root", char* fileout = "AliMuonAOD.root" ){
     
    gSystem->Load("libTree.so"); 
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics");
    
    // for analysis .par file based
    
    setupPar("STEERBase");
    setupPar("ESD");
    setupPar("AOD");
    setupPar("ANALYSIS");
    setupPar("ANALYSISalice");
    setupPar("PWG3base");
    
    // Input ESD file
    TChain* chain = new TChain("esdTree");  
    chain->Add(filein);
    AliESDInputHandler* esdHandler = new AliESDInputHandler();

    // Make aod output handler
    AliAODHandler* aodHandler = new AliAODHandler();
    aodHandler->SetOutputFileName(fileout);
    
    // Make the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("Muon AOD Manager", "Muon AOD Manager");
    mgr->SetInputEventHandler(esdHandler);
    mgr->SetOutputEventHandler(aodHandler);
    mgr-> SetDebugLevel(10);
    
    // Task for MUON AOD generation 
    AliAnalysisTaskESDMuonFilter *esdfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
    mgr->AddTask(esdfilter);
  
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cESD",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cAOD", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");

    mgr->ConnectInput  (esdfilter,  0, cinput1  );
    mgr->ConnectOutput (esdfilter,  0, coutput1 );
    
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
