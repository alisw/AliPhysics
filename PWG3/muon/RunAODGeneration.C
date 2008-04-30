//
//  Macro for AOD generation 
//  Gines Martinez, Subatech, October. 2007
//  Generated from Andrea Morsch macro JetAnalysisManagerLoc.C
//
//  In this example the libraries PWG0base and PWG3base (not really needed in this version) 
//  are supposed to be loaded from a par file via RunAnalysis
//
void RunAODGeneration(char* filein = "AliESDs.root", char* fileout = "AliAOD.root" )
{
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libAOD");
    gSystem->Load("libESD");  
  
    // Input ESD files
    TChain* chain = new TChain("esdTree");
    chain->Add(filein);
    
    // Make aod handler
    AliAODHandler* aodHandler = new AliAODHandler();
    aodHandler->SetOutputFileName(fileout);
    
    // Make the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("AOD Manager", "AOD Manager");
    mgr->SetOutputEventHandler(aodHandler);
    mgr-> SetDebugLevel(10);
    
    // Task for AOD generation from PWG0base directory
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);
  
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");

    mgr->ConnectInput  (esdfilter,  0, cinput1  );
    mgr->ConnectOutput (esdfilter,  0, coutput1 );

    
    // Run the analysis    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
}
