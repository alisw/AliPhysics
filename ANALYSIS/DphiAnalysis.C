void DphiAnalysis()
{
      gSystem->Load("libTree.so");
      gSystem->Load("libPhysics.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libSTEERBase.so");
      gSystem->Load("libAOD.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libANALYSISalice.so");
      gSystem->Load("libJETAN.so");

     //
    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //
    TChain* chain = new TChain("aodTree");
    chain->Add("./AliAODs.root");

    /////////////////////////////////////////////////////////////////////////////////// 
    // Create the analysis manager
    //
    // Input 
    AliMultiAODInputHandler* inpHandler = new AliMultiAODInputHandler(2);
    // Pool
    AliEventPoolOTF* pool = new AliEventPoolOTF("event pool", "event pool");
    pool->SetTagDirectory(".");
    pool->SetMultiplicityBin(0, 100, 1);
    pool->Init();
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetInputEventHandler  (inpHandler);
    mgr->SetEventPool(pool);
    inpHandler->SetEventPool(pool);
    

    mgr->SetDebugLevel(10);
    /////////////////////////////////////////////////////////////////////////////////// 
    AliAnalysisTaskPhiCorr *dphiana = new AliAnalysisTaskPhiCorr("Phi Correlation Analysis");
    dphiana->SetDebugLevel(10);
    mgr->AddTask(dphiana);
    
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

//    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
//							      AliAnalysisManager::kOutputContainer, "default");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");


    mgr->ConnectInput  (dphiana,  0, cinput1  );
    mgr->ConnectOutput (dphiana,  1, coutput2 );

    // 
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("mix",chain, 1000);
}
