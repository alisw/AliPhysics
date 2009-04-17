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

     //
    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //

    //TString path("/afs/cern.ch/user/m/morsch/public/");
    TString path("./");
    TChain* chain = new TChain("aodTree");
    chain->Add(Form("%s/%s",path.Data(),"AliAOD.root"));

    /////////////////////////////////////////////////////////////////////////////////// 
    // Create the analysis manager
    //
    // Input 
    AliMultiEventInputHandler* inpHandler = new AliMultiEventInputHandler(2, 1);
    // Pool
    AliEventPoolOTF* pool = new AliEventPoolOTF("event pool", "AOD");

    pool->SetTagDirectory(path.Data());
    pool->SetMultiplicityBin(0, 100, 100);
    pool->Init();
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetInputEventHandler  (inpHandler);
    mgr->SetEventPool(pool);
    inpHandler->SetEventPool(pool);
    

    mgr->SetDebugLevel(10);
    /////////////////////////////////////////////////////////////////////////////////// 
    gROOT->LoadMacro("AliAnalysisTaskPhiCorr.cxx++g");
    AliAnalysisTaskPhiCorr *dphiana = new AliAnalysisTaskPhiCorr("Phi Correlation Analysis");
    dphiana->SetDebugLevel(10);
    mgr->AddTask(dphiana);
    
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kExchangeContainer, "default");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");


    mgr->ConnectInput  (dphiana,  0,  mgr->GetCommonInputContainer());
    mgr->ConnectOutput (dphiana,  0, coutput1 );
    mgr->ConnectOutput (dphiana,  1, coutput2 );

    // 
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("mix",chain, 1000);
}
