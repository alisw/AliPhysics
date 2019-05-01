void JetAnalysisManager()
{
    //
    // Load relevant libraries
    //
    gSystem->Load("libTree");
    gSystem->Load("libNetx");
    gSystem->Load("libProof");
    gSystem->Load("libProofPlayer");
    gSystem->Load("libGeom");
    gSystem->Load("libEG");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libESD");
    gSystem->Load("libJETAN");
    //
    // Connect to alien
    //
    TGrid::Connect("alien://"); 
    
    AliTagAnalysis *TagAna = new AliTagAnalysis(); 
    AliEventTagCuts *EvCuts = new AliEventTagCuts();
    AliRunTagCuts   *RuCuts = new AliRunTagCuts();
    //EvCuts->SetNChargedAbove1GeVRange(1, 1000);
    //EvCuts->SetMultiplicityRange(11,120);
    //EvCuts->SetNPionRange(2,10000);
     TGridCollection* coll = gGrid->OpenCollection("tag100.xml");
     TGridResult* TagResult = coll->GetGridResult("", 0, 0);
     TagResult->Print();
     TagAna->ChainGridTags(TagResult);

  //////////////////////////////////////////////////////////////////
  //Get the chain
     printf("*******************************\n");
     printf("*** Getting the Chain       ***\n");
     printf("*******************************\n");
     TChain* chain1 = 0x0;
     chain1 = TagAna->QueryTags(RuCuts, EvCuts);
     chain1->ls();
     
    //
    // Make the analysis manager
    //
    AliAnalysisManager *mgr = new AliAnalysisManager("Manager", "Manager");
    mgr-> SetDebugLevel(10);
    
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    
    mgr->AddTask(jetana);
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),
							      AliAnalysisManager::kOutputContainer);

    mgr->ConnectInput (jetana,0,cinput1);
    mgr->ConnectOutput(jetana,0,coutput1);
    cinput1->SetData(chain1);

//
// Run the analysis
//    

    if (mgr->InitAnalysis()) {
	mgr->PrintStatus();
	mgr->StartAnalysis("local", chain1);
    }
}
