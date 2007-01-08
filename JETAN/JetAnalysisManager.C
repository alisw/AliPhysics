void JetAnalysisManager()
{
    //
    // Load relevant libraries
    //
    gSystem->Load("libEG.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libANALYSIS_NEW.so");
    gSystem->Load("libESD.so"); 
    gSystem->Load("libJETAN.so");
    //
    // Connect to alien
    //
    TGrid::Connect("alien://"); 
    //
    // Prepare the input event chain
    //
    AliTagAnalysis *TagAna = new AliTagAnalysis(); 
    // create an EventTagCut object
    AliEventTagCuts *EvCuts = new AliEventTagCuts();
    AliRunTagCuts   *RuCuts = new AliRunTagCuts();
    TAlienCollection* coll = TAlienCollection::Open("tags/tag100.xml");
    TGridResult* TagResult = coll->GetGridResult("");
    TagAna->ChainGridTags(TagResult);
    TChain* chain1 = 0x0;
    chain1 = TagAna->QueryTags(RuCuts, EvCuts);
    //
    // Make the analysis manager
    //
    AliAnalysisManager *mgr = new AliAnalysisManager();
    AliAnalysisTask *jetana = new AliAnalysisTaskJets("JetAnalysis");
    
    mgr->AddTask(jetana);
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),
							      AliAnalysisManager::kOutputContainer);
    printf("Connectiong I/O \n");
    
    mgr->ConnectInput (jetana,0,cinput1);
    mgr->ConnectOutput(jetana,0,coutput1);
    cinput1->SetData(chain1);

//
// Run the analysis
//    

    if (mgr->InitAnalysis()) {
	mgr->PrintStatus();
	chain1->Process(mgr);
    }
}
