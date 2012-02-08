void ConfigLegoTrainPWGJE(int iFlag = 0){

  // Adds the Global Variables depending on the flag
  // 

  if(iFlag==1008){ // 10h
    AliAnalysisManager::SetGlobalStr("kJetRunPeriod","LHC10h");
    AliAnalysisManager::SetGlobalInt("kGridRunRangeLo",110000);
    AliAnalysisManager::SetGlobalInt("kGridRunRangeUp",160000);

    AliAnalysisManager::SetGlobalStr("kDeltaAODJetName","");

    AliAnalysisManager::SetGlobalInt("kPhysicsSelectionFlag",AliVEvent::kMB);
    AliAnalysisManager::SetGlobalInt("kNTrigger",1);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit0",AliVEvent::kMB);    
    AliAnalysisManager::SetGlobalStr("kTriggerName0","kMB");    

    // jet selection
    AliAnalysisManager::SetGlobalDbl("kJetEtaWindow",0.5);

    // track selection 
    AliAnalysisManager::SetGlobalDbl("kTrackEtaWindow",0.9);
    AliAnalysisManager::SetGlobalDbl("kVertexWindow",10);
    AliAnalysisManager::SetGlobalInt("kHighPtFilterMask",272);


  }
  else if(iFlag==1108){ // 11h


    AliAnalysisManager::SetGlobalStr("kJetRunPeriod","LHC11h");
    AliAnalysisManager::SetGlobalInt("kGridRunRangeLo",166000);
    AliAnalysisManager::SetGlobalInt("kGridRunRangeUp",171000);


    AliAnalysisManager::SetGlobalStr("kDeltaAODJetName","");


    AliAnalysisManager::SetGlobalInt("kPhysicsSelectionFlag",AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral);    
    AliAnalysisManager::SetGlobalInt("kNTrigger",1);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit0",AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit1",AliVEvent::kMB);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit2",AliVEvent::kCentral);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit3",AliVEvent::kSemiCentral);
    AliAnalysisManager::SetGlobalStr("kTriggerName0","kMB|kCentral|kSemiCentral");    
    AliAnalysisManager::SetGlobalStr("kTriggerName1","kMB");    
    AliAnalysisManager::SetGlobalStr("kTriggerName2","kCentral");    
    AliAnalysisManager::SetGlobalStr("kTriggerName3","kSemiCentral");    


    // jet selection
    AliAnalysisManager::SetGlobalDbl("kJetEtaWindow",0.5);

    // track selection 
    AliAnalysisManager::SetGlobalDbl("kTrackEtaWindow",0.9);
    AliAnalysisManager::SetGlobalDbl("kVertexWindow",10);
    AliAnalysisManager::SetGlobalInt("kHighPtFilterMask",272);

  }
  else{
    // catch all
    AliAnalysisManager::SetGlobalInt("kGridRunRangeLo",110000);
    AliAnalysisManager::SetGlobalInt("kGridRunRangeUp",160000);


    AliAnalysisManager::SetGlobalStr("kDeltaAODJetName","");

    AliAnalysisManager::SetGlobalInt("kPhysicsSelectionFlag",AliVEvent::kMB);
    AliAnalysisManager::SetGlobalInt("kNTrigger",1);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit0",AliVEvent::kMB);    
    AliAnalysisManager::SetGlobalStr("kTriggerName0","kMB");    


    // jet selection
    AliAnalysisManager::SetGlobalDbl("kJetEtaWindow",0.5);

    // track selection 
    AliAnalysisManager::SetGlobalDbl("kTrackEtaWindow",0.9);
    AliAnalysisManager::SetGlobalDbl("kVertexWindow",10);
    AliAnalysisManager::SetGlobalInt("kHighPtFilterMask",272);



  }


}
