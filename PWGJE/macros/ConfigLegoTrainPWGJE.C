void ConfigLegoTrainPWGJE(int iFlag = 0){

  Printf("%s:%d with flag %d",(char*)__FILE__,__LINE__,iFlag);
  // Adds the Global Variables depending on the flag
  // 

  if(iFlag==1008){ // 10h
    AliAnalysisManager::SetGlobalStr("kJetRunPeriod","LHC10h");
    AliAnalysisManager::SetGlobalInt("kGridRunRangeLo",110000);
    AliAnalysisManager::SetGlobalInt("kGridRunRangeUp",160000);

    AliAnalysisManager::SetGlobalStr("kDeltaAODJetName","");

    // event selection
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


    AliAnalysisManager::SetGlobalInt("kPhysicsSelectionFlag",AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kEMCEJE);    
    AliAnalysisManager::SetGlobalInt("kNTrigger",5);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit0",AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kEMCEJE);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit1",AliVEvent::kMB);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit2",AliVEvent::kCentral);    
    AliAnalysisManager::SetGlobalInt("kTriggerBit3",AliVEvent::kSemiCentral);
    AliAnalysisManager::SetGlobalInt("kTriggerBit4",AliVEvent::kEMCEJE);
    AliAnalysisManager::SetGlobalStr("kTriggerName0","kMB|kCentral|kSemiCentral|kEMCEJE");    
    AliAnalysisManager::SetGlobalStr("kTriggerName1","kMB");    
    AliAnalysisManager::SetGlobalStr("kTriggerName2","kCentral");    
    AliAnalysisManager::SetGlobalStr("kTriggerName3","kSemiCentral");    
    AliAnalysisManager::SetGlobalStr("kTriggerName4","kEMCEJE");    


    // aceptance windows for Spectrum task

    Float_t fid = 0.2;
    Float_t fEMCPhiLo = 2.17 - 1.87/2. + fid; // some fiducial region
    Float_t fEMCPhiUp = 2.17 + 1.87/2. - fid; // some fiducial region
    Float_t fEMCEtaLo = - 0.7 + fid;
    Float_t fEMCEtaUp =   0.7 - fid;
    AliAnalysisManager::SetGlobalInt("kNAcceptanceSpec",2);    
    AliAnalysisManager::SetGlobalDbl("kAcceptancePhiMinSpec0",fEMCPhiLo);
    AliAnalysisManager::SetGlobalDbl("kAcceptancePhiMaxSpec0",fEMCPhiUp);
    AliAnalysisManager::SetGlobalDbl("kAcceptanceEtaMinSpec0",fEMCEtaLo);
    AliAnalysisManager::SetGlobalDbl("kAcceptanceEtaMaxSpec0",fEMCEtaUp);

    // iroc 13 
    AliAnalysisManager::SetGlobalDbl("kAcceptancePhiMinSpec1",4.7-0.4);
    AliAnalysisManager::SetGlobalDbl("kAcceptancePhiMaxSpec1",4.7+0.4);
    AliAnalysisManager::SetGlobalDbl("kAcceptanceEtaMinSpec1",-0.9);
    AliAnalysisManager::SetGlobalDbl("kAcceptanceEtaMaxSpec1",0);


    // jet selection
    AliAnalysisManager::SetGlobalDbl("kJetEtaWindow",0.5);

    // track selection 
    AliAnalysisManager::SetGlobalDbl("kTrackEtaWindow",0.9);
    AliAnalysisManager::SetGlobalDbl("kVertexWindow",10);
    AliAnalysisManager::SetGlobalInt("kHighPtFilterMask",768);
    AliAnalysisManager::SetGlobalInt("kHighPtFilterMaskBest",256);

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
