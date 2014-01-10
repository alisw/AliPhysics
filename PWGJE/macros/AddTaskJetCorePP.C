

AliAnalysisTaskJetCorePP* AddTaskJetCorePP(
   const Char_t* branchPrefix="clustersAOD",
   const Char_t* jetAlgo="ANTIKT",   
   Float_t jetParameterR = 0.4,  //jet R 
   UInt_t  trkFilterMask = 272, 
   Float_t trackLowPtCut = 0.15,
   const Char_t* jetbgAlgo="ANTIKT",    //background jet algo
   Float_t bgjetParameterR = 0.3,  //R of jet to be removed while bg calc 
   Float_t bgMaxJetPt = 8.0, //max jet pt to be accepted to bg 
   Int_t  rndTrials = 2000, //number of trials to get jet free cell area
   Float_t jetFreeAreaFrac = 0.75, //cell area free of jets
   Float_t bgConeR = 0.4,  //R of perp cone jet R 
   Int_t   collisionSystem = 0, //pp=0, pPb=1
   Int_t   offlineTriggerMask=AliVEvent::kMB, //MinBias=0 
   Int_t   minContribVtx = 1,
   Float_t vtxZMin = -10.0,
   Float_t vtxZMax = 10.0,
   Float_t centMin = 0.0,
   Float_t centMax = 100.0,
   Float_t triggerEtaCut = 0.9,
   Float_t trackEtaCut = 0.9,
   const Char_t* nonStdFile="",
   const Char_t* mcFullFlag="",  // real="", all jets= "MC"    
   const Char_t* mcChargFlag="",  // real="", charged jets = "MC2" 
   Bool_t bfillrespmx=0,  // 0=dont fill resp mx histos, 1=fill histos
   Bool_t bDiceEff=0,  // 0=leave efficiency as it is,  1= reduce efficiency by constant amount via SetFixedEfficiency
   Int_t triggerType=0,  //0=single incl trigger, 1=leading track, 2=hadron pt>10 
   Int_t evtRangeLow=0,   //last digit of range of ESD event number
   Int_t evtRangeHigh=9,  //first digit of range of ESD event number
   Float_t trigRangeLow=0,  //trigger pT low bin boreder works with triggType=0
   Float_t trigRangeHigh=50  //trigger pT high border works with triggType=0
  ){ 

   Printf("adding task jet response\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskJetCorePP", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskJetCorePP", "This task requires an input event handler.");
      return NULL;
   }

   Float_t jetEtaMin = -0.9 + jetParameterR;
   Float_t jetEtaMax =  0.9 - jetParameterR; 
    
   TString analBranch(branchPrefix);
   TString analBranchBg(branchPrefix);  // kT bg jets
   TString analBranchFullMC="";             //full jets MC to be used with charged jets MC 
   TString analBranchChargMC="";            //charged jets
   TString analBranchBgChargMC="";          //charged jets kt background 


   TString stJetAlgo(jetAlgo);
   TString stJetBgAlgo(jetbgAlgo);
   stJetAlgo.ToUpper();
   stJetBgAlgo.ToUpper();

 
   TString jet="";
   TString jetbg="";
   TString otherparams="";
   jet   = jet   + "_" + stJetAlgo   + Form("%02d",(Int_t) (10*jetParameterR));
   jetbg = jetbg + "_" + stJetBgAlgo + Form("%02d",(Int_t) (10*bgjetParameterR));
   
   otherparams = otherparams + "_B0"; //bg mode
   otherparams = otherparams + Form("_Filter%05d",(UInt_t) trkFilterMask);
   otherparams = otherparams + Form("_Cut%05d",(Int_t) (1000*trackLowPtCut));

   if(analBranch.BeginsWith("clustersAOD")){
      otherparams = otherparams + Form("_Skip%02d",0);
   }
   
   analBranch   = analBranch   + jet   + otherparams; //antikt jet 
   analBranchBg = analBranchBg + jetbg + otherparams; //kt bg jet

   if(bDiceEff){ //dicing efficiency relates rec only
      analBranch   = analBranch   + "Detector10Fr0"; //dice=1, smear=0, change eff fraction =0
      analBranchBg = analBranchBg + "Detector10Fr0"; 
   }

   //clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00   
   //clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_Detector10Fr0   
   //Skip00 none of the most energetic jets is ommited
   //Cut00150  pT min cut on track
   //Filter00272

   TString mcFullSuffix(mcFullFlag);    //MC = all jets
   TString mcChargSuffix(mcChargFlag);  //MC2= charged jets,  MC = all jets

   if(mcChargSuffix.Length()>0 && mcChargSuffix=="MC2"){  //charged jets generator level
      analBranchChargMC   = branchPrefix + mcChargSuffix + jet   + otherparams; 
      analBranchBgChargMC = branchPrefix + mcChargSuffix + jetbg + otherparams; 
   }

   if(mcFullSuffix.Length()>0 && mcFullSuffix=="MC"){ //full jets generator level
      analBranchFullMC    = branchPrefix + mcFullSuffix  + jet   + otherparams; 
   }


 
   AliAnalysisTaskJetCorePP *task = new AliAnalysisTaskJetCorePP(Form("JetCorePP_%s_%s_%d",analBranch.Data(), mcChargSuffix.Data(), offlineTriggerMask));

   task->SetBranchName(analBranch.Data());
   task->SetBranchNameFullMC(analBranchFullMC.Data());
   task->SetBranchNameChargMC(analBranchChargMC.Data());
   task->SetBranchNameBg(analBranchBg.Data()); //kt bg jets
   task->SetBranchNameBgChargMC(analBranchBgChargMC.Data()); //kt bg jets
   task->SetNonStdFile(nonStdFile);
   task->SetSystem(collisionSystem); 
   task->SetJetR(jetParameterR);
   task->SetBgJetR(bgjetParameterR);
   task->SetBgMaxJetPt(bgMaxJetPt);
   task->SetRndTrials(rndTrials);
   task->SetFreeAreaFrac(jetFreeAreaFrac);
   task->SetBgConeR(bgConeR); 
   task->SetOfflineTrgMask(offlineTriggerMask);
   task->SetMinContribVtx(minContribVtx);
   task->SetVtxZMin(vtxZMin);
   task->SetVtxZMax(vtxZMax);
   task->SetFilterMask(trkFilterMask);
   task->SetCentMin(centMin);
   task->SetCentMax(centMax);
   task->SetJetEtaMin(jetEtaMin);
   task->SetJetEtaMax(jetEtaMax);
   task->SetTriggerEtaCut(triggerEtaCut);
   task->SetTrackEtaCut(trackEtaCut);
   task->SetTrackLowPtCut(trackLowPtCut);
   task->SetTriggerType(triggerType); 
   task->SetEventNumberRangeLow(evtRangeLow);
   task->SetEventNumberRangeHigh(evtRangeHigh);
   task->SetTriggerPtRangeLow(trigRangeLow);
   task->SetTriggerPtRangeHigh(trigRangeHigh); 
   task->SetFillResponseMatrix(bfillrespmx);

   task->SetDebugLevel(0); //No debug messages 0
   mgr->AddTask(task);
   //E=  range of last two decimal numbers in event numbers
   //Ptt range of the cosidered trigger bin
   AliAnalysisDataContainer *coutputJetCorePP = mgr->CreateContainer(
      Form("pwgjejetcorepp_%s_%s%02d_%s_%d_T%d_E%d_%d_Ptt%.0f_%.0f",analBranch.Data(),jetbgAlgo,TMath::Nint(10*bgjetParameterR),mcChargSuffix.Data(),offlineTriggerMask,triggerType,evtRangeLow,evtRangeHigh,trigRangeLow,trigRangeHigh), 
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:PWGJE_jetcorepp_%s_%s%02d_%s_%d_T%d_E%d_%d_Ptt%.0f_%.0f",AliAnalysisManager::GetCommonFileName(),analBranch.Data(),jetbgAlgo,TMath::Nint(10*bgjetParameterR),mcChargSuffix.Data(),offlineTriggerMask,triggerType,evtRangeLow,evtRangeHigh,trigRangeLow,trigRangeHigh)
   );

   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetCorePP);

   return task;
}
