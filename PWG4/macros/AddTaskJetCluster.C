AliAnalysisTaskJetCluster *AddTaskJetCluster(char* bRec = "AOD",char* bGen = "",UInt_t filterMask = 16, UInt_t iPhysicsSelectionFlag = AliVEvent::kMB,Char_t *jf = "KT", Float_t radius = 0.4,Int_t nSkip = 0,Int_t kWriteAOD = kFALSE,char* deltaFile = "",Float_t ptTrackCut = 0.15, Float_t etaTrackWindow = 0.9,Int_t nSkipCone = 2);

Int_t kBackgroundMode = 0;
Float_t kPtTrackCut = 0.15;
Float_t kTrackEtaWindow = 0.8;

AliAnalysisTaskJetCluster *AddTaskJetClusterDelta(UInt_t filterMask = 16,Bool_t kUseAODMC = kFALSE,UInt_t iPhysicsSelectionFlag = AliVEvent::kMB,Char_t *jf = "KT", UInt_t iFlag){
   AliAnalysisTaskJetCluster *js = 0;
   if(kUseAODMC&&false){// do not use the MC info yet
     if(iFlag&(1<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,0.00001);
     if(iFlag&(1<<1))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,0.1);
     if(iFlag&(1<<2))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,0.2);
     if(iFlag&(1<<4))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,0.4);
     if(iFlag&(1<<6))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,0.6);
     if(iFlag&(1<<8))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,0.8);
     if(iFlag&(1<<10))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelectionFlag,jf,1.0);
   }
   else{
     if(iFlag&(1<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,0.00001);
     if(iFlag&(1<<1))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,0.1);
     if(iFlag&(1<<2))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,0.2);
     if(iFlag&(1<<4))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,0.4);
     if(iFlag&(1<<6))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,0.6);
     if(iFlag&(1<<8))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,0.8);
     if(iFlag&(1<<10))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelectionFlag,jf,1.0);
   }

   return js;
 }


AliAnalysisTaskJetCluster *AddTaskJetCluster(char* bRec,char* bGen ,UInt_t filterMask,UInt_t iPhysicsSelectionFlag,Char_t *jf,Float_t radius,Int_t nSkip,Int_t kWriteAOD,char *deltaFile,Float_t ptTrackCut,Float_t etaTrackWindow,Int_t nSkipCone)
 {
 // Creates a jet fider task, configures it and adds it to the analysis manager.
   kPtTrackCut = ptTrackCut;
   kTrackEtaWindow = etaTrackWindow;

   TString outputFile(deltaFile);
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
       ::Error("AddTaskJetCluster", "No analysis manager to connect to.");
       return NULL;
    }  

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetCluster", "This task requires an input event handler");
       return NULL;
    }

    TString type = mgr->GetInputEventHandler()->GetDataType();
    TString typeRec(bRec);
    TString typeGen(bGen);
    if(!typeRec.Contains("AODextra")) {
      typeGen.ToUpper();
      typeRec.ToUpper();
    }
    cout << "typeRec: " << typeRec << endl;
    // Create the task and configure it.
    //===========================================================================




    TString cAdd = "";
    cAdd += Form("%02d_",(int)((radius+0.01)*10.));
    cAdd += Form("B%d",(int)kBackgroundMode);
    cAdd += Form("_Filter%05d",filterMask);
    cAdd += Form("_Cut%05d",(int)(1000.*kPtTrackCut));
    cAdd += Form("_Skip%02d",nSkip);
    Printf("%s %E",cAdd.Data(),kPtTrackCut);
    AliAnalysisTaskJetCluster* pwg4spec = new  AliAnalysisTaskJetCluster(Form("JetCluster%s_%s%s",bRec,jf,cAdd.Data()));
      
   // or a config file
   // pwg4spec->SetAnalysisType(AliAnalysisTaskJetCluster::kAnaMC);
   // if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   // pwg4spec->SetDebugLevel(11); 
   pwg4spec->SetFilterMask(filterMask); 
   //   pwg4spec->SetUseGlobalSelection(kTRUE); 

   if(type == "AOD"){
     // Assume all jet are produced already
     pwg4spec->SetAODTrackInput(kTRUE);
     pwg4spec->SetAODMCInput(kTRUE);
   }

   if(typeRec.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCChargedAcceptance);
     pwg4spec->SetTrackPtCut(kPtTrackCut);
     pwg4spec->SetTrackEtaWindow(kTrackEtaWindow);
   }
   else if (typeRec.Contains("AODMC2")){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCCharged);
     pwg4spec->SetTrackPtCut(kPtTrackCut);
     pwg4spec->SetTrackEtaWindow(5);
   }
   else if (typeRec.Contains("AODMC")){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCAll);
     pwg4spec->SetTrackPtCut(kPtTrackCut);
     pwg4spec->SetTrackEtaWindow(5);
   }
   else if (typeRec.Contains("AODextraonly")) {
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODextraonly);
     pwg4spec->SetTrackPtCut(kPtTrackCut);
     pwg4spec->SetTrackEtaWindow(kTrackEtaWindow);
   }
   else if (typeRec.Contains("AODextra")) {
     cout << "AliAnalysisTaskJetCluster::kTrackAODextra: " << AliAnalysisTaskJetCluster::kTrackAODextra << endl;
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODextra);
     pwg4spec->SetTrackPtCut(kPtTrackCut);
     pwg4spec->SetTrackEtaWindow(kTrackEtaWindow);
   }
   else if (typeRec.Contains("AOD")) {
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAOD);
     pwg4spec->SetTrackPtCut(kPtTrackCut);
     pwg4spec->SetTrackEtaWindow(kTrackEtaWindow);
   }

   pwg4spec->SetRparam(radius);

   switch (jf) {
   case "ANTIKT":
     pwg4spec->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
     break;
   case "CA":
     pwg4spec->SetAlgorithm(1); // CA from fastjet/JetDefinition.hh
     break;
   case "KT":
     pwg4spec->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
     break;
   default:
     ::Error("AddTaskJetCluster", "Wrong jet finder selected\n");
     return 0;
   }

   
   if(kWriteAOD){
     if(outputFile.Length())pwg4spec->SetJetOutputFile(outputFile);
     pwg4spec->SetJetOutputBranch(Form("clusters%s_%s%s",bRec,jf,cAdd.Data()));
     pwg4spec->SetJetOutputMinPt(0); // store only jets / clusters above a certain threshold
   }

   pwg4spec->SetNSkipLeadingRan(nSkip);
   pwg4spec->SetNSkipLeadingCone(nSkipCone);

   if(iPhysicsSelectionFlag)pwg4spec->SelectCollisionCandidates(iPhysicsSelectionFlag);

   mgr->AddTask(pwg4spec);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer(Form("pwg4cluster_%s_%s_%s%s",bRec,bGen,jf,cAdd.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_cluster_%s_%s_%s%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen,jf,cAdd.Data()));

   mgr->ConnectInput  (pwg4spec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4spec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );
   
   return pwg4spec;
}
