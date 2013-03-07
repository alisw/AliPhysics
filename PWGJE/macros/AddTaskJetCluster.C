AliAnalysisTaskJetCluster *AddTaskJetCluster(char* bRec = "AOD",char* bGen = "",UInt_t filterMask = 16, UInt_t iPhysicsSelectionFlag = AliVEvent::kAny,Char_t *jf = "KT", Float_t radius = 0.4,Int_t nSkip = 0,Int_t kWriteAOD = kFALSE,char* deltaFile = "",Float_t ptTrackCut = 0.15, Float_t etaTrackWindow = 0.9,Float_t vertexWindow = 10.,Int_t nSkipCone = 2,Int_t dice=0,Int_t smear=0,Bool_t useOADB=kFALSE,Double_t changeEfficiencyFraction=0.);
AliAnalysisTaskJetCluster *AddTaskJetClusterDelta(UInt_t filterMask = 16,Bool_t kUseAODMC = kFALSE,UInt_t iPhysicsSelectionFlag = AliVEvent::kMB,Char_t *jf = "KT", UInt_t iFlag);

Int_t kBackgroundModeCl = 0;
Float_t kPtTrackCutCl = 0.15;
Float_t kTrackEtaWindowCl = 0.8;
Float_t kVertexWindowCl = 10;

AliAnalysisTaskJetCluster *AddTaskJetClusterDelta(UInt_t filterMask,Bool_t kUseAODMC,UInt_t iPhysicsSelectionFlag,Char_t *jf, UInt_t iFlag){
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


AliAnalysisTaskJetCluster *AddTaskJetCluster(char* bRec,char* bGen ,UInt_t filterMask,UInt_t iPhysicsSelectionFlag,Char_t *jf,Float_t radius,Int_t nSkip,Int_t kWriteAOD,char *deltaFile,Float_t ptTrackCut,Float_t etaTrackWindow,Float_t vertexWindow,Int_t nSkipCone,Int_t dice,Int_t smear,Bool_t useOADB,Double_t changeEfficiencyFraction)
 {
 // Creates a jet fider task, configures it and adds it to the analysis manager.
   kPtTrackCutCl = ptTrackCut;
   kTrackEtaWindowCl = etaTrackWindow;
   kVertexWindowCl = vertexWindow;

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
    if(!typeRec.Contains("AODextra") || !typeRec.Contains("AODMCextra")) {
      typeGen.ToUpper();
      typeRec.ToUpper();
    }
    cout << "typeRec: " << typeRec << endl;
    // Create the task and configure it.
    //===========================================================================




    TString cAdd = "";
    cAdd += Form("%02d_",(int)((radius+0.01)*10.));
    cAdd += Form("B%d",(int)kBackgroundModeCl);
    cAdd += Form("_Filter%05d",filterMask);
    cAdd += Form("_Cut%05d",(int)(1000.*kPtTrackCutCl));
    cAdd += Form("_Skip%02d",nSkip);
    if(dice>0 || smear>0)
      cAdd += Form("Detector%d%dFr%d",dice,smear,(int)(changeEfficiencyFraction*100.));
    

    Printf("%s %E %d %d",cAdd.Data(),kPtTrackCutCl,dice,smear);
    AliAnalysisTaskJetCluster* clus = new  AliAnalysisTaskJetCluster(Form("JetCluster%s_%s%s",bRec,jf,cAdd.Data()));
      
   // or a config file
   // clus->SetAnalysisType(AliAnalysisTaskJetCluster::kAnaMC);
   // if(iAODanalysis)clus->SetAODInput(kTRUE);
   // clus->SetDebugLevel(11); 
   clus->SetFilterMask(filterMask); 
   //   clus->SetUseGlobalSelection(kTRUE); 
   clus->SetVtxCuts(kVertexWindowCl,1);
   if(type == "AOD"){
     // Assume all jet are produced already
     clus->SetAODTrackInput(kTRUE);
     clus->SetAODMCInput(kTRUE);
   }

   if(typeRec.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCChargedAcceptance);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AODMC2")){
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCCharged);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(5);
   }
   else if (typeRec.Contains("AODMCextraonly")) {
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCextraonly);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AODMCextra")) {
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCextra);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AODMC")){
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCAll);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(5);
   }
   else if (typeRec.Contains("AODextraonly")) {
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODextraonly);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AODextra")) {
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODextra);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AOD")) {
     clus->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAOD);
     clus->SetTrackPtCut(kPtTrackCutCl);
     clus->SetTrackEtaWindow(kTrackEtaWindowCl);
   }

   clus->SetRparam(radius);
   clus->SetGhostArea(0.005);
   clus->SetGhostEtamax(kTrackEtaWindowCl);

   switch (jf) {
   case "ANTIKT":
     clus->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
     break;
   case "CA":
     clus->SetAlgorithm(1); // CA from fastjet/JetDefinition.hh
     break;
   case "KT":
     clus->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
     break;
   default:
     ::Error("AddTaskJetCluster", "Wrong jet finder selected\n");
     return 0;
   }

   
   if(kWriteAOD){
     if(outputFile.Length())clus->SetJetOutputFile(outputFile);
     clus->SetJetOutputBranch(Form("clusters%s_%s%s",bRec,jf,cAdd.Data()));
     clus->SetJetOutputMinPt(0); // store only jets / clusters above a certain threshold
   }

   clus->SetNSkipLeadingRan(nSkip);
   clus->SetNSkipLeadingCone(nSkipCone);

   if(iPhysicsSelectionFlag)clus->SelectCollisionCandidates(iPhysicsSelectionFlag);

   if(useOADB) {
     clus->SetUseTrResolutionFromOADB();
     clus->SetUseTrEfficiencyFromOADB();
     clus->SetChangeEfficiencyFraction(changeEfficiencyFraction);
   }

   mgr->AddTask(clus);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_clus = mgr->CreateContainer(Form("cluster_%s_%s_%s%s",bRec,bGen,jf,cAdd.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_cluster_%s_%s_%s%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen,jf,cAdd.Data()));

   mgr->ConnectInput  (clus, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (clus, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (clus,  1, coutput1_clus );
   
   return clus;
}
