AliAnalysisTaskJetCluster *AddTaskJetCluster(char* bRec = "AOD",char* bGen = "",UInt_t filterMask = 16, Int_t iPhysicsSelection = 1,Char_t *jf = "KT", Float_t radius = 0.4);


AliAnalysisTaskJetCluster *AddTaskJetClusterDelta(UInt_t filterMask = 16,Bool_t kUseAODMC = kFALSE,Int_t iPhysicsSelection = 1,Char_t *jf = "KT", UInt_t iFlag){
  AliAnalysisTaskJetCluster *js = 0;
  if(kUseAODMC&&false){// do not use the MC info yet
    if(iFlag&(1<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,0.00001);
    if(iFlag&(2<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,0.1);
    if(iFlag&(3<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,0.2);
    if(iFlag&(4<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,0.4);
    if(iFlag&(5<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,0.6);
    if(iFlag&(6<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,0.8);
    if(iFlag&(7<<0))js = AddTaskJetCluster("AOD","AODMC",filterMask,iPhysicsSelection,jf,1.0);
  }
  else{
    if(iFlag&(1<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,0.00001);
    if(iFlag&(2<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,0.1);
    if(iFlag&(3<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,0.2);
    //    if(iFlag&(4<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,0.4);
    if(iFlag&(5<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,0.6);
    if(iFlag&(6<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,0.8);
    if(iFlag&(7<<0))js = AddTaskJetCluster("AOD","",filterMask,iPhysicsSelection,jf,1.0);
  }
  
  return js;
}


AliAnalysisTaskJetCluster *AddTaskJetCluster(char* bRec,char* bGen ,UInt_t filterMask,Int_t iPhysicsSelection,Char_t *jf,Float_t radius)
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

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
   typeGen.ToUpper();
   typeRec.ToUpper();
   // Create the task and configure it.
   //===========================================================================


   char *cRadius = "";
   if(radius>0)cRadius = Form("%02d",(int)((radius+0.01)*10.));

   
   AliAnalysisTaskJetCluster* pwg4spec = new  AliAnalysisTaskJetCluster(Form("JetCluster_%s_%s",jf,cRadius));
      
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
   }
   else if (typeRec.Contains("AODMC2")){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCCharged);
   }
   else if (typeRec.Contains("AODMC")){
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAODMCAll);
   }
   else if (typeRec.Contains("AOD")) {
     pwg4spec->SetTrackTypeRec(AliAnalysisTaskJetCluster::kTrackAOD);
   }

   if(typeGen.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetCluster::kTrackAODMCChargedAcceptance);
   }
   else if (typeGen.Contains("AODMC2")){
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetCluster::kTrackAODMCCharged);
   }
   else if (typeGen.Contains("AODMC")){
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetCluster::kTrackAODMCAll);
   }
   else if (typeGen.Contains("AOD")) {
     pwg4spec->SetTrackTypeGen(AliAnalysisTaskJetCluster::kTrackAOD);
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

   
   if(iPhysicsSelection)pwg4spec->SelectCollisionCandidates();

   mgr->AddTask(pwg4spec);
     
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer(Form("pwg4cluster_%s_%s_%s_%s",bRec,bGen,jf,cRadius), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_cluster_%s_%s_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen,jf,cRadius));

   mgr->ConnectInput  (pwg4spec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4spec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );
   
   return pwg4spec;
}
