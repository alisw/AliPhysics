AliAnalysisTaskJetServices* serv = 0;
Bool_t ConfigWithFlagsJetServices();
AliAnalysisTaskJetServices *AddTaskJetServices()
{
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetServices", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
     ::Error("AddTaskJetServices", "This task requires an input event handler");
      return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   type.ToUpper();
   // Create the task and configure it.
   //===========================================================================
   
   serv = new AliAnalysisTaskJetServices("JetServices");
      

   if(type == "AOD"){
     serv->SetAODInput(kTRUE);
   }
   mgr->AddTask(serv);

   // evaluate global variables 
   Bool_t bGood1 = false;
   Bool_t bGood2 = false;


   if(!ConfigWithFlagsJetServices())return 0;
   serv->SetUsePhysicsSelection(kTRUE);
     
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Serv = mgr->CreateContainer("pwgje_services", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_services",AliAnalysisManager::GetCommonFileName()));

   mgr->ConnectInput  (serv, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (serv, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (serv,  1, coutput1_Serv );
   
   return serv;
}

Bool_t ConfigWithFlagsJetServices(){
    
  Bool_t bGood1 = kFALSE;
  Bool_t bGood2 = kFALSE;


  serv->SetRunRange(AliAnalysisManager::GetGlobalInt("kGridRunRangeLo",bGood1),
		      AliAnalysisManager::GetGlobalInt("kGridRunRangeLo",bGood2));
  
  if(!bGood1||!bGood2){
    Printf("%s:%d Run range not set",(char*)__FILE__,__LINE__);
    serv->SetRunRange(110000,160000);
  }
   

  Int_t nTrigger = AliAnalysisManager::GetGlobalInt("kNTrigger",bGood1);
  
  if(bGood1){
    serv->SetNTrigger(nTrigger);
    for(int it = 0;it < nTrigger;it++){
      serv->SetTrigger(it,
		       AliAnalysisManager::GetGlobalInt(Form("kTriggerBit%d",it),bGood1),
		       AliAnalysisManager::GetGlobalStr(Form("kTriggerName%d",it),bGood2));
     }
   }
     
   AliAnalysisManager::GetGlobalInt("kPhysicsSelectionFlag",bGood1);if(bGood1)serv->SetPhysicsSelectionFlag(AliAnalysisManager::GetGlobalInt("kPhysicsSelectionFlag",bGood1));
   else {Printf("%s%d: kPhysicsSelectionFlag not defined",(char*)__FILE__,__LINE__); return kFALSE; }
   AliAnalysisManager::GetGlobalStr("kDeltaAODJetName",bGood1);if(bGood1)serv->SetNonStdFile(AliAnalysisManager::GetGlobalStr("kDeltaAODJetName",bGood1));
   else {Printf("%s%d: kDeltaAODJetName not defined",(char*)__FILE__,__LINE__); return kFALSE; }

   AliAnalysisManager::GetGlobalDbl("kTrackEtaWindow",bGood1);if(bGood1)serv->SetTrackEtaWindow(AliAnalysisManager::GetGlobalDbl("kTrackEtaWindow",bGood1));
   else {Printf("%s%d: kTrackEtaWindow not defined",(char*)__FILE__,__LINE__); return kFALSE; }
   AliAnalysisManager::GetGlobalDbl("kVertexWindow",bGood1);if(bGood1)serv->SetZVertexCut(AliAnalysisManager::GetGlobalDbl("kVertexWindow",bGood1));
   else {Printf("%s%d: kVertexWindow not defined",(char*)__FILE__,__LINE__); return kFALSE; }
   
   AliAnalysisManager::GetGlobalInt("kHighPtFilterMask",bGood1);if(bGood1)serv->SetFilterMask(AliAnalysisManager::GetGlobalInt("kHighPtFilterMask",bGood1));
   else {Printf("%s%d: kHighPtFilterMask not defined",(char*)__FILE__,__LINE__); return kFALSE; }   

   TString cRun(AliAnalysisManager::GetGlobalStr("kJetRunPeriod",bGood1));
   if(cRun.Contains("10h")||cRun.Contains("11h")){
     serv->SetCollisionType(AliAnalysisTaskJetServices::kPbPb);
   }
   else{
     serv->SetCollisionType(AliAnalysisTaskJetServices::kPP);
   }
   
   return kTRUE;

}
