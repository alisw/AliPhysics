Bool_t ConfigWithFlagsJetSpectrum2();
AliAnalysisTaskJetSpectrum2 *jetspec = 0;

AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(const char* bRec = "jets",const char* bGen = "jetsAODMC_UA104",const char* nonStdFile="",UInt_t filterMask = 32, Int_t iPhysicsSelectionFlag = AliVEvent::kMB,Int_t hjet=0, UInt_t iEventSelectionMask = 0,Int_t iCl = 0, Bool_t bRCSparseDimensions = kFALSE);


AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2Delta(UInt_t filterMask = 32,Bool_t kUseAODMC = kFALSE,Int_t iPhysicsSelectionFlag = AliVEvent::kMB,Int_t hjet=0,UInt_t iFlag = 0xfffffff, UInt_t iEventSelectionMask = 0,char* back = ""){

  TString cBack = back;

  AliAnalysisTaskJetSpectrum2 *js = 0;
  if(kUseAODMC){
    if(iFlag&(1<<0)){ // UA104
      js = AddTaskJetSpectrum2("jets","jetsAODMC_UA104",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
      js = AddTaskJetSpectrum2("jets","jetsAODMC2_UA104",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
    }
    if(iFlag&(1<<1)){ // ANTIKT 04
      js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","jetsAODMC_FASTJET04",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
      js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","jetsAODMC2_FASTJET04",cBack.Data(),filterMask,iPhysicsSelectionFlag,  iEventSelectionMask);
      // cross check MC only background subtration
      js = AddTaskJetSpectrum2("jetsAODMC2_FASTJET04","jetsAODMC_FASTJET04",cBack.Data(),filterMask,iPhysicsSelectionFlag,  iEventSelectionMask);
    }
    if(iFlag&(1<<2)){ // KT 04
      js = AddTaskJetSpectrum2("jetsAOD_FASTKT04","jetsAODMC_FASTKT04",cBack.Data(),filterMask,iPhysicsSelectionFlag,iEventSelectionMask);
      js = AddTaskJetSpectrum2("jetsAOD_FASTKT04","jetsAODMC2_FASTKT04",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
    }
    if(iFlag&(1<<3)){ // SISCONE 04
      js = AddTaskJetSpectrum2("jetsAOD_SISCONE04","jetsAODMC_SISCONE04",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
      js = AddTaskJetSpectrum2("jetsAOD_SISCONE04","jetsAODMC2_SISCONE04",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
    }
    // here can go other radii
  }
  else { // only the data ... no MC
    if(iFlag&(1<<0)){ // UA104
      js = AddTaskJetSpectrum2("jets","",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask,1);
    }
    if(iFlag&(1<<1)){ // ANTIKT 04
      js = AddTaskJetSpectrum2("jetsAOD_FASTJET04","",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
    }
    if(iFlag&(1<<2)){ // KT 04
      js = AddTaskJetSpectrum2("jetsAOD_FASTKT04","",cBack.Data(),filterMask,iPhysicsSelectionFlag,iEventSelectionMask);
    }
    if(iFlag&(1<<3)){ // SISCONE 04
      js = AddTaskJetSpectrum2("jetsAOD_SISCONE04","",cBack.Data(),filterMask,iPhysicsSelectionFlag, iEventSelectionMask);
    }
  }
  return js;
}


AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2(const char* bRec,const char* bGen ,const char* nonStdFile,UInt_t filterMask,Int_t iPhysicsSelectionFlag,Int_t hjet, UInt_t iEventSelectionMask,Int_t iCl, Bool_t bRCSparseDimensions)
{
  // Creates a jet fider task, configures it and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
     ::Error("AddTaskJetSpectrum2", "No analysis manager to connect to.");
    return NULL;
   }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetSpectrum2", "This task requires an input event handler");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  TString typeRec(bRec);
  TString typeGen(bGen);
  typeGen.ToUpper();
  typeRec.ToUpper();
  

   // Create the task and configure it.
   //===========================================================================

   jetspec = new  AliAnalysisTaskJetSpectrum2(Form("JetSpectrum2%s-%s_%010d_Class%02d",bRec,bGen,iEventSelectionMask,iCl));
   if(iCl)jetspec->SetEventClass(iCl);

   // add the filter mask for non default jets branches
   TString cAdd("");
   cAdd += Form("_Filter%05d",filterMask);

   
   jetspec->SetTRP(hjet);
   jetspec->SetBranchGen(bGen); 
   //  if(typeGen.Contains("JETSAOD")&&!typeGen.Contains("MC"))jetspec->SetBranchGen(Form("%s%s",bGen,cAdd.Data())); 

   jetspec->SetBranchRec(bRec); 
   // if(typeRec.Contains("JETSAOD")&&!typeRec.Contains("MC"))     jetspec->SetBranchRec(Form("%s%s",bRec,cAdd.Data())); 



   if(type == "AOD"){
     // Assume all jets are not yet produced 
     //     jetspec->SetAODJetInput(kTRUE);
     jetspec->SetAODTrackInput(kTRUE);
     jetspec->SetAODMCInput(kTRUE);
   }
   else{
     if(mgr->GetMCtruthEventHandler()){
       jetspec-> SetAnalysisType(AliAnalysisTaskJetSpectrum2::kAnaMCESD);
     }
   }

   if(typeRec.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     jetspec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODMCChargedAcceptance);
   }
   else if (typeRec.Contains("AODMC2")){
     jetspec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODMCCharged);
   }
   else if (typeRec.Contains("AODMC")){
     jetspec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAODMCAll);
   }
   else { // catch akk use AOD
     jetspec->SetTrackTypeRec(AliAnalysisTaskJetSpectrum2::kTrackAOD);
   }

   if(typeGen.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     jetspec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCChargedAcceptance);
   }
   else if (typeGen.Contains("AODMC2")){
     jetspec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCCharged);
   }
   else if (typeGen.Contains("AODMC")){
     jetspec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAODMCAll);
   }
   else if (typeGen.Length()>0){ // catch all use AOD
     jetspec->SetTrackTypeGen(AliAnalysisTaskJetSpectrum2::kTrackAOD);
   }

   if(iEventSelectionMask)jetspec->SetEventSelectionMask(iEventSelectionMask);

   //   jetspec->SetDebugLevel(10);
   if(!ConfigWithFlagsJetSpectrum2())return 0;

   jetspec->SetUseGlobalSelection(kTRUE); 
   jetspec->SetMinJetPt(-1.);//keep all jets 

   if(bRCSparseDimensions) {
     jetspec->SetNRPBins(100);
     jetspec->SetNMultBins(1);
     jetspec->SetNPtLeadingBins(1);
   }
   
   // to fetch the AOD from the AOD extension ouput 
   if(strlen(nonStdFile))jetspec->SetNonStdFile(nonStdFile);
   mgr->AddTask(jetspec);
     
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = ;
   if(hjet>0)
     coutput1_Spec = mgr->CreateContainer(Form("pwgje_spec2_%s_%s_%010d_Class%02d_HJ%d",bRec,bGen,iEventSelectionMask,iCl,hjet),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_spec2_%s_%s_%010d_Class%02d_HJ%d",AliAnalysisManager::GetCommonFileName(),bRec,bGen,iEventSelectionMask,iCl,hjet));
   else
     coutput1_Spec = mgr->CreateContainer(Form("pwgje_spec2_%s_%s_%010d_Class%02d",bRec,bGen,iEventSelectionMask,iCl),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_spec2_%s_%s_%010d_Class%02d",AliAnalysisManager::GetCommonFileName(),bRec,bGen,iEventSelectionMask,iCl));

   mgr->ConnectInput  (jetspec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (jetspec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (jetspec,  1, coutput1_Spec );
   
   return jetspec;
}

void SetAODInput(AliAnalysisTaskJetSpectrum2 *taskJetSpectrum){
  taskJetSpectrum->SetAODJetInput(kTRUE);
  taskJetSpectrum->SetAODTrackInput(kTRUE);
  // taskJetSpectrum->SetUseGlobalSelection(kFALSE);
}

Bool_t ConfigWithFlagsJetSpectrum2(){
    
  Bool_t bGood1 = kFALSE;
  Bool_t bGood2 = kFALSE;

  
  Int_t nTrigger = AliAnalysisManager::GetGlobalInt("kNTrigger",bGood1);
  
  if(bGood1){
    jetspec->SetNTrigger(nTrigger);
    for(int it = 0;it < nTrigger;it++){
      jetspec->SetTrigger(it,
			  AliAnalysisManager::GetGlobalInt(Form("kTriggerBit%d",it),bGood1),
			  AliAnalysisManager::GetGlobalStr(Form("kTriggerName%d",it),bGood2));
    }
   }
   else {
     Printf("%s%d: kNTrigger not defined",(char*)__FILE__,__LINE__); return kFALSE; 
   }


  Int_t nAcceptance = AliAnalysisManager::GetGlobalInt("kNAcceptanceSpec",bGood1);
  
  if(bGood1){
    jetspec->SetNAcceptance(nAcceptance);
    for(int ia = 0;ia < nAcceptance;ia++){
      jetspec->SetAcceptance(ia,
			     AliAnalysisManager::GetGlobalDbl(Form("kAcceptancePhiMinSpec%d",ia),bGood1),
			     AliAnalysisManager::GetGlobalDbl(Form("kAcceptancePhiMaxSpec%d",ia),bGood1),
			     AliAnalysisManager::GetGlobalDbl(Form("kAcceptanceEtaMinSpec%d",ia),bGood1),
			     AliAnalysisManager::GetGlobalDbl(Form("kAcceptanceEtaMaxSpec%d",ia),bGood1));
    }
   }
   else {
     Printf("%s%d: kNAcceptance not defined",(char*)__FILE__,__LINE__); 
   }

     
   AliAnalysisManager::GetGlobalInt("kPhysicsSelectionFlag",bGood1);if(bGood1){
     jetspec->SelectCollisionCandidates(AliAnalysisManager::GetGlobalInt("kPhysicsSelectionFlag",bGood1));
   }
   else {Printf("%s%d: kPhysicsSelectionFlag not defined",(char*)__FILE__,__LINE__); return kFALSE; }


   AliAnalysisManager::GetGlobalStr("kDeltaAODJetName",bGood1);if(bGood1)jetspec->SetNonStdFile(AliAnalysisManager::GetGlobalStr("kDeltaAODJetName",bGood1));
   else {Printf("%s%d: kDeltaAODJetName not defined",(char*)__FILE__,__LINE__); return kFALSE; }

   AliAnalysisManager::GetGlobalDbl("kTrackEtaWindow",bGood1);if(bGood1)jetspec->SetTrackEtaWindow(AliAnalysisManager::GetGlobalDbl("kTrackEtaWindow",bGood1));
   else {Printf("%s%d: kTrackEtaWindow not defined",(char*)__FILE__,__LINE__); return kFALSE; }

   AliAnalysisManager::GetGlobalDbl("kJetEtaWindow",bGood1);if(bGood1)jetspec->SetJetEtaWindow(AliAnalysisManager::GetGlobalDbl("kJetEtaWindow",bGood1));
   else {Printf("%s%d: kJetEtaWindow not defined",(char*)__FILE__,__LINE__); return kFALSE; }

   
   AliAnalysisManager::GetGlobalInt("kHighPtFilterMask",bGood1);if(bGood1)jetspec->SetFilterMask(AliAnalysisManager::GetGlobalInt("kHighPtFilterMask",bGood1));
   else {Printf("%s%d: kHighPtFilterMask not defined",(char*)__FILE__,__LINE__); return kFALSE; }   
   
   return kTRUE;


}
