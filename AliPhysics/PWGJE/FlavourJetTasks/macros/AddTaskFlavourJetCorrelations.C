AliAnalysisTaskFlavourJetCorrelations *AddTaskFlavourJetCorrelations(
  Bool_t theMCon = kFALSE,
  Bool_t reco = kTRUE /*must be true if theMCon is false*/,
  TString jetArrname = "",
  UInt_t pSel=AliVEvent::kAny,
  TString trigClass="",
  TString suffix = "",
  TString trackArrname = "PicoTracks",
  Bool_t triggerOnLeadingJet = kFALSE,
  Int_t leadingHadType = 0 /*0 = charged, 1 = neutral, 2 = both*/,
  Float_t R = 0.4,
  Float_t jptcut = 10.,
  const char *cutType = "TPC",
  Double_t percjetareacut = 1.,
  Bool_t bJetOnly=kTRUE
  )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
     ::Error("AddTaskFlavourJetCorrelations::AddTaskFlavourJetCorrelations", "No analysis manager to connect to.");
     return NULL;
  } 
  if(bJetOnly){//only jet histograms, no correlation with D mesons
     
     //used for event selection
     AliRDHFCuts *dummyDcut=new AliRDHFCutsD0toKpi("dummy");
     dummyDcut->SetTriggerMask(pSel);
     dummyDcut->SetTriggerClass(trigClass);
     AliAODPidHF *pidHF=new AliAODPidHF();
     dummyDcut->SetPidHF(pidHF);
     
     Int_t dummycand=AliAnalysisTaskFlavourJetCorrelations::kD0toKpi;
     
     // create the task
     AliAnalysisTaskFlavourJetCorrelations *task = new AliAnalysisTaskFlavourJetCorrelations("AnaTaskFlavourJetCorrelations", 
     	dummyDcut, dummycand, bJetOnly);
     task->SetJetArrayName(jetArrname);
     task->SetTrackArrayName(trackArrname);
     task->SetMC(theMCon);
     task->SetUseReco(reco);
     task->SetTriggerOnLeadingJet(triggerOnLeadingJet);
     task->SetJetAcceptanceType(cutType);
     task->SetJetPtCut(jptcut);
     task->SetPercAreaCut(percjetareacut);
     mgr->AddTask(task);

     if(theMCon) {
     	suffix+="MC";
     	if(reco) suffix+="rec";  
     }
     
     // Create and connect containers for input/output
     TString outputfile = AliAnalysisManager::GetCommonFileName();
     outputfile += ":PWGJE_HF_DEmcalJet";
     outputfile += jetArrname;
     outputfile += trigClass;
     TString nameContainer1="listJets";
     TString nameContainer2="evselCut";
     //nameContainer1 += candname;
     
     nameContainer1 += suffix;
     nameContainer2 += suffix;
     
     // ------ input data ------
     AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
     
     // ----- output data -----
     AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(nameContainer1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
     AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(nameContainer2, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

     mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(task,1,coutput1);
     mgr->ConnectOutput(task,2,coutput2);

     
  } else{
     AliFatal("ERROR: Use AddTaskDFilterAndCorrelations instead");
     
  }

  Printf("Input and Output connected to the manager");
  return task ;
}
