//used to instantiate an object of the task,define input and output and connect it to manager

AliAnalysisTaskCorrelationhhK0s* AddTaskCorrelationhhK0s(Float_t minpt=3, Float_t maxpt=15, Bool_t ishhCorr,bool isLocal=kTRUE, bool isMC=kTRUE, bool isEff=kTRUE, bool isHybridMCTr = kTRUE, Int_t EvtToMix=50, Float_t EtaTrigger=0.8, Float_t EtahAssoc=0.8,Float_t EtaV0Assoc=0.8, Int_t FilterBitValue=128, Int_t year=2010, TString AssocParticle="K0s", TString suffix = ""){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }

  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  TString fileName;  
  if(isLocal) {
    fileName = "AnalysisResultsLocal.root";
  }
  else {
    fileName = AliAnalysisManager::GetCommonFileName();
  }
  fileName += ":MyTask";     
  if (AssocParticle!="K0s")  fileName += AssocParticle;
  if (isMC && !isEff && !isHybridMCTr) fileName += "_MCTruth";
  if (isMC && isHybridMCTr) fileName += "_MCHybrid";
  fileName += Form("_PtTrigMin%.1f", minpt);
  fileName += Form("_PtTrigMax%.1f", maxpt);

  TString combinedName=Form("hK0s_Task_%s", suffix.Data());
  AliAnalysisTaskCorrelationhhK0s* task = new AliAnalysisTaskCorrelationhhK0s(combinedName);   
  if(!task) return 0x0;
  //task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  // add your task to the manager
   task->SetMinPt(minpt); 
   task->SetMaxPt(maxpt);
   task->SetCorr(ishhCorr);
   task->SetMC(isMC);
   task->SetEff(isEff);
   task->SetHybridMCTruth(isHybridMCTr);
   task->SetEvtToMix(EvtToMix);
   task->SetEtaTrigger(EtaTrigger);
   task->SetEtahAssoc(EtahAssoc);
   task->SetEtaV0Assoc(EtaV0Assoc);
   task->SetFilterBit(FilterBitValue);
   task->SetYear(year);
   task->SetAssocParticle(AssocParticle);

  mgr->AddTask(task);

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("MyOutputContainer_%s",combinedName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("MyOutputContainer1_%s",combinedName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("MyOutputContainer2_%s",combinedName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("MyOutputContainer3_%s",combinedName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("MyOutputContainer4_%s",combinedName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,6,mgr->CreateContainer(Form("Risoluzione_%s",combinedName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  return task;
}
