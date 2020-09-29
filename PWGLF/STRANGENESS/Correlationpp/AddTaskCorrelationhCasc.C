//used to instantiate an object of the task,define input and output and connect it to manager

AliAnalysisTaskCorrelationhCasc* AddTaskCorrelationhCasc(TString name = "name", Float_t minpt=3, Float_t maxpt=15, bool isLocal=kTRUE, bool isMC=kTRUE, bool isEff=kTRUE, Int_t EvtToMix=50, Float_t EtaTrigger=0.8, Float_t EtahAssoc=0.8,Float_t EtaV0Assoc=0.8, Int_t FilterBitValue=128, Int_t year=2010, TString AssocParticle="Xi"){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }

  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  TString fileName;  
  if(isLocal) {
    fileName = "AnalysisResultsLocalhCasc.root";
  }
  else {
    fileName = AliAnalysisManager::GetCommonFileName();
  }
  fileName += ":MyTask";     
  fileName += AssocParticle;

  AliAnalysisTaskCorrelationhCasc* task = new AliAnalysisTaskCorrelationhCasc(name.Data());   
  if(!task) return 0x0;
  //task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  // add your task to the manager
   task->SetMinPt(minpt); 
   task->SetMaxPt(maxpt);
   task->SetMC(isMC);
   task->SetEff(isEff);
   task->SetEvtToMix(EvtToMix);
   task->SetEtaTrigger(EtaTrigger);
   task->SetEtahAssoc(EtahAssoc);
   task->SetEtaV0Assoc(EtaV0Assoc);
   task->SetFilterBit(FilterBitValue);
   task->SetYear(year);
   task->SetAssocParticle(AssocParticle);

   //   Float_t fMassLimInf=0;
   //   Float_t fMassLimSup=0;

  mgr->AddTask(task);

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer("MyOutputContainer1", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,3,mgr->CreateContainer("MyOutputContainer2", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,4,mgr->CreateContainer("MyOutputContainer3", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,5,mgr->CreateContainer("MyOutputContainer4", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,6,mgr->CreateContainer("Risoluzione", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  return task;
}
