AliAnalysisTaskCheckEvSel *AddTaskCheckEvSel(TString suffix="",
					     Int_t system=1,
					     Bool_t readMC=kFALSE,
					     TString filecutName="",
					     TString cutObjname="",
					     ULong64_t trigMask=AliVEvent::kMB |AliVEvent::kINT7,
					     Double_t minCent=0,
					     Double_t maxCent=100,
					     Int_t cutOnZVertexSPD=0,
					     Int_t optPileup=AliRDHFCuts::kRejectPileupEvent,
					     Int_t minContPileup=3,
					     Double_t minDzPileup=0.6,
					     Bool_t multDepPileup=kFALSE,
					     Int_t optForOOBPileupInPbPb=0,
					     Bool_t useCentrCorrelCuts=kFALSE,
					     Bool_t doEvPrNtuple=kFALSE,
					     Bool_t doVtxNtuple=kFALSE)
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCombinHF", "No analysis manager to connect to.");
  }
  
  //Analysis Task
  AliRDHFCutsD0toKpi* evselCuts = 0x0;

  if(filecutName.EqualTo("")){
    evselCuts=new AliRDHFCutsD0toKpi("EvSelCuts");
    evselCuts->SetTriggerClass("");
    evselCuts->SetTriggerMask(trigMask);
    if(minCent>=0 && minCent<=110 && maxCent>=0 && maxCent<=110){
      evselCuts->SetUseCentrality(1);
      evselCuts->SetMinCentrality(minCent);
      evselCuts->SetMaxCentrality(maxCent);
    }
    evselCuts->SetOptPileup(optPileup);
    if(optPileup==AliRDHFCuts::kRejectPileupEvent){
      evselCuts->ConfigureSPDPileupCuts(minContPileup,minDzPileup);
      evselCuts->SetUseMultDepPileupCut(multDepPileup);
    }
    else if(optPileup==AliRDHFCuts::kRejectMVPileupEvent){
      evselCuts->SetMinContribPileupMV(minContPileup);
    }
    evselCuts->SetCutOnzVertexSPD(cutOnZVertexSPD);
    evselCuts->SetUseCentralityCorrelationCuts(useCentrCorrelCuts);
    evselCuts->SetUsePbPbOutOfBunchPileupCut(optForOOBPileupInPbPb);
  }else{
    TFile* filecuts=TFile::Open(filecutName.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      ::Fatal("AddTaskCheckEvSel","Input file with cuts not found");
    }
    evselCuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutObjname.Data());
    if(!evselCuts){
      ::Fatal("AddTaskCheckEvSel","Cut object not found");
    }
  }
  
  AliAnalysisTaskCheckEvSel *dTask = new AliAnalysisTaskCheckEvSel(readMC,system,evselCuts);
  dTask->SetCutOnzVertexSPD(cutOnZVertexSPD);
  dTask->SetEnableVertexNtuple(doVtxNtuple);
  dTask->SetEnableEvPropNtuple(doEvPrNtuple);
  mgr->AddTask(dTask);
  
  
  // Create containers for input/output
  
  TString inname = "cinpEvSelCheck";
  TString outname = "coutEvSelCheck";
  TString normname = "coutNorm";
  TString nt1name = "coutEvPropNtuple";
  TString nt2name = "coutZvertNtuple";

  inname += suffix.Data();
  outname += suffix.Data();
  normname += suffix.Data();
  nt1name += suffix.Data();
  nt2name += suffix.Data();
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":D2H_CheckEvSel";
  outputfile += suffix.Data();
  
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname,TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  
  AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputNt1 = mgr->CreateContainer(nt1name,TNtuple::Class(),
							      AliAnalysisManager::kOutputContainer,
							      outputfile.Data());
  coutputNt1->SetSpecialOutput();

  AliAnalysisDataContainer *coutputNt2 = mgr->CreateContainer(nt2name,TNtuple::Class(),
							      AliAnalysisManager::kOutputContainer,
							      outputfile.Data());
  coutputNt2->SetSpecialOutput();
  
  mgr->ConnectInput(dTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dTask,1,coutput);
  mgr->ConnectOutput(dTask,2,coutputNorm);
  mgr->ConnectOutput(dTask,3,coutputNt1);
  mgr->ConnectOutput(dTask,4,coutputNt2);
  
  return dTask;
}
