AliAnalysisTaskHaHFECorrel * ConfigHaHFECorrel(Bool_t UseTender, Double_t period, Double_t AssPtCut, Int_t ITSnCut, Int_t AssTPCnCut, Int_t TPCnCut, Int_t HTPCnCut, Bool_t AssITSrefitCut, Bool_t HITSrefitCut, Bool_t HTPCrefitCut, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, Bool_t rejectKinkMother, Bool_t CorrHadron, Bool_t CorrLP, Bool_t OpeningAngleCut, Double_t InvmassCut)
{

  //AliHFEcuts *hfecuts = new AliHFECuts("name", "title");
  //hfecuts->

  AliAnalysisTaskHaHFECorrel *task = new AliAnalysisTaskHaHFECorrel("HaHFECorrel");
  printf("task ------------------------ %p\n ", task);
  task->SetTender(UseTender);
  task->SetPeriod(period);
  task->SetAssPtCut(AssPtCut);
  task->SetITSnCut(ITSnCut);
  task->SetAssTPCnCut(AssTPCnCut);
  task->SetTPCnCut(TPCnCut);
  task->SetHTPCnCut(HTPCnCut);
  task->SetAssITSrefitCut(AssITSrefitCut);
  task->SetHITSrefitCut(HITSrefitCut);
  task->SetHTPCrefitCut(HTPCrefitCut);
  task->SetSigmaITScut(SigmaITScut);
  task->SetSigmaTOFcut(SigmaTOFcut);
  task->SetSigmaTPCcut(SigmaTPCcut);
  task->SetRejectKinkMother(rejectKinkMother);
  task->SetHadronCorrelation(CorrHadron);
  task->SetLPCorrelation(CorrLP);
  task->SetOpeningAngleCut(OpeningAngleCut);
  task->SetInvmassCut(InvmassCut);
 // task->SetHFECuts(hfecuts);

  // Define PID
  //AliHFEpid *pid = task->GetPID();
  //if(useMC) pid->SetHasMCData(kTRUE);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  //task->PrintStatus();
  //pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
