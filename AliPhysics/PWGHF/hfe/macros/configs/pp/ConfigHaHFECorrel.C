
AliAnalysisTaskHaHFECorrel * ConfigHaHFECorrel(Int_t period, Int_t MinNTr, Int_t MaxNTr, Bool_t TRDQA, Bool_t TagEff, Bool_t RecEff, Bool_t OneTimeCheck, Bool_t CorrHadron, Bool_t CorrLP, Bool_t MCTruth,  Bool_t IsMC,Bool_t IsAOD, Bool_t IsHFE, Bool_t UseTender, Bool_t UseEventWeights, Double_t EtaMax, Int_t ITSnCut, Float_t ITSSharedCluster, Int_t TPCnCut, Int_t TPCnCutdEdx, Double_t PhotElecPtCut, Int_t PhotElecTPCnCut, Bool_t PhotElecITSrefitCut, Int_t PhotCorrCase, Double_t InvmassCut, Int_t HTPCnCut,Bool_t HITSrefitCut, Bool_t HTPCrefitCut,Bool_t UseITSsa, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, Int_t VarOptE, Int_t VarOptH, Int_t VarOptPhot, const char * ID="")
{

  //AliHFEcuts *hfecuts = new AliHFECuts("name", "title");
  //hfecuts->

  TString Name;
  Name = Form("HaHFECorrel_%s", ID);
  //Name+=Form("_%i", VarOptE);
  printf("NAME ---------- %s", Name.Data());
  AliAnalysisTaskHaHFECorrel *task = new AliAnalysisTaskHaHFECorrel(Name);
  printf("task ------------------------ %p\n ", task);
  task->SetPeriod(period);
  printf("\nRunning over period %i\n", period);

  printf("\nEventSelection based on NTracklets: %i <= NTr <= %i", MinNTr, MaxNTr);
  task->SetMinNTr(MinNTr);
  task->SetMaxNTr(MaxNTr);
  
  printf("\nConfigList: "); 
  //GeneralSettings
  TString Config;
  Config= Form("%i, %i, %i, %s, %s, %s, %s", period, MinNTr, MaxNTr, TRDQA ? "kTRUE" : "kFALSE",  TagEff  ? "kTRUE" : "kFALSE",  RecEff ? "kTRUE" : "kFALSE",  OneTimeCheck ? "kTRUE" : "kFALSE");
  Config+=Form(", %s, %s, %s, %s, %s, %s, %s, %s",  CorrHadron ? "kTRUE" : "kFALSE",  CorrLP ? "kTRUE" : "kFALSE",  MCTruth ? "kTRUE" : "kFALSE",  IsMC ? "kTRUE" : "kFALSE",  IsAOD ? "kTRUE" : "kFALSE",  IsHFE ? "kTRUE" : "kFALSE", UseTender ? "kTRUE" : "kFALSE", UseEventWeights ? "kTRUE" : "kFALSE");
  Config+=Form(", %4.2f, %i, %4.2f, %i, %i, %4.2f, %i, %s, %i, %4.2f, %i, %s, %s, %s, %4.2f, %4.2f, %4.2f, %i, %i, %i,  %s", EtaMax, ITSnCut, ITSSharedCluster, TPCnCut, TPCnCutdEdx, PhotElecPtCut, PhotElecTPCnCut, PhotElecITSrefitCut ? "kTRUE" : "kFALSE", PhotCorrCase, InvmassCut, HTPCnCut, HITSrefitCut ? "kTRUE" : "kFALSE", HTPCrefitCut ? "kTRUE" : "kFALSE", UseITSsa ? "kTRUE" : "kFALSE", SigmaITScut, SigmaTOFcut, SigmaTPCcut, VarOptE, VarOptH, VarOptPhot, ID );
  printf("%s",Config.Data()); 

  if (TRDQA) printf("\nPerforming TRDQA");
  task->SetTRDQA(TRDQA);
  if (TagEff) printf ("\nGenerate TagEff Histograms");
  task->SetTagEff(TagEff);
 if (RecEff) printf ("\nGenerate RecEff histograms");
  task->SetRecEff(RecEff);
  if (OneTimeCheck) printf ("\nGenerate OneTimeCheck histograms for the chossen setting");
  task->SetOneTimeCheck(OneTimeCheck);
  if (CorrHadron) printf("\nCorrelating Hadrons");
  task->SetHadronCorrelation(CorrHadron);
  if (CorrLP) printf("\nCorrelating LP");
  task->SetLPCorrelation(CorrLP);
  if (MCTruth) printf("\nCorrelating MCTruh");
  task->SetMCTruthCorrelation(MCTruth);

  task->SetMC(IsMC);
  printf("\nIs MC? %i", IsMC);
  task->SetAODanalysis(IsAOD);
  printf("\nIs AOD? %i", IsAOD);
  task->SetTender(UseTender);
  printf("\nUse Tender? %i", UseTender);
  task->SetUseEventWeights(UseEventWeights);
  printf("\nUse EventWeights? %i", UseEventWeights);
  task->SetEtaMax(EtaMax);
  printf("\nElectron |EtaMax| %f", EtaMax);

  task->SetITSnCut(ITSnCut);
  printf("\nElectron ITSNclsCut: %i", ITSnCut);

  task->SetITSSharedClusterCut(ITSSharedCluster);
  printf("\nElectron ITSSharedClusterCut: %f", ITSSharedCluster); 

  task->SetTPCnCut(TPCnCut);
  printf("\nElectron TPCNclsCut: %i", TPCnCut);

  task->SetTPCnCutdEdx(TPCnCutdEdx);
  printf("\nElectron TPCNclsdEdxCut: %i", TPCnCutdEdx);

  task->SetPhotElecPtCut(PhotElecPtCut);
  printf("\nPhotElec PtCut: %f", PhotElecPtCut);

  task->SetPhotElecTPCnCut(PhotElecTPCnCut);
  printf("\nPhotElec TPCNclsCut: %i", PhotElecTPCnCut);
  
  task->SetPhotElecITSrefitCut(PhotElecITSrefitCut);
  printf("\nPhotElec ITSrefit? %i", PhotElecITSrefitCut);

  task->SetPhotCorrCase(PhotCorrCase);
  printf("\nPhotCorrCase? %i", PhotCorrCase);

  task->SetInvmassCut(InvmassCut);
  printf("\nInvMass Cut: %f", InvmassCut);

  task->SetHTPCnCut(HTPCnCut);
  printf("\nHadron TPCNclsCut: %i", HTPCnCut);

  task->SetHITSrefitCut(HITSrefitCut);
  printf("\nHadron ITSrefit? %i", HITSrefitCut);

  task->SetHTPCrefitCut(HTPCrefitCut);
  printf("\nHadron TPCrefit? %i", HTPCrefitCut);

 
  task->SetUseITSsa(UseITSsa);
  printf("\nUse ITS? %i", UseITSsa);

  task->SetSigmaITScut(SigmaITScut);
  printf("\nITS nSigma: %f", SigmaITScut);
  task->SetSigmaTOFcut(SigmaTOFcut);
  printf("\nTOF nSigma: %f", SigmaTOFcut);
  task->SetSigmaTPCcut(SigmaTPCcut);
  printf("\nTPC nSigma min: %f\n", SigmaTPCcut);


  // task->SetOpeningAngleCut(OpeningAngleCut);
 // task->SetHFECuts(hfecuts);

  // Define PID
  //AliHFEpid *pid = task->GetPID();
  //if(useMC) pid->SetHasMCData(kTRUE);

  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  //printf("MC is %i", IsMC);
  //task->PrintStatus();
  //pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
