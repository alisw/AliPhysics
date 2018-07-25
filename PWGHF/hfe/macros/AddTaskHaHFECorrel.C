AliAnalysisTaskHaHFECorrel *AddTaskHaHFECorrel(Double_t period, Double_t MinPtEvent, Double_t MaxPtEvent, Bool_t TRDQA, Bool_t TagEff, Bool_t RecEff, Bool_t CorrHadron, Bool_t CorrLP, Bool_t MCTruth,  Bool_t IsMC, Bool_t IsAOD, Bool_t IsHFE, Bool_t UseTender, Double_t EtaMax, Int_t ITSnCut,  Int_t TPCnCut, Int_t TPCnCutdEdx,   Double_t PhotElecPtCut, Int_t PhotElecTPCnCut,Bool_t PhotElecITSrefitCut,Double_t InvmassCut, Int_t HTPCnCut,   Bool_t HITSrefitCut, Bool_t HTPCrefitCut, Bool_t UseITS, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, const char * ID="")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHaHFECorrel", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHaHFECorrel", "This task requires an input event handler");
    return 0x0;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();

  /*
  AliMCEventHandler* mcHand = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHand);
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if (!mcH) {
    MCthere=kFALSE;
  }
  */



  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHaHFECorrel.C");
  AliAnalysisTaskHaHFECorrel *taskMB = 
    ConfigHaHFECorrel(period, MinPtEvent, MaxPtEvent, TRDQA, TagEff, RecEff,  CorrHadron, CorrLP, MCTruth, IsMC, IsAOD, UseTender,EtaMax, ITSnCut, TPCnCut, TPCnCutdEdx, PhotElecPtCut,PhotElecTPCnCut, PhotElecITSrefitCut,  InvmassCut,  HTPCnCut,  HITSrefitCut, HTPCrefitCut, UseITS, SigmaITScut, SigmaTOFcut, SigmaTPCcut, ID);
  if (!taskMB) {
    Error("AddTaskHaHFECorrel", "No task found.");
  }
  taskMB->SelectCollisionCandidates(AliVEvent::kINT7);
  
  // Load correction weights for pi0, eta
  if (IsMC) {
    TH1::AddDirectory(kFALSE);
    printf("Loading Pi0EtaCorrectionFiles\n");
    TString CorrectPi0EtaFile;
    if (!IsHFE) CorrectPi0EtaFile="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/Pi0EtaWeightsMinBias.root";
    else if (IsHFE) CorrectPi0EtaFile="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/Pi0EtaWeightsHFE.root";
    TFile *CorrectPi0Eta = TFile::Open(CorrectPi0EtaFile.Data());
    if (CorrectPi0Eta) {    
      TH1F * Pi0W = (TH1F*)CorrectPi0Eta->Get("Pi0WeightsIncl");
      TH1F * EtaW = (TH1F*)CorrectPi0Eta->Get("EtaWeightsIncl");
      if (Pi0W) taskMB->SetPi0WeightToData(*Pi0W);
      else printf("Could not load Pi0Weights\n");
      if (EtaW)  taskMB->SetEtaWeightToData(*EtaW);
      else printf("Could not load EtaWeights\n");
    }
    else printf("Could not open Pi0Eta correction file \n");
    TH1::AddDirectory(kTRUE);
  }




  TH1::AddDirectory(kFALSE);
  printf("Loading SPDnTrAvg\n");
  TString SPDnTrFileName;
  if (IsMC) SPDnTrFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/SPDnTrAvg_MC.root";
  else SPDnTrFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/SPDnTrAvg_Data.root";
  TFile *SPDnTrFile  = TFile::Open(SPDnTrFileName.Data());
  if (SPDnTrFile) {    
    TProfile* SPDnTrAvg = (TProfile*)SPDnTrFile->Get("SPDnTrAvg");
    if (SPDnTrAvg) taskMB->SetSPDnTrAvg(*SPDnTrAvg);
    else printf("Could not load SPDnTrAvg\n");
    }
    else printf("Could not open SPDnTrAvg correction file \n");
  TH1::AddDirectory(kTRUE);

  TH1::AddDirectory(kFALSE);
  printf("Loading RecEffFiles\n");
  TString RecEffFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/RecEff.root";
  TFile *RecEffFile = TFile::Open(RecEffFileName.Data());
  RecEffFile->ls();
  //RecEffFile->ls();
  if (RecEffFile) {    
    TH3F * HadRecEff = (TH3F*)RecEffFile->Get("HadRecEff");
    TH2F * EleRecEff = (TH2F*)RecEffFile->Get("EleRecEff");
    if (HadRecEff) taskMB->SetHadRecEff(*HadRecEff);
    else printf("Could not load HadRecEff\n");
    if (EleRecEff) taskMB->SetEleRecEff(*EleRecEff);
    else printf("Could not load EleRecEff\n");
  }
  else printf("Could not open RecEff correction file \n");
  TH1::AddDirectory(kTRUE);




  mgr->AddTask(taskMB);

  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_HaHFECorrel_kINT7_";
  containerName1 += ID;
        
  TString name1 = "histMB_";
  name1 += ID;
        
  TString name2 = "Main_";
  name2 += ID;
  TString name3 = "LP_";
  name3 += ID;
  TString name4 = "Hadron_";
  name4 += ID;
  TString name5 = "QA_";
  name5 += ID;


  AliAnalysisDataContainer *cinput      = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1    = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputMain = mgr->CreateContainer(name2.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputLP   = mgr->CreateContainer(name3.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputHa   = mgr->CreateContainer(name4.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputQA   = mgr->CreateContainer(name5.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());

  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);
  mgr->ConnectOutput(taskMB, 2, coutputMain);
  mgr->ConnectOutput(taskMB, 3, coutputLP);
  mgr->ConnectOutput(taskMB, 4, coutputHa);
  mgr->ConnectOutput(taskMB, 5, coutputQA);



  return NULL;
}
