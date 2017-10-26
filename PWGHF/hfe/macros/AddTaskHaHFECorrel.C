AliAnalysisTask *AddTaskHaHFECorrel(Double_t period,Bool_t CorrHadron, Bool_t CorrLP,  Bool_t IsMC, Bool_t UseTender, Int_t ITSnCut,  Int_t TPCnCut, Int_t TPCnCutdEdx,   Double_t PhotElecPtCut, Int_t PhotElecTPCnCut,Bool_t PhotElecITSrefitCut,Double_t InvmassCut, Int_t HTPCnCut,   Bool_t HITSrefitCut, Bool_t HTPCrefitCut, Bool_t UseITS, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, TString ID="")
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

  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if (!mcH) {
    MCthere=kFALSE;
  }


  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHaHFECorrel.C");
  AliAnalysisTaskHaHFECorrel *taskMB = 
    ConfigHaHFECorrel(period, CorrHadron, CorrLP, IsMC, UseTender, ITSnCut, TPCnCut, TPCnCutdEdx, PhotElecPtCut,PhotElecTPCnCut, PhotElecITSrefitCut,  InvmassCut,  HTPCnCut,  HITSrefitCut, HTPCrefitCut, UseITS, SigmaITScut, SigmaTOFcut, SigmaTPCcut, ID);
  if (!taskMB) {
    Error("AddTaskHaHFECorrel", "No task found.");
  }
  taskMB->SelectCollisionCandidates(AliVEvent::kINT7);
  
  if (IsMC) {
    TH1::AddDirectory(kFALSE);
    printf("Loading Pi0EtaCorrectionFiles\n");
    TString CorrectPi0EtaFile="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/Pi0EtaWeights.root";
    TFile *CorrectPi0Eta = TFile::Open(CorrectPi0EtaFile.Data());
    if (CorrectPi0Eta) {    
      TH1F * Pi0W = (TH1F*)CorrectPi0Eta->Get("Pi0Weights");
      TH1F * EtaW = (TH1F*)CorrectPi0Eta->Get("EtaWeights");
      if (Pi0W) taskMB->SetPi0WeightToData(*Pi0W);
      else printf("Could not load Pi0Weights\n");
      if (EtaW)  taskMB->SetEtaWeightToData(*EtaW);
      else printf("Could not load EtaWeights\n");
    }
    else printf("Could not open Pi0Eta correction file \n");
    TH1::AddDirectory(kTRUE);
  }

  mgr->AddTask(taskMB);

  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_HaHFECorrel_kINT7_";
  containerName1 += ID;
        
  TString name1 = "histMB_";
  name1 += ID;
        
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);




  return NULL;
}
