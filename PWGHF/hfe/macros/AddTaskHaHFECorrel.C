
AliAnalysisTask *AddTaskHaHFECorrel(Bool_t UseTender, Double_t period, Double_t AssPtCut, Int_t ITSnCut, Int_t AssTPCnCut, Int_t TPCnCut, Int_t HTPCnCut, Bool_t AssITSrefitCut, Bool_t HITSrefitCut, Bool_t HTPCrefitCut, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, Bool_t rejectKinkMother, Bool_t CorrHadron, Bool_t CorrLP, Bool_t OpeningAngleCut, Double_t InvmassCut, TString ID="ContName")
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
    ConfigHaHFECorrel(UseTender, period, AssPtCut, ITSnCut, AssTPCnCut, TPCnCut, HTPCnCut, AssITSrefitCut, HITSrefitCut, HTPCrefitCut, SigmaITScut, SigmaTOFcut, SigmaTPCcut, rejectKinkMother, CorrHadron, CorrLP, OpeningAngleCut, InvmassCut);
  if (!taskMB) {
    Error("AddTaskHaHFECorrel", "No task found.");
  }
  taskMB->SelectCollisionCandidates(AliVEvent::kINT7);
    
  mgr->AddTask(taskMB);

  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_HaHFECorrel_kINT7";
  containerName1 += ID;
        
  TString name1 = "histMB";
  name1 += ID;
        
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);




  return NULL;
}
