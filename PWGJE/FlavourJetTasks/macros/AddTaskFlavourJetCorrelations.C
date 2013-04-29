// $Id$

AliAnalysisTaskFlavourJetCorrelations *AddTaskFlavourJetCorrelations(
  TString sUsedTrks    = "",
  TString sUsedClus    = "",
  TString sUsedJets    = "",
  TString sUsedRho     = "",
  TString sTaskName    = "AliAnalysisTaskFlavourJetCorrelations",
  Bool_t  bIsMakeHisto = kTRUE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskFlavourJetCorrelations.C::AddTaskFlavourJetCorrelations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlavourJetCorrelations.C::AddTaskFlavourJetCorrelations", "This task requires an input event handler");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskFlavourJetCorrelations.C::AddTaskFlavourJetCorrelations", "Task manager to have an ESD or AOD input handler.");
    return NULL;
  }
//=============================================================================

  AliAnalysisTaskFlavourJetCorrelations *taskFlavourCJ = new AliAnalysisTaskFlavourJetCorrelations(sTaskName.Data(),bIsMakeHisto);
  taskFlavourCJ->SetClusName(sUsedClus.Data());
  taskFlavourCJ->SetTracksName(sUsedTrks.Data());
  taskFlavourCJ->SetJetsName(sUsedJets.Data());
  taskFlavourCJ->SetRhoName(sUsedRho.Data());

  mgr->AddTask(taskFlavourCJ);
  mgr->ConnectInput( taskFlavourCJ, 0, mgr->GetCommonInputContainer());
//mgr->ConnectOutput(taskFlavourCJ, 0, mgr->GetCommonOutputContainer());

  TString sOutput0 = sTaskName + "_" + sUsedJets;
  TString sOutput1 = sOutput0  + "_GeneralHistograms";
  TString sOutput2 = sOutput0  + "_ControlHistograms";
  TString sOutput3 = sOutput0  + "_AnDzeroHistograms";
  TString sOutput4 = sOutput0  + "_AnDstarHistograms";
  TString sCommon  = AliAnalysisManager::GetCommonFileName();
  mgr->ConnectOutput(taskFlavourCJ, 1, mgr->CreateContainer(sOutput1.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,sCommon.Data()));
  mgr->ConnectOutput(taskFlavourCJ, 2, mgr->CreateContainer(sOutput2.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,sCommon.Data()));
  mgr->ConnectOutput(taskFlavourCJ, 3, mgr->CreateContainer(sOutput3.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,sCommon.Data()));
  mgr->ConnectOutput(taskFlavourCJ, 4, mgr->CreateContainer(sOutput4.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,sCommon.Data()));

  return taskFlavourCJ;
}
