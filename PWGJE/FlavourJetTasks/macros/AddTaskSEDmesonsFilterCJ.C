// $Id$

AliAnalysisTaskSEDmesonsFilterCJ *AddTaskSEDmesonsFilterCJ(TString sCutDzero = "", TString sCutDstar = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskSEDmesonsFilterCJ.C::AddTaskSEDmesonsFilterCJ", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskSEDmesonsFilterCJ.C::AddTaskSEDmesonsFilterCJ", "Task manager to have an ESD or AOD input handler.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskSEEAddTaskSEDmesonsFilterCJ.C::AddTaskSEEAddTaskSEDmesonsFilterCJ", "This task requires an input event handler");
    return NULL;
  }
//=============================================================================

  TFile *file = 0;
  AliRDHFCutsD0toKpi      *cutDzero = 0;
  AliRDHFCutsDStartoKpipi *cutDstar = 0;

  if (!sCutDzero.IsNull()) {
    file = TFile::Open(sCutDzero.Data(), "READ");
    cutDzero = dynamic_cast<AliRDHFCutsD0toKpi*>(file->Get("D0toKpiCuts"));
    file->Close();
  }

  if (!sCutDstar.IsNull()) {
    file = TFile::Open(sCutDstar.Data(), "READ");
    cutDstar = dynamic_cast<AliRDHFCutsDStartoKpipi*>(file->Get("DStartoKpipiCuts"));
    file->Close();
  }

  AliAnalysisTaskSEDmesonsFilterCJ *taskDemsonFilterCJ = new AliAnalysisTaskSEDmesonsFilterCJ("AliAnalysisTaskSEDmesonsFilterCJ");
  taskDemsonFilterCJ->SetCutDzero(cutDzero);
  taskDemsonFilterCJ->SetCutDstar(cutDstar);
  mgr->AddTask(taskDemsonFilterCJ);

  mgr->ConnectInput( taskDemsonFilterCJ, 0, mgr->GetCommonInputContainer());
//mgr->ConnectOutput(taskDemsonFilterCJ, 0, mgr->GetCommonOutputContainer());

  TString sCommon  = AliAnalysisManager::GetCommonFileName();
  TString sCutName = "AnalysisCutsDmeson.root";
  TString sOutput1 = "AliAnalysisTaskSEDmesonsFilterCJ_ListDmesonCuts";
  TString sOutput2 = "AliAnalysisTaskSEDmesonsFilterCJ_Dzero_ControlHistograms";
  TString sOutput3 = "AliAnalysisTaskSEDmesonsFilterCJ_Dstar_ControlHistograms";
  mgr->ConnectOutput(taskDemsonFilterCJ, 1, mgr->CreateContainer(sOutput1.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,sCutName.Data()));
  mgr->ConnectOutput(taskDemsonFilterCJ, 2, mgr->CreateContainer(sOutput2.Data(),TList::Class(),AliAnalysisManager::kOutputContainer, sCommon.Data()));
  mgr->ConnectOutput(taskDemsonFilterCJ, 3, mgr->CreateContainer(sOutput3.Data(),TList::Class(),AliAnalysisManager::kOutputContainer, sCommon.Data()));
  return taskDemsonFilterCJ;
}
