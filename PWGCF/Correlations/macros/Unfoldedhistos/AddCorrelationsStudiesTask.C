#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TROOT.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCorrelationsStudies.h"
#endif

AliAnalysisTaskSE *AddCorrelationsStudiesTask(const char *suffix, const char *configstring, const char *corrconfigstring, const char *corrbinning) {

  TString szContainerPrefix;

  AliAnalysisTaskCorrelationsStudies* taskCS;
  AliAnalysisManager    *mgr = AliAnalysisManager::GetAnalysisManager();
  TString outfilename = mgr->GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  if (mgr != NULL) {

    taskCS = new AliAnalysisTaskCorrelationsStudies(Form("Correlation Studies %s", suffix));
    taskCS->Configure(configstring);
    if (TString(corrconfigstring).Contains("simulate"))
      taskCS->ConfigureCorrelations(corrconfigstring,Form("density_%s_", suffix),szContainerPrefix);
    else
      taskCS->ConfigureCorrelations(corrconfigstring,Form("correction_%s_", suffix),szContainerPrefix);
    taskCS->ConfigureCorrelationsBinning(corrbinning);
    mgr->AddTask(taskCS);

    // create containers for input/output
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("CorrelationStudies_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("CSDptDptCorrelations_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("CSDptDptCorrelationsMCRecOptions_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(Form("CSDptDptTrueCorrelations_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);

    // connect input/output
    mgr->ConnectInput(taskCS, 0, cinput);
    mgr->ConnectOutput(taskCS, 1, coutput1);
    mgr->ConnectOutput(taskCS, 2, coutput2);
    mgr->ConnectOutput(taskCS, 3, coutput3);
    mgr->ConnectOutput(taskCS, 4, coutput4);

    return taskCS;
  }
  else {
    AliFatalGeneral("AddCorrelationsStudiesTask", "Add task needs a previously started analysis manager");
    return NULL;
  }
}
