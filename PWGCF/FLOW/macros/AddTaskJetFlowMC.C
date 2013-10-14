///////////////////////////////////////////////////////////////////////////////
//                         AddTaskJetFlowToyMC                               //
//   Author: Redmer A. Bertens, Utrecht University, 2013, rbertens@cern.ch   //
///////////////////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliAnalysisTaskJetFlowMC;

AliAnalysisTaskJetFlowMC* AddTaskJetFlowMC(
  const char *outputTracks      = "JetFlowToyMC",
  const char *inputTracks       = "PicoTracks",
  const char *name              ="AliAnalysisTaskJetFlowMC"
  )
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)                             return 0x0;
  if (!mgr->GetInputEventHandler())     return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":";
  fileName += name;
        // create the task
  AliAnalysisTaskJetFlowMC *eTask = new AliAnalysisTaskJetFlowMC(name);
  eTask->SetTracksOutName(outputTracks);
  eTask->SetTracksInName(inputTracks);
        // connect input and output
  mgr->AddTask(eTask);
  mgr->ConnectInput  (eTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (eTask, 1, mgr->CreateContainer(Form("%s_container", fileName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return eTask;
}
