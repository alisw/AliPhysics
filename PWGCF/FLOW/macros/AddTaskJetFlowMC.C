///////////////////////////////////////////////////////////////////////////////
//                         AddTaskJetFlowToyMC                               //
//   Author: Redmer A. Bertens, Utrecht University, 2013, rbertens@cern.ch   //
///////////////////////////////////////////////////////////////////////////////

/* AddTask macro for jet flow toy mc task
 * task uses an afterburner to tune vn in the pico track 
 * selection which can be used by a jet finder 
 * note that this task does not generate MC particles, it changes
 * the azimuthal distribution of already available tracks
*/

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
