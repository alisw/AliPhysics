/// \file AddTaskEmcalTriggerQA.C
/// \brief This is the AddTask macro of AliEmcalTriggerQATask
///
/// This macro is used in a LEGO train to add an instance
/// of AliEmcalTriggerQATask in the analysis manager.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr 4, 2016

AliEmcalTriggerQATask* AddTaskEmcalTriggerQA(
    const char* triggerPatchesName  = "EmcalTriggers",
    const char* cellsName           = 0,
    const char* triggersName        = 0,
    AliEmcalTriggerQATask::EBeamType_t beamType = AliEmcalTriggerQATask::kpp,
    AliEmcalTriggerQATask::ETriggerAnalysisType_t anaType = AliEmcalTriggerQATask::kTriggerOfflineExpertAnalysis,
    const char* suffix              = "")
{
  return AliEmcalTriggerQATask::AddTaskEmcalTriggerQA(triggerPatchesName, cellsName, triggersName, beamType, anaType, "", suffix);
}
