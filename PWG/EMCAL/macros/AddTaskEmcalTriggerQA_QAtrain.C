/// \file AddTaskEmcalTriggerQA_QAtrain.C
/// \brief This is the AddTask macro of AliEmcalTriggerQATask
///
/// This macro is used in a LEGO train to add an instance
/// of AliEmcalTriggerQATask in the analysis manager.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
///\modif M. Germin
/// \date July 4 , 2017

AliEmcalTriggerQATask* AddTaskEmcalTriggerQA_QAtrain(const Int_t runnumber=253681)
{
  return AliEmcalTriggerQATask::AddTaskEmcalTriggerQA_QAtrain(runnumber);
}
