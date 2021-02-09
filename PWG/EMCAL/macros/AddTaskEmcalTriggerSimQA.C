/// \file AddTaskEmcalTriggerSimQA.C
/// \brief This is the AddTask macro of AliEmcalTriggerSimQATask
///
/// This macro is used in a LEGO train to add an instance
/// of AliEmcalTriggerSimQATask in the analysis manager.
///
/// \author Michael Oliver <michael.oliver@cern.ch>, Yale University
/// \date Dec 11, 2018

AliEmcalTriggerSimQATask* AddTaskEmcalTriggerSimQA()
{
  return AliEmcalTriggerSimQATask::AddTaskEmcalTriggerSimQA();
}
