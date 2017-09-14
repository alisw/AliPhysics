AliAnalysisTaskPHOSTimeCalib* AddTaskPHOSTimeCalib(
    const char* name     = "PHOSCalibration",
    const UInt_t trigger = AliVEvent::kAny
    )
{
  //Add a task AliPHOSCalibration to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTimeCalib", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTimeCalib", "This task requires an input event handler");
    return NULL;
  }

	TString TriggerName="";
	if     (trigger == (UInt_t)AliVEvent::kAny)      TriggerName = "kAny";
	else if(trigger == (UInt_t)AliVEvent::kINT7)     TriggerName = "kINT7";
	else if(trigger == (UInt_t)AliVEvent::kPHI7)     TriggerName = "kPHI7";
	else                                             TriggerName = "ELSE";


  AliAnalysisTaskPHOSTimeCalib* task = new AliAnalysisTaskPHOSTimeCalib(Form("%sTask_%s",name,TriggerName.Data()));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  TString cname(Form("hist_PHOSCalibration_%s",TriggerName.Data()));
  TString outputFile = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), THashList::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}

