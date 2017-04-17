/// \file AddTaskEMCALTriggerQA.C
/// \ingroup EMCALPerformanceMacros
/// \brief Configuration analysis task for EMCal trigger QA.
///
/// Configuration analysis task for EMCal trigger QA. The only parameter is the check on
/// MC or real data.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

AliAnalysisTaskEMCALTriggerQA * AddTaskEMCALTriggerQA(Bool_t kSimulation = kFALSE, TString outputFile = "")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTaskEMCALTriggerQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTaskEMCALTriggerQA", "This task requires an input event handler");
    return NULL;
  }
  
  TString taskName = "QAEMCALTrigger";
  
  AliAnalysisTaskEMCALTriggerQA * qatrigger = new AliAnalysisTaskEMCALTriggerQA(taskName);
    
  // Configuration
  
  Bool_t kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
    
  if(kUseKinematics) kSimulation = kTRUE;
    
  if(kSimulation) qatrigger->SwitchOnMCData();
  else            qatrigger->SwitchOffMCData();
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
  if(inputDataType.Contains("AOD")) qatrigger->SwitchOffV0SignalHistograms();
  else                              qatrigger->SwitchOnV0SignalHistograms();
  
  // Define the input/output containers
    
  qatrigger->GetRecoUtils()->SwitchOnBadChannelsRemoval ();
  
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    
  if(outputFile.Length()==0) outputFile = AliAnalysisManager::GetCommonFileName();
 
  printf("*** Task Name %s; Output file name: %s ***\n",taskName.Data(),outputFile.Data());
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s",taskName.Data()), 
                                                           TList::Class(), AliAnalysisManager::kOutputContainer,
                                                           Form("%s",outputFile.Data()));

  mgr->AddTask(qatrigger);
  mgr->ConnectInput  (qatrigger, 0, cinput1);
  mgr->ConnectOutput (qatrigger, 1, coutput);
  
  return qatrigger;
}
