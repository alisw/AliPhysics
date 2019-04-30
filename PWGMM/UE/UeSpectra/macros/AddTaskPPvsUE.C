///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskPPvsUE                                                 //
// Author: Valentina Zaccolo: valentina.zaccolo@cern.ch          //
//                                                               //
///////////////////////////////////////////////////////////////////
 class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskPPvsUE()
{

  Bool_t AnalysisMC= (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  

  // Get the pointer to the existing analysis manager via the static access method. 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
    ::Error("AddTaskPPvsUE", "No analysis manager to connect to.");
      return NULL; 
    }
// Check the analysis type using the event handlers connected to the analysis manager.  
//==============================================================================
  if(!mgr->GetInputEventHandler()){
   ::Error("AddTaskPPvsUE", "This task requires an input event handler");
      return NULL;   
   }


  gROOT->LoadMacro("CreateTrackCutsPWGJE.C");

   AliAnalysisTaskPPvsUE * task = new AliAnalysisTaskPPvsUE("AnalysisPPvsUE");
   if(!task) return 0x0;
   TString kInputDataType = mgr->GetInputEventHandler()->GetDataType(); 

   Bool_t is13TeV = kTRUE;
   task->SetAnalysisMC(AnalysisMC);
   task->SetAnalysisType(kInputDataType);
   task->SetDebugLevel(0);
   task->SetEtaCut(0.8);
   task->SetVtxCut(10.0);
   task->SetPileUpRej(kTRUE);	
   task->SetAveMultiInTrans(6.026); // VZ: to be checked

   mgr->AddTask(task);

   TString fileName = AliAnalysisManager::GetCommonFileName();
   TString kName  = "outputPPvsUE";
 
   mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task,1,mgr->CreateContainer(kName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
   return task;
}
