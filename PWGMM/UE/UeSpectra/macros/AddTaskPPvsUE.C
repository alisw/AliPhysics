///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskPPvsUE                                                 //
// Author: Valentina Zaccolo: valentina.zaccolo@cern.ch          //
//                                                               //
///////////////////////////////////////////////////////////////////
 class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskPPvsUE()
{
  // set authomatically to true if MC
  Bool_t AnalysisMC= (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  // set to true if you wish to build corrections (memory consuming)
  Bool_t AnalysisCorr = kTRUE; 

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
   task->SetAnalysisCorr(AnalysisCorr);
   task->SetAnalysisType(kInputDataType);
   task->SetDebugLevel(0);
   task->SetEtaCut(0.8);
   task->SetVtxCut(10.0);
   task->SetPileUpRej(kTRUE);	
   task->SetAveMultiInTrans(4.939);
   task->SetAveGenMultiInTrans(7.392); // this is PYTHIA 8

   mgr->AddTask(task);

   TString fileName = AliAnalysisManager::GetCommonFileName();
   TString kName  = "outputPPvsUE";
 
   mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task,1,mgr->CreateContainer(kName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
   return task;
}
