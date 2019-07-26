/////////////////////////////////////////////////////////////////////////
//                                                                     //            
// AddTaskUeRTDep                                                      //
// Author:Aditya NAth Mishra (aditya.mishra@correo.nucleares.unam.mx)  //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
 class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskUeRTDep()
{
  // set authomatically to true if MC
  Bool_t AnalysisMC= (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  // set to true if you wish to build corrections (memory consuming)
  Bool_t AnalysisCorr = kTRUE;
  UInt_t AnaTrigger = AliVEvent::kINT7;

  // Get the pointer to the existing analysis manager via the static access method. 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
    ::Error("AddTaskUeRTDep", "No analysis manager to connect to.");
      return NULL; 
    }
// Check the analysis type using the event handlers connected to the analysis manager.  
//==============================================================================
  if(!mgr->GetInputEventHandler()){
   ::Error("AddTaskUeRTDep", "This task requires an input event handler");
      return NULL;   
   }


  gROOT->LoadMacro("CreateTrackCutsPWGJE.C");
  
  AliAnalysisTaskUeRTDep * Mytask = new AliAnalysisTaskUeRTDep("AnalysisUeRTDep");
   if(!Mytask) return 0x0;
   TString kInputDataType = mgr->GetInputEventHandler()->GetDataType(); 

   Bool_t is13TeV = kTRUE;
   Mytask->SetAnalysisType(kInputDataType);
   Mytask->SetTrigger(AnaTrigger);
   Mytask->SetDebugLevel(0);
   Mytask->SetEtaCut(0.8);
   Mytask->SetVtxCut(10.0);
   Mytask->SetPileUpRej(kTRUE);	

   mgr->AddTask(Mytask);

   TString fileName = AliAnalysisManager::GetCommonFileName();
   TString kName  = "outputUeRTDep";
 
   mgr->ConnectInput(Mytask,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(Mytask,1,mgr->CreateContainer(kName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
   return Mytask;
}
