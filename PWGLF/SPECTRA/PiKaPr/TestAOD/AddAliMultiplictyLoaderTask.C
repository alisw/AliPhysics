#if defined(__CLING__)
#include "AliMultiplictyLoaderTask.h"
#endif

AliMultiplictyLoaderTask* AddAliMultiplictyLoaderTask(Bool_t useAliPPVsMultUtils=kTRUE, TString cent="V0M")
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 	if (!mgr) 
   	{
   	   ::Error("AliMultiplictyLoaderTask", "No analysis manager to connect to.");
   	   return NULL;
   	 }   
  
  	  if (!mgr->GetInputEventHandler()) 
  	  {
  	    ::Error("AliMultiplictyLoaderTask", "This task requires an input event handler");
  	    return NULL;
  	  } 
   	  TString type = mgr->GetInputEventHandler()->GetDataType(); 
   	  if(type.Contains("AOD"))
    	  {
      		::Error("AliMultiplictyLoaderTask", "This task requires to run on ESD");
      		return NULL;
    	  }
	
 	 AliMultiplictyLoaderTask *task = new  AliMultiplictyLoaderTask("AliMultiplictyLoaderTask");
 	 task->SetUseAliPPVsMultUtils(useAliPPVsMultUtils);
 	 task->SetCentEstimator(cent);	
	 TString outputFileName = AliAnalysisManager::GetCommonFileName();
 	 AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
   	 mgr->ConnectInput(task, 0, cinput);
 	 mgr->AddTask(task);

 	 return task;

}
