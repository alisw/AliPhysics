#include "AliAnalysisTaskBSEmbedding.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"   
#include "TClonesArray.h" 
#include "AliAnalysisManager.h"                                                                                         
AliAnalysisTaskBSEmbedding::AliAnalysisTaskBSEmbedding()                                   
  : AliAnalysisTaskSE ("AliAnalysisTaskBSEmbedding")                                           
{                                                                                                  
}                                                                                                  
                                                                                                   
AliAnalysisTaskBSEmbedding::AliAnalysisTaskBSEmbedding (const char *name, const char *option)
  : AliAnalysisTaskSE(name)                                                                        
  , fOption(option)                                                                                
{                                                                                                  
}                                                                                                  
AliAnalysisTaskBSEmbedding::AliAnalysisTaskBSEmbedding(const AliAnalysisTaskBSEmbedding& ap) 
	: AliAnalysisTaskSE(ap.fName)                           
  , fOption(ap.fOption)                                               
{                                                                     
}                                                                     
AliAnalysisTaskBSEmbedding& AliAnalysisTaskBSEmbedding::operator = (const AliAnalysisTaskBSEmbedding& ap) 
{                                                                     
                                                                      
  this->~AliAnalysisTaskBSEmbedding();                                            
  new(this) AliAnalysisTaskBSEmbedding(ap);                                       
  return *this;                                                       
}                                                                     
AliAnalysisTaskBSEmbedding::~AliAnalysisTaskBSEmbedding()   
{                                   
}                                   


AliAnalysisTaskBSEmbedding* AliAnalysisTaskBSEmbedding::AddTaskBSEmbedding()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalEmbeddingHelper", "No analysis manager to connect to.");
    return 0;
  }
                                                                                                    
  // Check the analysis type using the event handlers connected to the analysis manager.            
  //==============================================================================                  
  AliVEventHandler* handler = mgr->GetInputEventHandler();                                          
  if (!handler)                                                                                     
  {                                                                                                 
    ::Error("AddTaskEmcalEmbeddingHelper", "This task requires an input event handler");            
    return 0;                                                                                       
  }                                                                                                 
                                                                                                    
  TString name = "AliAnalysisTaskBSEmbedding";                                             
                                                                                                    
  AliAnalysisTaskBSEmbedding * mgrTask = static_cast<AliAnalysisTaskBSEmbedding*>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;                                                                      
  
  // Create the task that manages                                                                   
  AliAnalysisTaskBSEmbedding* embeddingHelper = new AliAnalysisTaskBSEmbedding(name.Data(),name.Data());
  
  //-------------------------------------------------------                                         
  // Final settings, pass to manager and set the containers                                         
  //-------------------------------------------------------                                         
  
  mgr->AddTask(embeddingHelper);
                                                                                                    
  // Create containers for input/output                                                             
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();                                
  
      
  mgr->ConnectInput(embeddingHelper, 0, cInput);
  
  return embeddingHelper;
}                                                                                                   


                                                                                                   
void AliAnalysisTaskBSEmbedding::UserExec(Option_t *option){                                   
                                                                                                   
  AliVEvent *event = InputEvent();                                                                 
	AliVEvent *event_emb = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent(); 
	TClonesArray* fMCArray = (TClonesArray*) event_emb->FindListObject("mcparticles");             
	event -> AddObject(fMCArray);                                                                  
                                                                                                   
}                                                                                                  

