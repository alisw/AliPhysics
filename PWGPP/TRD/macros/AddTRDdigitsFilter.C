//=============================================================================
//
// *** AddTRDdigitsFilter
//
// This macro sets up the TRD digits filter to write the digits for a
// PID reference sample of electrons and pions from V0 decays.
//
//=============================================================================

#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisTask.h>
#include <AliTRDdigitsFilter.h>
#endif

AliAnalysisTask  *AddTRDdigitsFilter(Int_t runNumber)
{
  //gSystem->Load("libTRDrec");
  // pointer to the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTRDdigitsFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }  

  /////////////////////////
  // The TRD filter Task
  /////////////////////////
  AliTRDdigitsFilter *filterTask = new AliTRDdigitsFilter();

  mgr->AddTask(filterTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  if (!cinput) cinput = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput =mgr->CreateContainer("TRDdigitsFilter",TList::Class(), AliAnalysisManager::kOutputContainer, "DigitsFilter.root");  


  mgr->ConnectInput(filterTask,0,cinput);
  mgr->ConnectOutput(filterTask,1,coutput);
  return filterTask;

}


