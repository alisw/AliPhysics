// $Id$

AliAnalysisTaskRhoBase* AddTaskRhoBase(
   const char *rhoname        = "Rho",
   TF1        *rfunc          = 0
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoBase", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRho", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("Rho_Base"));
  AliAnalysisTaskRho *rhotask = new AliAnalysisTaskRho(name);
  rhotask->SetRhoName(rhoname);
  rhotask->SetRhoFunction(rfunc);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput (rhotask, 0, mgr->GetCommonInputContainer() );

  return rhotask;
}
