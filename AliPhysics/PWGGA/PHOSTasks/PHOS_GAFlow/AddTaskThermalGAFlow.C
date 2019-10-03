AliAnalysisTaskThermalGAFlow* AddTaskThermalGAFlow(
   const char *outfilename    = "AnalysisOutput.root",

//Default Cuts
   const Int_t fDebug = 0,
   const Int_t fMinCells = 3,
   const Double_t fMinE = 0.3,
   const Double_t fMinTrackDr = 0,
   const Double_t fMaxVertexx = 10,
   const Double_t fMinCentrality = -1,
   const Double_t fMaxCentrality = 100,
   const Double_t fCoreRadius = 3.5,
   const Double_t fMinCoreEnergyRatio = 0.4,
   const Double_t fMinLambdaDisp = 0.3,
   const Double_t fMinCPVStd = 2.5,

   const Int_t fMixVertxbins = 1,
   const Int_t fMixCentbins = 1,
   const Int_t fMixEvbins = 1,
   const Int_t fNptbins = 150,
//End Default Cuts

   const char *tag            = ""

)

{    
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskPHOSThermalGAFlow", "No analysis manager to connect to.");
    return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskPHOSThermalGAFlow", "This task requires an input event handler");
    return NULL;
  }
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name(Form("ash_%s", tag));
  printf("Adding\n");
  AliAnalysisTaskThermalGAFlow *PHOSGAtask = new AliAnalysisTaskThermalGAFlow(name);

  PHOSGAtask->SetDebug(fDebug);
  PHOSGAtask->SetMinCells(fMinCells);
  PHOSGAtask->SetMinE(fMinE);
  PHOSGAtask->SetMinTrackDr(fMinTrackDr);
  PHOSGAtask->SetMaxVertexx(fMaxVertexx);
  PHOSGAtask->SetMinCentrality(fMinCentrality);
  PHOSGAtask->SetMaxCentrality(fMaxCentrality);
  PHOSGAtask->SetCoreRadius(fCoreRadius);
  PHOSGAtask->SetMinCoreEnergyRatio(fMinCoreEnergyRatio);
  PHOSGAtask->SetMinLambdaDisp(fMinLambdaDisp);
  PHOSGAtask->SetMinCPVStd(fMinCPVStd);

  PHOSGAtask->SetMixVertxbins(fMixVertxbins);
  PHOSGAtask->SetMixCentbins(fMixCentbins);
  PHOSGAtask->SetMixEvbins(fMixEvbins);
  PHOSGAtask->SetNptbins(fNptbins);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(PHOSGAtask);
  // Create containers for input/output
  mgr->ConnectInput (PHOSGAtask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coGAFlow = mgr->CreateContainer(name,
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    outfilename);

  mgr->ConnectOutput(PHOSGAtask,1,coGAFlow);
  return PHOSGAtask;
}

//Note: Order of global variables matter!  They must match header, addtask, and cxx.
