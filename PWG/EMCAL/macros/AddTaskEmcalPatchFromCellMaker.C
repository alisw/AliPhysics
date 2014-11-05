// $Id$

AliEmcalPatchFromCellMaker* AddTaskEmcalPatchFromCellMaker(
  Int_t       patchDim            = 32,
  Double_t    cellMinE            = 0.05,
  Bool_t      bSlideL1            = kFALSE,
  const char *patchOutName        = "EmcalPatches",
  const char *cellsName           = 0,
  const char *taskName            = "PatchFromCellMaker")
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalPatchFromCellMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalPatchFromCellMaker", "This task requires an input event handler");
    return NULL;
  }

  TString strCellsName(cellsName);
  if(strCellsName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strCellsName = "EMCALCells";
      ::Info("AddTaskEmcalPatchFromCellMaker", Form( "ESD analysis, cellsName = \"%s\"", strCellsName.Data() ));
    }
    else {
      strCellsName = "emcalCells";
      ::Info("AddTaskEmcalPatchFromCellMaker", Form( "AOD analysis, cellsName = \"%s\"", strCellsName.Data() ));
    }
  }

  TString strPatchOutName = Form("%s%dx%d",patchOutName,patchDim,patchDim);
  TString name = Form("%s_%s",taskName,strPatchOutName.Data());
 
   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalPatchFromCellMaker *eTask = new AliEmcalPatchFromCellMaker(name);
  eTask->SetCaloTriggersOutName(strPatchOutName.Data());
  eTask->SetCaloCellsName(strCellsName.Data());
  eTask->SetPatchDimension(patchDim);
  eTask->SetMinCellE(cellMinE);
  eTask->ActivateSlidingPatch(bSlideL1);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput  (eTask, 0,  cinput1 );
  mgr->ConnectOutput (eTask, 1, coutput1 );

  return eTask;
}
