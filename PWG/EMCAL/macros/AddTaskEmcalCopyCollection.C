AliEmcalCopyCollection* AddTaskEmcalCopyCollection(AliEmcalCorrectionTask::InputObject_t inputObjectType = AliEmcalCorrectionTask::kNoDefinedInputObject,
                          const char * collectionToCopyName = "",
                          const char * newCollectionName = "",
                          bool isEmbedding = false)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalCopyCollection", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalCopyCollection", "This task requires an input event handler");
    return 0;
  }
  
  TString name = "AliEmcalCopyCollection";
  name += TString::Format("_%s_%s", collectionToCopyName, newCollectionName);
  if (isEmbedding) {
    name += "_Embedded";
  }
  
  AliEmcalCopyCollection* mgrTask = static_cast<AliEmcalCopyCollection *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;
  
  // Create the task that manages copying collections
  AliEmcalCopyCollection* copyCollection = new AliEmcalCopyCollection(name.Data(),
                                    inputObjectType,
                                    collectionToCopyName,
                                    newCollectionName,
                                    isEmbedding);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(copyCollection);

  // Create containers for input
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(copyCollection, 0, cInput);
  
  /*TString outputContainerName(name);
  outputContainerName += "_histos";
  
  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectOutput(copyCollection, 1, cOutput);*/
  
  //TObjArray* cnt = mgr->GetContainers();

  return copyCollection;
}
