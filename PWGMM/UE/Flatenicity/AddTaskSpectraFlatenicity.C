///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskSpectraFlatenicity macro                    //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskSpectraFlatenicity*
AddTaskSpectraFlatenicity(  const Char_t *taskname = "Flat", 
                            Bool_t woTrivialscaling = kFALSE, 
                            Bool_t useMC = kTRUE,
                            Double_t minpT = 0.5,
                            Int_t tracksyst = 0
                         )
{
  // get the manager via the static access member. since it's static, you don't
  // need an instance of the class to call the function

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and
  // feeds events to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  // now you create an instance of your task
  AliAnalysisTaskSpectraFlatenicity *taskFlat = new AliAnalysisTaskSpectraFlatenicity("taskFlat");

  if (!taskFlat)
    return 0x0;
  
  taskFlat->SetUseMC(useMC);
  taskFlat->SetPtMin(minpT);
  taskFlat->SetRemoveTrivialScaling(woTrivialscaling);
  taskFlat->SetSysVarTrkCuts(tracksyst);
  
  mgr->AddTask(taskFlat);

  const Char_t *complement;
  if (woTrivialscaling) {
    complement = "wotrivialscal";
  } else {
    complement = "wtrivialscal";
  }
  
  mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput( taskFlat, 1, mgr->CreateContainer( Form("cList%s_%s", taskname, complement), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));
  return taskFlat;
}
