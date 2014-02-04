// $Id$

AliEmcalAodTrackFilterTask* AddTaskEmcalAodTrackFilter(
  const char *name         = "FilterTracks",
  const char *inname       = "tracks",
  const char *runperiod    = "",
  const char *taskName     = "AliEmcalAodTrackFilterTask"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskAodTrackFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskAodTrackFilter", "This task requires an input event handler");
    return NULL;
  }

  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (inputDataType != "AOD")) {
    ::Error("AddTaskAodTrackFilter", "This task works only on AOD analysis");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalAodTrackFilterTask *aodTask = new AliEmcalAodTrackFilterTask(taskName);
  aodTask->SetTracksOutName(name);
  aodTask->SetTracksInName(inname);
  aodTask->SetMC(kFALSE);

  Bool_t includeNoITS = kFALSE;
  Bool_t doProp = kFALSE;
  TString runPeriod(runperiod);
  runPeriod.ToLower();
  if (runPeriod == "lhc11h" || runPeriod == "lhc13b" || runPeriod == "lhc13c" || 
      runPeriod == "lhc13d" || runPeriod == "lhc13e" || runPeriod == "lhc13f" || 
      runPeriod == "lhc13g" || runPeriod == "lhc12g" || runPeriod == "lhc10h" ||
      runPeriod == "lhc10d" || runPeriod == "lhc10e" || runPeriod == "lhc12d") {
    aodTask->SetAODfilterBits(256,512); // hybrid tracks
    if (runPeriod == "lhc10h" || runPeriod == "lhc10d" || runPeriod == "lhc10e")
      includeNoITS = kTRUE;
  } else if (runPeriod == "lhc12a15e" || runPeriod == "lhc13b4" || runPeriod == "lhc13b4_fix" || runPeriod == "lhc12a15f") {
    aodTask->SetAODfilterBits(256,512); // hybrid tracks
    aodTask->SetMC(kTRUE);
  } else if (runPeriod == "lhc11a" || runPeriod == "lhc10hold") {
    aodTask->SetAODfilterBits(256,16); // hybrid tracks
    includeNoITS = kTRUE;
  } else if (runPeriod.Contains("lhc12a15a")) {
    aodTask->SetAODfilterBits(256,16); // hybrid tracks
    aodTask->SetMC(kTRUE);
    includeNoITS = kTRUE;
  } else if (runPeriod.Contains(":")) {
    TObjArray *arr = runPeriod.Tokenize(":");
    TString arg1(arr->At(0)->GetName());
    TString arg2("-1");
    if (arr->GetEntries()>1)
      arg2 = arr->At(1)->GetName();
    if (arr->GetEntries()>2) {
      TString arg3 = arr->At(2)->GetName();
      if (arg3.Contains("includeNoITS=kTRUE"))
	includeNoITS=kTRUE;
      if (arg3.Contains("doProp=kTRUE"))
	doProp=kTRUE;
    }
    aodTask->SetAODfilterBits(arg1.Atoi(),arg2.Atoi());
    delete arr;
  } else {
    if (!runPeriod.IsNull())
      ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data());
  }
  aodTask->SetIncludeNoITS(includeNoITS);
  aodTask->SetDoPropagation(doProp);
  //aodTask->SetAttemptProp(kTRUE);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(aodTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(aodTask, 0,  cinput1 );
  
  return aodTask;
}
