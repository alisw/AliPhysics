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
  if (inputDataType != "AOD") {
    ::Info("AddTaskEmcalAodTpcTrack", "This task is only needed for AOD analysis. No task added.");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalAodTrackFilterTask *aodTask = new AliEmcalAodTrackFilterTask(taskName);
  aodTask->SetTracksOutName(name);
  aodTask->SetTracksInName(inname);
  aodTask->SetMC(kFALSE);

  Bool_t includeNoITS  = kFALSE;
  Bool_t doProp        = kFALSE; //force propagation of all tracks to EMCal
  Bool_t doAttemptProp = kTRUE;  //only propagate the tracks which were not propagated during AOD filtering
  Bool_t isMC          = kFALSE;

  TString runPeriod(runperiod);
  runPeriod.ToLower();
  if (runPeriod == "lhc10b" || runPeriod == "lhc10c" || runPeriod == "lhc10d" ||
      runPeriod == "lhc10e" || runPeriod == "lhc10h" || 
      runPeriod == "lhc11h" || runPeriod == "lhc12a" || runPeriod == "lhc12b" ||
      runPeriod == "lhc12c" || runPeriod == "lhc12d" || runPeriod == "lhc12e" ||
      runPeriod == "lhc12f" || runPeriod == "lhc12g" || runPeriod == "lhc12h" ||
      runPeriod == "lhc12i" || runPeriod == "lhc13b" || runPeriod == "lhc13c" ||
      runPeriod == "lhc13d" || runPeriod == "lhc13e" || runPeriod == "lhc13f" ||
      runPeriod == "lhc13g"
      ) {
    aodTask->SetAODfilterBits(256,512); // hybrid tracks
    if (runPeriod == "lhc10b" || runPeriod == "lhc10c" || runPeriod == "lhc10d" || runPeriod == "lhc10e" || runPeriod == "lhc10h") {
      includeNoITS = kTRUE;
    }
  } else if (runPeriod == "lhc10f7a"    || runPeriod == "lhc12a15e"   || runPeriod.Contains("lhc12a17") || runPeriod == "lhc13b4" ||
	     runPeriod == "lhc13b4_fix" || runPeriod == "lhc13b4_plus"    || runPeriod.Contains("lhc14a1") || runPeriod.Contains("lhc13b2_efix")
	     ) {
    aodTask->SetAODfilterBits(256,512); // hybrid tracks
    isMC = kTRUE;
    if (runPeriod == "lhc10f7a") {
      includeNoITS = kTRUE;
    }
  } else if (runPeriod == "lhc11a" || runPeriod == "lhc10hold") {
    aodTask->SetAODfilterBits(256,16); // hybrid tracks
    includeNoITS = kTRUE;
  } else if(runPeriod == "lhc11d") {
    aodTask->SetAODfilterBits(256,16); // hybrid tracks (MV: not 100% sure)
    includeNoITS = kFALSE;
  } else if (runPeriod.Contains("lhc12a15a") || runPeriod == "lhc12a15f" || runPeriod == "lhc12a15g") {
    aodTask->SetAODfilterBits(256,16); // hybrid tracks
    isMC = kTRUE;
    includeNoITS = kTRUE;
  } else if (runPeriod.Contains("lhc11a1")){
    aodTask->SetAODfilterBits(256, 16);
    isMC = kTRUE;
    includeNoITS=kTRUE;
  }  else if (runPeriod.Contains(":")) {
    TString runPeriodToken(runperiod);
    TObjArray *arr = runPeriodToken.Tokenize(":");
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
      if (arg3.Contains("doAttemptProp=kFALSE"))
	doAttemptProp=kFALSE;
      if (arg3.Contains("isMC=kTRUE"))
	isMC = kTRUE;
    }
    aodTask->SetAODfilterBits(arg1.Atoi(),arg2.Atoi());
    delete arr;
  } else {
    if (!runPeriod.IsNull())
      ::Warning("AddTaskEmcalAodTrackFilter", Form("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data()));
  }
  aodTask->SetIncludeNoITS(includeNoITS);
  aodTask->SetDoPropagation(doProp);
  aodTask->SetAttemptProp(doAttemptProp);
  aodTask->SetMC(isMC);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(aodTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(aodTask, 0,  cinput1 );
  
  return aodTask;
}
