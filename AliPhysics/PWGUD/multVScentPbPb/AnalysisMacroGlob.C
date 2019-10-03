void AnalysisMacroGlob(TString dataset="/alice/sim/LHC10f8f_130844",
		       TString outFName="glovar.root",
		       Int_t  nEvents     = 5000,
		       Int_t  nEventsSkip = 0) 
{
  //  
  TString format = GetFormatFromDataSet(dataset);
  //
  // ALICE stuff
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Test train");
  //
  Bool_t isMC = dataset.Contains("sim");
  InputHandlerSetup(format,kFALSE,isMC);
  gProof->Load("AliTaskGlobVar.cxx++g");
  //
  /*
  gROOT->LoadMacro("AddTaskGlobVar.C");
  AliTaskGlobVar *task = AddTaskGlobVar(outFName.Data());
  */
  AliTaskGlobVar *task = new AliTaskGlobVar("AliTaskGlobVar");
  // create output container
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("clist", TList::Class(),AliAnalysisManager::kOutputContainer,outFName.Data());
  // add our task to the manager
  mgr->AddTask(task);
  //
  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  //
  if (!isMC) {
    AddPhysicsSelection();
    task->SelectCollisionCandidates( AliVEvent::kUserDefined );
    task->SetUseMC(kFALSE);
  }
  else task->SetUseMC(kTRUE);
  //  else task->SelectCollisionCandidates( AliVEvent::kMB);
  //
  // Run analysis
  mgr->InitAnalysis();
  // process dataset  
  mgr->StartAnalysis("proof", dataset.Data(), nEvents, nEventsSkip); 
  //
  TString evstCmd = "if [ -e event_stat.root ]; then \nmv event_stat.root evstat_"; 
  evstCmd += outFName;  evstCmd += " \nfi";
  gSystem->Exec( evstCmd.Data() );
  
}


TString GetFormatFromDataSet(TString dataset) {
  
//   Info("runAAF.C","Detecting format from dataset (may take while, depends on network connection)...");
  TString dsTreeName;
  if (dataset.Contains("#")) {
    Info("runAAF.C",Form("Detecting format from dataset name '%s' ...",dataset.Data()));
    dsTreeName=dataset(dataset.Last('#'),dataset.Length());
  } else {
    Info("runAAF.C",Form("Detecting format from dataset '%s' (may take while, depends on network connection) ...",dataset.Data()));
    TFileCollection *ds = gProof->GetDataSet(dataset.Data());
    if (!ds) {
      Error(Form("Dataset %s doesn't exist on proof cluster!!!!",dataset.Data()));
      return "";
    }
    dsTreeName = ds->GetDefaultTreeName();
  }

  if (dsTreeName.Contains("esdTree")) {
    Info("runAAF.C","ESD input format detected ...");
    return "ESD";
  } else if (dsTreeName.Contains("aodTree"))  {
    Info("runAAF.C","AOD input format detected ...");
    return "AOD";
  } else {
    Error("runAAF.C",Form("Tree %s is not supported !!!",dsTreeName.Data()));
    Error("runAAF.C",Form("Maybe set your DS to %s#esdTree or %s#aodTree",dataset.Data(),dataset.Data()));
  }
  
  return "";
}

Bool_t InputHandlerSetup(TString format = "esd", Bool_t useRP=kFALSE, Bool_t useKine = kFALSE)
{
  format.ToLower();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisDataContainer *cin = mgr->GetCommonInputContainer();

  if (cin) return;

  if (!format.CompareTo("esd"))
  {
    AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdInputHandler)
    {
      Info("CustomAnalysisTaskInputSetup", "Creating esdInputHandler ...");
      if (useRP) esdInputHandler = new AliESDInputHandlerRP();
      else       esdInputHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdInputHandler);
    }
    
    if (useKine)
    {
      AliMCEventHandler* mcInputHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

      if (!mcInputHandler)
      {
        Info("CustomAnalysisTaskInputSetup", "Creating mcInputHandler ...");
        AliMCEventHandler* mcInputHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcInputHandler);
      }
    }

  }
  else if (!format.CompareTo("aod"))
  {
    AliAODInputHandler *aodInputHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!aodInputHandler)
    {
      Info("CustomAnalysisTaskInputSetup", "Creating aodInputHandler ...");
      aodInputHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodInputHandler);
    }
  }
  else
  {
    Info("Wrong input format!!! Only ESD and AOD are supported. Skipping Task ...");
    return kFALSE;
  }

  return kTRUE;
}

void AddPhysicsSelection(Bool_t isMC=kFALSE)
{
  // physics selection a la Michele
  printf("Requesting physics selection in %s mode\n",isMC ? "MC":"Data");
  gROOT->ProcessLine(".L $ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  //isMC is true when processing monte carlo, the second 0 disables the cluster vs tracklets
  AliPhysicsSelectionTask* physicsSelectionTask = AddTaskPhysicsSelection(isMC,0);
  if(!isMC) {
    AliPhysicsSelection * physSel = physicsSelectionTask->GetPhysicsSelection();
    physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");
    /*
    physSel->AddCollisionTriggerClass("+CMBS1C-B-NOPF-ALL");
    physSel->AddCollisionTriggerClass("+CMBS1A-B-NOPF-ALL");
    */
    //
    physSel->AddCollisionTriggerClass("+CMBS2C-B-NOPF-ALL");
    physSel->AddCollisionTriggerClass("+CMBS2A-B-NOPF-ALL");
    //
    // This are needed only to fill the statistics tables
    physSel->AddBGTriggerClass("+CMBAC-C-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBAC-A-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBAC-E-NOPF-ALL");
    //
    /*
    physSel->AddBGTriggerClass("+CMBS1C-C-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS1C-A-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS1C-E-NOPF-ALL");
    //
    physSel->AddBGTriggerClass("+CMBS1A-C-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS1A-A-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS1A-E-NOPF-ALL");
    //
    */
    /*
    //
    physSel->AddBGTriggerClass("+CMBS2C-C-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS2C-A-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS2C-E-NOPF-ALL");
    //
    physSel->AddBGTriggerClass("+CMBS2A-C-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS2A-A-NOPF-ALL");
    physSel->AddBGTriggerClass("+CMBS2A-E-NOPF-ALL");
    */
  } 
  // if you use the following line, your task only gets the selected events
  //  task->SelectCollisionCandidates(AliVEvent::kUserDefined);
  //
  //Alternatively, in the UserExec of your task:
  //Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kUserDefined);
  //
}
