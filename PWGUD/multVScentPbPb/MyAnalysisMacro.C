void MyAnalysisMacro(TString dataset="/alice/sim/LHC10f8f_130844",
		     TString outFName="trbg.root",
		     Bool_t doRec  = kTRUE,
		     Bool_t doInj  = kTRUE,
		     Bool_t doRot  = kTRUE,
		     Bool_t doMix  = kFALSE,
		     Bool_t useMC  = kTRUE,
		     Bool_t checkReconstructables = kTRUE,
		     // 
		     Float_t etaCut     = 2.5,
		     //
		     // specific parameters for reconstruction
		     //----------------------- Zv selection parameters important for mixing, to be tuned
		     Float_t zMin       = -20,
		     Float_t zMax       =  20,
		     Float_t zMixBinSz  =  5, //0.1,
		     //---------------------------------------------------------------------------------
		     //
		     //----------------------- Ntracklets selection parameters important for mixing, to be tuned
		     Float_t ntMin      =   1,
		     Float_t ntMax      = 15000,
		     Float_t ntMixBinSz =  5000, 
		     //---------------------------------------------------------------------------------
		     //
		     float  phiRot      = 3.14159e+00,
		     float  injScale    = 0.7,
		     Bool_t scaleDTheta = kTRUE,
		     float  nStdDev     = 25.,
		     float  dphi        = 0.08, 
		     float  dtht        = 0.025, 
		     float  phishift    = 0.0045, 
		     Bool_t remOvl      = kTRUE, 
		     float  ovlPhiCut   = 0.005, 
		     float  ovlZetaCut  = 0.05,
		     Int_t  nEvents     = 50000,
		     Int_t  nEventsSkip = 0) 
{
  //  
  TString format = GetFormatFromDataSet(dataset);
  //
  // ALICE stuff
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Test train");
  //
  InputHandlerSetup(format,useMC);
  if (doMix) MixHandlerSetup(ntMin,ntMax,ntMixBinSz, zMin,zMax,zMixBinSz);
  // compile our task
  gProof->Load("AliITSMultRecBg.cxx++");
  gProof->Load("AliTrackletTaskUni.cxx++");
  //
  // load and run AddTask macro
  gROOT->LoadMacro("AddMultTaskRS.C");
  //
  // create our task
  AliTrackletTaskUni *mltTask = AddMultTaskRS(outFName.Data());
  //
  mltTask->SetDoNormalReco(doRec);
  mltTask->SetDoInjection(doInj);
  mltTask->SetDoRotation(doRot);
  mltTask->SetDoMixing(doMix);  
  //
  mltTask->SetUseMC(useMC);
  mltTask->SetCheckReconstructables(checkReconstructables);
  //
  mltTask->SetEtaCut(etaCut);
  mltTask->SetZVertexMin(zMin);
  mltTask->SetZVertexMax(zMax);
  mltTask->SetMultCutMin(ntMin);
  mltTask->SetMultCutMax(ntMax);
  //
  mltTask->SetScaleDThetaBySin2T(scaleDTheta);
  mltTask->SetNStdDev(nStdDev);
  mltTask->SetPhiWindow(dphi);
  mltTask->SetThetaWindow(dtht);
  mltTask->SetPhiShift(phishift);
  mltTask->SetPhiOverlapCut(ovlPhiCut);
  mltTask->SetZetaOverlapCut(ovlZetaCut);
  mltTask->SetPhiRot(phiRot);
  mltTask->SetInjScale(injScale);
  mltTask->SetRemoveOverlaps(remOvl);
  //
  printf("new Task: %p\n",mltTask);
  //
  AddPhysicsSelection(useMC);
  mltTask->SelectCollisionCandidates(useMC ? (AliVEvent::kMB) : (AliVEvent::kUserDefined) );
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

Bool_t InputHandlerSetup(TString format = "esd", Bool_t useKine = kTRUE)
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
      esdInputHandler = new AliESDInputHandlerRP();
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

void MixHandlerSetup(float ntMin,float ntMax,float ntMixBinSz,
		     float zMin, float zMax, float zMixBinSz)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return;
  int bufferSize = 1;
  AliESDInputHandlerRP *esdH = dynamic_cast<AliESDInputHandlerRP*>(mgr->GetInputEventHandler());
  if (!esdH) return;
  //
  AliMixEventInputHandler *esdMixH = new AliMixEventInputHandler(bufferSize);
  esdMixH->SetInputHandlerForMixing(esdH);
  AliMixEventPool *evPool = new AliMixEventPool("MyPool");
  AliMixEventCutObj *tracklets = new AliMixEventCutObj(AliMixEventCutObj::kNumberTracklets, ntMin,ntMax,ntMixBinSz);
  AliMixEventCutObj *zvertex = new AliMixEventCutObj(AliMixEventCutObj::kZVertex, zMin,zMax, zMixBinSz);
  //  evPool->AddCut(tracklets);
  evPool->AddCut(zvertex);
  //evPool->Init();
  evPool->Print();
  esdMixH->SetEventPool(evPool);
  esdH->SetMixingHandler(esdMixH);
}

void AddPhysicsSelection(Bool_t isMC)
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
