void MyAnalysisMacroTrackletMulti
(TString dataset="/alice/sim/LHC10f8c_130844",
 TString outFName = "trbg.root",
 Int_t   nEvents   = -1,
 Float_t etaMin     =-0.5,        // min eta range to fill in histos
 Float_t etaMax     = 0.5,        // max eta range to fill in histos
 Float_t zMin       = -7,         // process events with Z vertex min
 Float_t zMax       =  7,         //                     max positions
 Int_t   useCentVar = 0,          // centrality variable to use: 0=V0, 1=SPD2corr
 //
 Float_t cutSigNStd  = 1.5,        // cut on weighed distance used to extract signal
 Float_t cutSigDPhiS = -1,        // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
 Bool_t  useMC  = kTRUE,          // fill MC info (doRec=kTRUE)
 Float_t scaleMCV0 = 0.7520,      // rescale MC V0 to match the data
 //
 Bool_t doRec  = kTRUE,           // fill data histos from new reco
 Bool_t doInj  = kTRUE,           // create Inj. bg
 Bool_t doRot  = kFALSE,          // create Rot. bg
 Bool_t doMix  = kFALSE,//kTRUE,  // create Mix. bg
 // 
 // specific parameters for reconstruction
 float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
 float  injScale    = 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
 Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
 float  nStdDev     = 25.,         // number of st.dev. for tracklet cut to keep
 float  dphi        = 0.06,        // dphi window (sigma of tracklet cut)
 float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
 float  phishift    = 0.0045,      // bending shift
 Bool_t remOvl      = kTRUE,       
 float  ovlPhiCut   = 0.005, 
 float  ovlZetaCut  = 0.05,
 Int_t  nEventsSkip = 0,
 //----------------------- Ntracklets selection parameters important for mixing, to be tuned
 Float_t ntMin      =   1,         // process events with ESDmult 
 Float_t ntMax      = 20000,       // within this range
 Float_t ntMixBinSz = 20000,       // ESDMult bin size for mixing
 //----------------------- Zv selection parameters important for mixing, to be tuned
 Float_t zMixBinSz  =  14,       //0.1,  // Zv. bin for mixing
 //---------------------------------------------------------------------------------
 //
 Bool_t checkReconstructables = kFALSE//kTRUE, // fill histos for reconstructable (needs useMC and doRec) 
 //
 )
{
  //  
  if (cutSigDPhiS<0) cutSigDPhiS = TMath::Sqrt(cutSigNStd)*dphi;
  //
  printf("Start Analysis for %s, max %d Events skipping %d, Event Cuts: %.1f<eta<%.1f, %.2f<Zv<%.2f\n",
	 dataset.Data(),nEvents,nEventsSkip,etaMin,etaMax,zMin,zMax);
  printf("Centrality variable: %d\n",useCentVar);
  printf("Tracklet cuts: dPhi:%.3f dTheta:%.3f phiShift:%.4f | Keep %.1f NstDev\n"
	 "Scale dTheta: %s | Signal Selection: NstDev:%.1f, dPhiS: %.3f\n", 
	 dphi,dtht,phishift,nStdDev,scaleDTheta ? "ON":"OFF",
	 cutSigNStd,cutSigDPhiS);
  //
  printf("UseMC: %s. V0 scale: %.4f\n",useMC ? "ON":"OFF",scaleMCV0);
  printf("Operations:  \n"
	 "Reco:%s (RemOvl:%s phi:%.3f zeta:%.3f)\n"
	 "Inj:%s (scale: %.2f)\n"
	 "Rot:%s (phi: %.4f)\n"
	 "Mix:%s (Nt:%.1f:%.1f:%.1f Zv:%.1f:%.1f:%.2f)\n", 
	 doRec ? "ON":"OFF",remOvl? "ON":"OFF",ovlPhiCut,ovlZetaCut,
	 doInj ? "ON":"OFF",injScale,
	 doRot ? "ON":"OFF",phiRot,
	 doMix ? "ON":"OFF",ntMin,ntMax,ntMixBinSz, zMin,zMax,zMixBinSz);
  //
  if (nEvents<0) nEvents = int(1e9);
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
  gProof->Load("AliTrackletTaskMulti.cxx++");
  //
  printf("Loading Centrality task\n");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  taskCentrality->SetDebugLevel(2);
  if (useMC) taskCentrality->SetMCInput();
  //  taskCentrality->Dump();
  //
  // load and run AddTask macro
  gROOT->LoadMacro("AddMultTaskTrackletMulti.C");
  //
  // create our task
  AliTrackletTaskMulti *mltTask = AddMultTaskTrackletMulti(outFName.Data());
  //
  mltTask->SetUseCentralityVar(useCentVar);
  mltTask->SetDoNormalReco(doRec);
  mltTask->SetDoInjection(doInj);
  mltTask->SetDoRotation(doRot);
  mltTask->SetDoMixing(doMix);  
  //
  mltTask->SetUseMC(useMC);
  mltTask->SetCheckReconstructables(checkReconstructables);
  //
  mltTask->SetEtaMin(etaMin);
  mltTask->SetEtaMax(etaMax);
  mltTask->SetZVertexMin(zMin);
  mltTask->SetZVertexMax(zMax);
  //
  mltTask->SetDPhiSCut(cutSigDPhiS);
  mltTask->SetNStdCut(cutSigNStd);
  mltTask->SetScaleMCV0(scaleMCV0);
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
  printf("Requesting physics selection in %s mode\n",useMC ? "MC":"Data");
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelectionTask = AddTaskPhysicsSelection(useMC,0);
  mltTask->SelectCollisionCandidates();//AliVEvent::kMB);
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
  if(!isMC) {
    //AliPhysicsSelection * physSel = physicsSelectionTask->GetPhysicsSelection();
    //    physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");
    /*
    physSel->AddCollisionTriggerClass("+CMBS1C-B-NOPF-ALL");
    physSel->AddCollisionTriggerClass("+CMBS1A-B-NOPF-ALL");
    */
    //
    //    physSel->AddCollisionTriggerClass("+CMBS2C-B-NOPF-ALL");
    //    physSel->AddCollisionTriggerClass("+CMBS2A-B-NOPF-ALL");
    //
    // This are needed only to fill the statistics tables
    //    physSel->AddBGTriggerClass("+CMBAC-C-NOPF-ALL");
    //    physSel->AddBGTriggerClass("+CMBAC-A-NOPF-ALL");
    //    physSel->AddBGTriggerClass("+CMBAC-E-NOPF-ALL");
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
  //  task->SelectCollisionCandidates();
  //
  //Alternatively, in the UserExec of your task:
  //Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kUserDefined);
  //
}
