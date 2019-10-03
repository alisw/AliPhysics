AliAnalysisTaskEmcalDiJetAna* AddTaskEmcalDiJetAna(TString     kTracksName         = "PicoTracks", 
						   TString     kClusName           = "caloClusterCorr",
						   Double_t    R                   = 0.4, 
						   Double_t    ptminTrack          = 0.15, 
						   Double_t    etminClus           = 0.3, 
						   Int_t       rhoType             = 0,
						   TString     trigClass           = "",
						   const char *CentEst             = "V0A",
						   Int_t       pSel                = AliVEvent::kINT7,
						   Int_t       matchFullCh         = AliAnalysisTaskEmcalDiJetBase::kNoMatching,
						   Double_t    ptTrackBias         = 0.,
						   Int_t       corrType            = AliAnalysisTaskEmcalDiJetBase::kCorrelateTwo,
						   Float_t     nefCut              = 10.,
						   Int_t       nCentBins           = 5,
						   Double_t    scaleFact           = 1.28
						   ) {
  
  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

  // #### Define manager and data container names
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalDiJet", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalDiJet", "This task requires an input event handler");
      return NULL;
    }

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

  AliEmcalJetTask* jetFinderTaskFull    = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
  AliEmcalJetTask* jetFinderTaskCharged = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
  jetFinderTaskFull->SelectCollisionCandidates(AliVEvent::kAny);
  jetFinderTaskCharged->SelectCollisionCandidates(AliVEvent::kAny);


  TString strJetsFull = jetFinderTaskFull->GetName();
  TString strJetsCh   = jetFinderTaskCharged->GetName();

  // Add kt jet finder and rho task in case we want background subtraction
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
  AliEmcalJetTask *jetFinderKt;
  AliEmcalJetTask *jetFinderAKt;
  AliAnalysisTaskRhoSparse *rhoTask;
  if(rhoType==1) {
    jetFinderKt   = AddTaskEmcalJet(kTracksName, kClusName, kKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",0.);
    jetFinderAKt  = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
    jetFinderKt->SelectCollisionCandidates(AliVEvent::kAny);
    jetFinderAKt->SelectCollisionCandidates(AliVEvent::kAny);

    TF1 *fScale = new TF1("fScale","[0]",0.,100.);
    fScale->SetParameter(0,scaleFact);
    TString rhoSparseName = Form("RhoSparseR%03d",(int)(100*R));
    rhoTask = AddTaskRhoSparse(jetFinderKt->GetName(),
			       jetFinderAKt->GetName(),
			       kTracksName.Data(),
			       kClusName.Data(),
			       rhoSparseName.Data(),
			       R,
			       "TPC",
			       0.01,
			       0.15,
			       0,
			       fScale,
			       0,
			       kTRUE,
			       rhoSparseName.Data(),
			       kTRUE
			       );
    rhoTask->SetCentralityEstimator(CentEst);
 
  }
  TString wagonName = Form("DiJet_%s_%s_Rho%dTC%sMatch%dHadTrig%d",strJetsFull.Data(),strJetsCh.Data(),rhoType,trigClass.Data(),matchFullCh,(Int_t)(ptTrackBias));

  //Configure DiJet task
  AliAnalysisTaskEmcalDiJetAna *taskDiJet = NULL;
  taskDiJet = new AliAnalysisTaskEmcalDiJetAna(wagonName.Data());

  taskDiJet->SetUseAliAnaUtils(kTRUE);
  taskDiJet->SetVzRange(-10.,10.);
  taskDiJet->SetTriggerClass(trigClass.Data());

  if(ptminTrack==0.) {
    taskDiJet->SetIsPythia(kTRUE);
    taskDiJet->SetDoFullFull(kTRUE);
  }
  taskDiJet->SetJetCorrelationType(corrType);


  Printf("strJetsFull: %s",strJetsFull.Data());
  Printf("strJetsCh: %s",strJetsCh.Data());

  taskDiJet->AddParticleContainer(kTracksName.Data());
  taskDiJet->AddClusterContainer(kClusName.Data());
   
  taskDiJet->SetContainerFull(0);
  taskDiJet->SetContainerCharged(1);
  taskDiJet->AddJetContainer(strJetsFull.Data(),"EMCAL",R);
  taskDiJet->AddJetContainer(strJetsCh.Data(),"TPC",R);

  taskDiJet->SetZLeadingCut(0.98,0.98,0);
  taskDiJet->SetNEFCut(0.,nefCut,0);

  for(Int_t i=0; i<2; i++) {
    taskDiJet->SetPercAreaCut(0.6, i);
    taskDiJet->SetPtBiasJetTrack(ptTrackBias,i);
  }

  taskDiJet->SetRhoType(rhoType);
  if(rhoType==1) {
    taskDiJet->SetRhoName(rhoTask->GetOutRhoScaledName(),0);
    taskDiJet->SetRhoName(rhoTask->GetOutRhoName(),1);
  }

  taskDiJet->SetCentralityEstimator(CentEst);
  taskDiJet->SetCentRange(0.,100.);
  taskDiJet->SetNCentBins(nCentBins);

  taskDiJet->SelectCollisionCandidates(pSel);

  taskDiJet->SetFullChargedMatchingType(matchFullCh);

  taskDiJet->SetDoChargedCharged(kTRUE);
  taskDiJet->SetDoFullCharged(kTRUE);
  taskDiJet->SetMatchFullCharged(kFALSE);

  mgr->AddTask(taskDiJet);

  //Connnect input
  mgr->ConnectInput (taskDiJet, 0, mgr->GetCommonInputContainer() );

  //Connect output
  AliAnalysisDataContainer *coutput1 = 0x0;

  TString contName(wagonName);
  contName += "_histos";

  //  TString outputfile = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),wagonName.Data());
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());

  coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  mgr->ConnectOutput(taskDiJet,1,coutput1);
  
  return taskDiJet;
}
