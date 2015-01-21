AliAnalysisTaskEmcalDiJetResponse* AddTaskEmcalDiJetResponse(TString     kTracksName         = "PicoTracks", 
							     TString     kClusName           = "caloClusterCorr",
							     TString     kMCTracksName       = "MCParticles",
							     Double_t    R                   = 0.4, 
							     Double_t    ptminTrack          = 0.15, 
							     Double_t    etminClus           = 0.3, 
							     Int_t       rhoType             = 0,
							     TString     trigClass           = "",
							     const char *CentEst             = "V0A",
							     Int_t       pSel                = AliVEvent::kINT7,
							     Int_t       matchFullCh         = AliAnalysisTaskEmcalDiJetBase::kNoMatching,
							     Double_t    ptTrackBias         = 0.,
							     Int_t       responseVar         = 0,
							     Int_t       corrType            = AliAnalysisTaskEmcalDiJetBase::kCorrelateTwo,
							     Float_t     nefCut              = 10.
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

  AliEmcalJetTask* jetFinderTaskFull    = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus);
  AliEmcalJetTask* jetFinderTaskCharged = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus);
  jetFinderTaskFull->SelectCollisionCandidates(AliVEvent::kAny);
  jetFinderTaskCharged->SelectCollisionCandidates(AliVEvent::kAny);

  AliEmcalJetTask* jetFinderTaskFullMC = AddTaskEmcalJet(kMCTracksName ,"", kANTIKT, R, kFULLJETS, kPartLevPtCut, kPartLevPtCut);
  AliEmcalJetTask* jetFinderTaskChargedMC = AddTaskEmcalJet(kMCTracksName ,"", kANTIKT, R, kCHARGEDJETS, kPartLevPtCut, kPartLevPtCut);
  jetFinderTaskFullMC->SelectCollisionCandidates(AliVEvent::kAny);
  jetFinderTaskChargedMC->SelectCollisionCandidates(AliVEvent::kAny);

  TString strJetsFull = jetFinderTaskFull->GetName();
  TString strJetsCh   = jetFinderTaskCharged->GetName();

  TString strJetsFullMC = jetFinderTaskFullMC->GetName();
  TString strJetsChMC   = jetFinderTaskChargedMC->GetName();


  TString wagonName = Form("DiJetResponse_%s_%s_Rho%dTC%sMatch%dHadTrig%dRV%d",strJetsFull.Data(),strJetsFullMC.Data(),rhoType,trigClass.Data(),matchFullCh,(Int_t)(ptTrackBias),responseVar);

  //Configure DiJet task
  AliAnalysisTaskEmcalDiJetResponse *taskDiJetResp = new AliAnalysisTaskEmcalDiJetResponse(wagonName.Data());
 
  Printf("strJetsFull: %s",strJetsFull.Data());
  Printf("strJetsCh: %s",strJetsCh.Data());

  taskDiJetResp->SetUseAliAnaUtils(kTRUE);
  taskDiJetResp->SetVzRange(-10.,10.);
  taskDiJetResp->SetIsPythia(kTRUE);

  taskDiJetResp->SetJetCorrelationType(corrType);

  taskDiJetResp->SetContainerFull(0);
  taskDiJetResp->SetContainerCharged(1);
  taskDiJetResp->SetContainerFullMC(2);
  taskDiJetResp->SetContainerChargedMC(3);

  taskDiJetResp->AddParticleContainer(kTracksName.Data());
  taskDiJetResp->AddClusterContainer(kClusName.Data());
   
  taskDiJetResp->AddJetContainer(strJetsFull.Data(),"EMCAL",R);
  taskDiJetResp->AddJetContainer(strJetsCh.Data(),"TPC",R);
  taskDiJetResp->AddJetContainer(strJetsFullMC.Data(),"EMCAL",R);
  taskDiJetResp->AddJetContainer(strJetsChMC.Data(),"TPC",R);

  taskDiJetResp->SetZLeadingCut(0.98,0.98,0);
  taskDiJetResp->SetZLeadingCut(0.98,0.98,2);

  taskDiJetResp->SetNEFCut(0.,nefCut,0);
  taskDiJetResp->SetNEFCut(0.,nefCut,2);

  for(Int_t i=0; i<4; i++) {
    taskDiJetResp->SetPercAreaCut(0.6, i);
    taskDiJetResp->SetPtBiasJetTrack(ptTrackBias,i);
  }

  taskDiJetResp->SetRhoType(rhoType);
  taskDiJetResp->SetCentralityEstimator(CentEst);
  taskDiJetResp->SelectCollisionCandidates(pSel);
  taskDiJetResp->SetFullChargedMatchingType(matchFullCh);
  taskDiJetResp->SetDoChargedCharged(kTRUE);
  taskDiJetResp->SetDoFullCharged(kTRUE);
  taskDiJetResp->SetMatchFullCharged(kTRUE);
  taskDiJetResp->SetResponseVar(responseVar);
  taskDiJetResp->SetPtMinTriggerJet(0.);

  mgr->AddTask(taskDiJetResp);

  //Connnect input
  mgr->ConnectInput (taskDiJetResp, 0, mgr->GetCommonInputContainer() );

  //Connect output
  AliAnalysisDataContainer *coutput1 = 0x0;
  TString contName(wagonName);
  contName += "_histos";
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(taskDiJetResp,1,coutput1);
  
  return taskDiJetResp;
}
