enum AlgoType {kKT, kANTIKT};
enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

JETriggerRejectionAna::AliAnalysisTaskTriggerRejection* AddTaskTriggerRejection(TString     kTracksName         = "PicoTracks", 
							  TString     kClusName           = "caloClustersCorr",
							  Double_t    R                   = 0.4, 
							  Double_t    ptminTrack          = 0.15, 
							  Double_t    etminClus           = 0.3, 
							  Int_t       rhoType             = 1,
							  const char *CentEst             = "V0M",
							  Int_t       pSel                = AliVEvent::kAny,
							  TString     kEmcalTriggers      = "",
							  TString     kPeriod             = "LHC11h",
							  TString     kBeamType           = "PbPb", //or pPb or pp
							  TString     tag                 = ""
							   ) {
  // The following three lines are added for backwards compatibility
  kPeriod.ToLower();
  if(kPeriod.EqualTo("lhc10h") || kPeriod.EqualTo("lhc11h")) kBeamType = "PbPb";
  if(kPeriod.EqualTo("lhc13b") || kPeriod.EqualTo("lhc13c") || kPeriod.EqualTo("lhc13d") || kPeriod.EqualTo("lhc13e") || kPeriod.EqualTo("lhc13f")) kBeamType = "pPb";

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskTriggerRejection","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskTriggerRejection", "This task requires an input event handler");
      return NULL;
    }

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

  AliEmcalJetTask* jetFinderTask1;
  AliEmcalJetTask* jetFinderTask2;
  if(kClusName.IsNull()) {  //particle level jets
    jetFinderTask1 = AddTaskEmcalJet(kTracksName, "", kANTIKT, R, kFULLJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
    jetFinderTask2 = AddTaskEmcalJet(kTracksName, "", kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
  }
  else if(kTracksName.IsNull()) { //neutral/calo jets
    jetFinderTask1 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kNEUTRALJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
    jetFinderTask2 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
  }
  else { //full jets
    jetFinderTask1 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
    jetFinderTask2 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
  }
  jetFinderTask1->SelectCollisionCandidates(AliVEvent::kAny);
  jetFinderTask2->SelectCollisionCandidates(AliVEvent::kAny);

  TString strJets1 = jetFinderTask1->GetName();
  TString strJets2 = jetFinderTask2->GetName();

  AliAnalysisTaskRhoBase *rhoTask;
  if(rhoType==1) {
    rhoTask = AttachRhoTask(kBeamType,kTracksName,kClusName,R,ptminTrack,etminClus);
    if(rhoTask) {
      rhoTask->SetCentralityEstimator(CentEst);  
      rhoTask->SelectCollisionCandidates(AliVEvent::kAny);
    }
    else {
      Warning("AddTaskTriggerRejection","Asked for rho task but configuration unknown. Continuing configuration without rho task.");
      rhoType = 0;
    }
  }

  TString wagonName = Form("TriggerRejectionQA_%s_%s_%s",strJets1.Data(),strJets2.Data(),tag.Data());

  //Configure TriggerQA task
  JETriggerRejectionAna::AliAnalysisTaskTriggerRejection *task = new JETriggerRejectionAna::AliAnalysisTaskTriggerRejection(wagonName);
  AliParticleContainer *trackCont  = task->AddParticleContainer(kTracksName.Data());
  AliClusterContainer *clusterCont = task->AddClusterContainer(kClusName.Data());

  task->SetContainerFull(0);
  task->SetContainerCharged(1);
  AliJetContainer *jetCont0 = task->AddJetContainer(strJets1.Data(),"EMCAL",R);
  if(rhoType==1) task->SetRhoName(rhoTask->GetOutRhoScaledName(),0);
  AliJetContainer *jetCont1 = NULL;
  if(strJets2.Contains("Charged")) {
    jetCont1 = task->AddJetContainer(strJets2.Data(),"TPC",R);
    if(rhoType==1) task->SetRhoName(rhoTask->GetOutRhoName(),1);
  }
  else { 
    jetCont1 = task->AddJetContainer(strJets2.Data(),"EMCAL",R);
    task->SetZLeadingCut(0.98,0.98,1);
  }

  task->SetZLeadingCut(0.98,0.98,0);
  for(Int_t i=0; i<2; i++)
    task->SetPercAreaCut(0.6, i);

  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);

  task->SetUseAliAnaUtils(kFALSE);
  task->SetVzRange(-10.,10.);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  AliAnalysisDataContainer *coutput1 = 0x0;
  TString containerName1 = Form("%s",wagonName.Data());
  TString outputfile = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),wagonName.Data());
  coutput1 = mgr->CreateContainer(containerName1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  
}

AliAnalysisTaskRhoBase *AttachRhoTask(TString     kBeamType           = "pp",
				      TString     kTracksName         = "PicoTracks", 
				      TString     kClusName           = "caloClustersCorr",
				      Double_t    R                   = 0.4, 
				      Double_t    ptminTrack          = 0.15, 
				      Double_t    etminClus           = 0.3 
				      ) {
  
  AliAnalysisTaskRhoBase *rhoTaskBase;

  // Add kt jet finder and rho task in case we want background subtraction
  Double_t minJetPt = 0.1;
  if(kBeamType == "pPb") minJetPt = 0.;
  AliEmcalJetTask *jetFinderKt;
  AliEmcalJetTask *jetFinderAKt;
  jetFinderKt   = AddTaskEmcalJet(kTracksName, kClusName, kKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",minJetPt);
  jetFinderAKt  = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,1,"Jet",1.);
  jetFinderKt->SelectCollisionCandidates(AliVEvent::kAny);
  jetFinderAKt->SelectCollisionCandidates(AliVEvent::kAny);

  if(kBeamType == "pPb") {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");  
    TF1 *fScale = new TF1("fScale","1.28",0.,100.); //scale factor for pPb
    AliAnalysisTaskRhoSparse *rhoTaskSparse = AddTaskRhoSparse(
			       jetFinderKt->GetName(),
			       jetFinderAKt->GetName(),
			       kTracksName,
			       kClusName,
			       Form("RhoSparseR%03d",(int)(100*R)),
			       R,
			       "TPC",
			       0.01,
			       0.15,
			       0,
			       fScale,
			       0,
			       kTRUE,
			       Form("RhoSparseR%03d",(int)(100*R)),
			       kTRUE
			       );
    rhoTaskSparse->SetUseAliAnaUtils(kTRUE);
    rhoTaskBase = dynamic_cast<AliAnalysisTaskRhoBase*>rhoTaskSparse;
  }
  else if(kBeamType=="PbPb") {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");

    TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
    sfunc->SetParameter(2,1.76458);
    sfunc->SetParameter(1,-0.0111656);
    sfunc->SetParameter(0,0.000107296);
    AliAnalysisTaskRho *rhoTask = AddTaskRho(
					     jetFinderKt->GetName(), 
					     kTracksName, 
					     kClusName, 
					     Form("RhoR%03d",(int)(100*R)), 
					     R, 
					     "TPC", 
					     0.01, 
					     0, 
					     sfunc, 
					     2, 
					     kTRUE);
    rhoTask->SetHistoBins(100,0,250);

    rhoTaskBase = dynamic_cast<AliAnalysisTaskRhoBase*>rhoTask;
  }

  return rhoTaskBase;
}
