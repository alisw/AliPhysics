AliAnalysisTaskEmcalJetTriggerQA* AddTaskEmcalJetTriggerQA(TString     kTracksName         = "PicoTracks", 
							   TString     kClusName           = "caloClusterCorr",
							   Double_t    R                   = 0.4, 
							   Double_t    ptminTrack          = 0.15, 
							   Double_t    etminClus           = 0.3, 
							   Int_t       rhoType             = 0,
							   UInt_t      type                = AliAnalysisTaskEmcal::kEMCAL,
							   TString     trigClass           = "",
							   TString     kEmcalCellsName     = "",
							   const char *CentEst             = "V0A",
							   Int_t       pSel                = AliVEvent::kINT7
							   ) {

  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetTriggerQA","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetTriggerQA", "This task requires an input event handler");
      return NULL;
    }

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

  AliEmcalJetTask* jetFinderTask1;
  AliEmcalJetTask* jetFinderTask2;
  if(kClusName.IsNull()) {  //particle level jets
    jetFinderTask1 = AddTaskEmcalJet(kTracksName, "", kANTIKT, R, kFULLJETS, ptminTrack, etminClus);
    jetFinderTask2 = AddTaskEmcalJet(kTracksName, "", kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus);
  }
  else if(kTracksName.IsNull()) { //neutral/calo jets
    jetFinderTask1 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kNEUTRALJETS, ptminTrack, etminClus);
    jetFinderTask2 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus);
  }
  else { //full jets
    jetFinderTask1 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus);
    jetFinderTask2 = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus);
  }

  TString strJets1 = jetFinderTask1->GetName();
  TString strJets2 = jetFinderTask2->GetName();

  // Add kt jet finder and rho task in case we want background subtraction
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
  AliEmcalJetTask *jetFinderKt;
  AliEmcalJetTask *jetFinderAKt;
  AliAnalysisTaskRhoSparse *rhoTask;
  if(rhoType==1) {
    jetFinderKt   = AddTaskEmcalJet(kTracksName, kClusName, kKT, R, kCHARGEDJETS, ptminTrack, etminClus);
    jetFinderKt->SetMinJetPt(0.);
    jetFinderAKt  = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus);
    TF1 *fScale = new TF1("fScale","1.42",0.,100.);
    rhoTask = AddTaskRhoSparse(
			       jetFinderKt->GetName(),
			       jetFinderAKt->GetName(),
			       kTracksName,
			       kClusName,
			       Form("RhoSparseR%03d",(int)(100*R)),
			       R,
			       AliAnalysisTaskEmcal::kTPC,
			       0.01,
			       0.15,
			       0,
			       fScale,
			       0,
			       kTRUE,
			       Form("RhoSparseR%03d",(int)(100*R)),
			       kTRUE
			       );
    rhoTask->SetCentralityEstimator(CentEst);
    
  }

  TString wagonName = Form("TriggerQA_%s_%s_TC%s",strJets1.Data(),strJets2.Data(),trigClass.Data());

  //Configure TriggerQA task
  AliAnalysisTaskEmcalJetTriggerQA *task = new AliAnalysisTaskEmcalJetTriggerQA(wagonName);
  task->SetTracksName(kTracksName.Data());
  task->SetClusName(kClusName.Data());
  task->SetJetsName(strJets1.Data());
  task->SetJetRadius(R);
  task->SetJetPtCut(0.15);
  task->SetPercAreaCut(0.6);
  task->SetTrackPtCut(ptminTrack);
  task->SetClusPtCut(etminClus);
  task->SetAnaType(type);
  task->SetTriggerClass(trigClass.Data());
  task->SetCaloCellsName(kEmcalCellsName.Data());

  task->SetJetsName2(strJets2.Data());
  if(strJets2.Contains("Charged")) {
    task->SetMinEtaJets2(-0.9+R);
    task->SetMaxEtaJets2(0.9-R);
    task->SetMinPhiJets2(-10.);
    task->SetMaxPhiJets2(10.);
  }
  else {
    task->SetMinEtaJets2(-0.7+R);
    task->SetMaxEtaJets2(0.7-R);
    task->SetMinPhiJets2(1.4+R);
    task->SetMaxPhiJets2(TMath::Pi()-R);
  }

  if(rhoType==1) {
    task->SetRhoName(rhoTask->GetRhoScaledName());
    task->SetRhoChName(rhoTask->GetRhoName());
  }
  task->SetCentralityEstimator(CentEst);

  task->SelectCollisionCandidates(pSel);

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
