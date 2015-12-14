enum AlgoType {kKT, kANTIKT};
enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

AliAnalysisTaskEmcalTriggerPatchJetMatch* AddTaskEmcalTriggerPatchJetMatch(
  TString     strJets1            = "jets1",
  TString     strJets2            = "jets2",
  TString     kTracksName         = "PicoTracks",
  TString     kClusName           = "caloClustersCorr",
  Double_t    R                   = 0.4,
  Double_t    ptminTrack          = 0.15,
  Double_t    etminClus           = 0.3,
  Int_t       rhoType             = 0,
  TString     trigClass           = "",
  TString     kEmcalCellsName     = "",
  const char *CentEst             = "V0A",
  Int_t       pSel                = AliVEvent::kINT7,
  Float_t     nefCut              = 10.,
  TString     kEmcalTriggers      = "",
  TString     kPeriod             = "LHC13b",
  TString     kBeamType           = "pp", //or pPb or PbPb
  Bool_t      comments            = kFALSE,
  Bool_t      UseAllRecalcPatches = kFALSE,
  Bool_t      newFramework        = 0,
  TString     tag                 = ""
) {
  // The following three lines are added for backwards compatibility
  kPeriod.ToLower();
  if(kPeriod.EqualTo("lhc10h") || kPeriod.EqualTo("lhc11h")) kBeamType = "PbPb";
  if(kPeriod.EqualTo("lhc13b") || kPeriod.EqualTo("lhc13c") || kPeriod.EqualTo("lhc13d") || kPeriod.EqualTo("lhc13e") || kPeriod.EqualTo("lhc13f")) kBeamType = "pPb";

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEmcalTriggerPatchJetMatch","No analysis manager found.");
    return 0;
  }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEmcalTriggerPatchJetMatch", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskRhoBase *rhoTask;
  if(rhoType==1) {
    rhoTask = AttachRhoTask(kBeamType,kTracksName,kClusName,R,ptminTrack,etminClus);
    if(rhoTask) {
      rhoTask->SetCentralityEstimator(CentEst);  
      rhoTask->SelectCollisionCandidates(AliVEvent::kAny);
    }
    else {
      Warning("AddTaskEmcalTriggerPatchJetMatch","Asked for rho task but configuration unknown. Continuing configuration without rho task.");
      rhoType = 0;
    }
  }

  TString wagonName = Form("TriggerQA_%s_TC%s_%s",strJets1.Data(),trigClass.Data(),tag.Data());

  //Configure TriggerQA task
  AliAnalysisTaskEmcalTriggerPatchJetMatch *task = new AliAnalysisTaskEmcalTriggerPatchJetMatch(wagonName);

  AliParticleContainer *trackCont  = task->AddParticleContainer(kTracksName.Data());
  AliClusterContainer *clusterCont = task->AddClusterContainer(kClusName.Data());

  task->SetJetContainer(0);
  AliJetContainer *jetCont0 = task->AddJetContainer(strJets1.Data(),"EMCAL",R);
//  AliClusterContainer *JETconstit = task->AddClusterContainer("EmcCaloClusters");
  AliClusterContainer *JETconstit = task->AddClusterContainer(kClusName);
  jetCont0->ConnectClusterContainer(JETconstit);

  if(rhoType==1) task->SetRhoName(rhoTask->GetOutRhoScaledName(),0);
//  task->SetZLeadingCut(0.98,0.98,0);
//  task->SetNEFCut(0.,nefCut,0);
  task->SetPercAreaCut(0.6, 0);

  task->SetTriggerClass(trigClass.Data());
  task->SetCaloCellsName(kEmcalCellsName.Data());
  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());

  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);

  task->SetUseAliAnaUtils(kTRUE);
  task->SetVzRange(-10.,10.);
  if(kPeriod.Contains("LHC13b4")) task->SetIsPythia(kTRUE);
  task->SetdoComments(comments);
  task->SetUseALLrecalcPatches(UseAllRecalcPatches);
  if(newFramework) task->SetCaloClustersName("caloClusters");

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
    AliAnalysisTaskRhoSparse *rhoTaskSparse = AddTaskRhoSparse(jetFinderKt->GetName(), jetFinderAKt->GetName(), kTracksName, kClusName, Form("RhoSparseR%03d",(int)(100*R)), R, "TPC", 0.01, 0.15, 0, fScale, 0, kTRUE, Form("RhoSparseR%03d",(int)(100*R)), kTRUE);
    rhoTaskSparse->SetUseAliAnaUtils(kTRUE);
    rhoTaskBase = dynamic_cast<AliAnalysisTaskRhoBase*>rhoTaskSparse;
  }
  else if(kBeamType=="PbPb") {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");

    TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
    sfunc->SetParameter(2,1.76458);
    sfunc->SetParameter(1,-0.0111656);
    sfunc->SetParameter(0,0.000107296);
    AliAnalysisTaskRho *rhoTask = AddTaskRho(jetFinderKt->GetName(), kTracksName, kClusName, Form("RhoR%03d",(int)(100*R)), R, "TPC", 0.01, 0, sfunc, 2, kTRUE);
    rhoTask->SetHistoBins(100,0,250);

    rhoTaskBase = dynamic_cast<AliAnalysisTaskRhoBase*>rhoTask;
  }

  return rhoTaskBase;
}
