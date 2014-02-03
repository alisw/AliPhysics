enum AlgoType {kKT, kANTIKT};
enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

AliAnalysisTaskEmcalJetMass* AddTaskEmcalJetMass(const char * njetsBase, 
						 const char * njetsTag, 
						 const Double_t R,
						 const char * nrhoBase, 
						 const char * nrhoTag, 
						 const char * ntracks, 
						 const char * nclusters,
						 const char *type,					     
						 const char *CentEst,
						 Int_t       pSel,
						 TString     trigClass      = "",
						 TString     kEmcalTriggers = "");

AliAnalysisTaskEmcalJetMass* AddTaskEmcalJetMass(TString     kTracksName         = "PicoTracks", 
						 TString     kClusName           = "caloClustersCorr",
						 Double_t    R                   = 0.4, 
						 Double_t    ptminTrack          = 0.15, 
						 Double_t    etminClus           = 0.3, 
						 Double_t    ptminTag            = 4.,
						 Int_t       rhoType             = 1,
						 const char *type                = "EMCAL",
						 const char *CentEst             = "V0M",
						 Int_t       pSel                = AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB,
						 TString     trigClass           = "",
						 TString     kEmcalTriggers      = "",
						 TString     kPeriod             = "LHC11h"
						 ) {
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetMass","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetMass", "This task requires an input event handler");
      return NULL;
    }

  // #### Add necessary jet finder tasks
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

  AliEmcalJetTask* jetFinderTaskBase = 0x0;
  if (strcmp(type,"TPC")==0)
    jetFinderTaskBase = AddTaskEmcalJet(kTracksName, "", kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,0);
  else if (strcmp(type,"EMCAL")==0)
    jetFinderTaskBase = AddTaskEmcalJet(kTracksName, kClusName, kANTIKT, R, kFULLJETS, ptminTrack, etminClus,0.005,0);

  TString strJetsBase = jetFinderTaskBase->GetName();

  AliAnalysisTaskRhoBase *rhoTaskBase;
  TString rhoNameBase = "";
  if(rhoType==1) {
    rhoTaskBase = AttachRhoTask(kPeriod,kTracksName,kClusName,R,ptminTrack,etminClus);
    if(rhoTaskBase) {
      rhoTaskBase->SetCentralityEstimator(CentEst);  
      rhoTaskBase->SelectCollisionCandidates(AliVEvent::kAny);
      if (strcmp(type,"TPC")==0)
	rhoNameBase = rhoTaskBase->GetOutRhoName();
      if (strcmp(type,"EMCAL")==0)
	rhoNameBase = rhoTaskBase->GetOutRhoScaledName();
    }
  }

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetMass *task = AddTaskEmcalJetMass(jetFinderTaskBase->GetName(),
							  R,
							  rhoNameBase,
							  kTracksName.Data(),
							  kClusName.Data(),
							  type,
							  CentEst,
							  pSel,
							  trigClass,
							  kEmcalTriggers
							  );
  
  return task;  

}

AliAnalysisTaskEmcalJetMass* AddTaskEmcalJetMass(const char * njetsBase, 
						 const Double_t R,
						 const char * nrhoBase, 
						 const char * ntracks, 
						 const char * nclusters,
						 const char * type,
						 const char * CentEst,
						 Int_t        pSel,
						 TString      trigClass,
						 TString      kEmcalTriggers
						 ) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetMass","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetMass", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("JetMass_%s_TC%s",njetsBase,trigClass.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetMass *task = new AliAnalysisTaskEmcalJetMass(wagonName.Data());

  task->SetNCentBins(4);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  task->SetJetContainerBase(0);

  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetZLeadingCut(0.98,0.98);
    jetContBase->SetPercAreaCut(0.6);
  }

  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());

  task->SetCentralityEstimator(CentEst);

  task->SelectCollisionCandidates(pSel);

  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  

}


AliAnalysisTaskRhoBase *AttachRhoTask(TString     kPeriod             ,// = "LHC13b",
				      TString     kTracksName         ,// = "PicoTracks", 
				      TString     kClusName           ,// = "caloClustersCorr",
				      Double_t    R                   ,// = 0.4, 
				      Double_t    ptminTrack          ,// = 0.15, 
				      Double_t    etminClus            // = 0.3 
				      ) {
  
  AliAnalysisTaskRhoBase *rhoTaskBase;

  kPeriod.ToLower();

  // Add kt jet finder and rho task in case we want background subtraction
  AliEmcalJetTask *jetFinderKt;
  AliEmcalJetTask *jetFinderAKt;
  jetFinderKt   = AddTaskEmcalJet(kTracksName, "", kKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,0);
  jetFinderAKt  = AddTaskEmcalJet(kTracksName, "", kANTIKT, R, kCHARGEDJETS, ptminTrack, etminClus,0.005,0);

  if(kPeriod.EqualTo("lhc13b") || kPeriod.EqualTo("lhc13c") || kPeriod.EqualTo("lhc13d") || kPeriod.EqualTo("lhc13e") || kPeriod.EqualTo("lhc13f")) {

    jetFinderKt->SetMinJetPt(0.);

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");  
    TF1 *fScale = new TF1("fScale","1.42",0.,100.); //scale factor for pPb
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
  else if(kPeriod.EqualTo("lhc10h") || kPeriod.EqualTo("lhc11h") ) {

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");

    TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
    sfunc->SetParameter(2,1.76458);
    sfunc->SetParameter(1,-0.0111656);
    sfunc->SetParameter(0,0.000107296);
    TString rhoname = Form("RhoR%03dptmin%3.0f%s",(int)(100*R),ptminTrack*1000.0,kTracksName.Data());
    AliAnalysisTaskRho *rhoTask = AddTaskRho(
					     jetFinderKt->GetName(), 
					     kTracksName, 
					     kClusName, 
					     rhoname, 
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
