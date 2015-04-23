// $Id$

AliAnalysisTaskJetJTJT* AddTaskJetJTJT(
  const char *runPeriod 	 = "LHC13b",
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Int_t       trigger            = AliVEvent::kMB,
  Int_t       nCentBins          = 1,
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.6,
  const char *type               = "EMCAL",
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskJetJTJT",
  Int_t       effMode		 = 1,	
  Int_t       debug 		 = 0	
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetJTJT", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetJTJT", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }
    name += "_R";
    name += jetradius*10;
    name += "_T";
    name += trigger;

  TString tracksName = "PicoTracks";
  TString clustersName = "EmcCaloClusters";
  TString clustersCorrName = "CaloClustersCorr";
  TString rhoName = "";

  enum AlgoType {kKT, kANTIKT};
  enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};
  int doBkg = 1;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(tracksName,clustersCorrName,kANTIKT,0.4,kFULLJETS,0.15,0.300); // anti-kt

  if (doBkg) {
    rhoName = "Rho";
    AliEmcalJetTask* jetFinderTaskKT = AddTaskEmcalJet(tracksName, clustersCorrName, kKT, 0.4, kFULLJETS, 0.150, 0.300);

    TString kTpcKtJetsName = jetFinderTaskKT->GetName();
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");
    rhotask = (AliAnalysisTaskRho*) AddTaskRho(kTpcKtJetsName, tracksName, clustersCorrName, rhoName, 0.4, "TPC", 0.01, 0, 0, 2, kTRUE);
    //rhotask__->SetScaleFunction(sfunc);
    //rhotask->SelectCollisionCandidates(kPhysSel);
    rhotask->SetHistoBins(100,0,250);
  }

  ntracks = tracksName;
  nclusters = clustersCorrName;
  nrho = rhoName;
  Printf("name: %s",name.Data());

  AliAnalysisTaskJetJTJT* jtTask = new AliAnalysisTaskJetJTJT(name);
  jtTask->SetCentRange(0.,100.);
  jtTask->SetNCentBins(nCentBins);
  if(debug > 1){
  cout << "SetTrackArrayName: " << ntracks << endl;
  }
  jtTask->SetTrackArrayName(ntracks);
  jtTask->setDebug(debug);
  jtTask->setEffMode(effMode);
  jtTask->setRunPeriod(runPeriod);
  Double_t borders[5] = {0,10,20,40,100};
  Double_t triggpt[8] = {0,20,30,40,60,80,100,150};
  //Double_t triggpta[2] = {0,100};
  Double_t triggpta[8] = {0,0.5,1,3,10,20,50,100};
  //cout << "Size of {0,10,20,40,100}: " << borders->size() << endl;
  jtTask->setCentBinBorders(5,borders);
  jtTask->setTriggPtBorders(8,triggpt);
  jtTask->setAssocPtBorders(8,triggpta);
  jtTask->SelectCollisionCandidates(trigger);

  AliParticleContainer *trackCont  = jtTask->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = jtTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jtTask->AddJetContainer(jetFinderTask->GetName(),strType,jetradius);
  if(jetCont) {
	  jetCont->SetRhoName(nrho);
	  jetCont->ConnectParticleContainer(trackCont);
	  jetCont->ConnectClusterContainer(clusterCont);
	  jetCont->SetZLeadingCut(0.98,0.98);
	  jetCont->SetPercAreaCut(0.6);
	  jetCont->SetJetPtCut(jetptcut);    
	  jetCont->SetLeadingHadronType(leadhadtype);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jtTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  cout << "Create container " << contname << endl;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
		  TList::Class(),AliAnalysisManager::kOutputContainer,
		  Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput  (jtTask, 0,  cinput1 );
  mgr->ConnectOutput (jtTask, 1, coutput1 );

  return jtTask;
}
