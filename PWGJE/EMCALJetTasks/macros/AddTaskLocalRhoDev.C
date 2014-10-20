/*
 * AddTask macro for class 
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands
 */

AliAnalysisTaskLocalRhoDev* AddTaskLocalRhoDev(
  const char *ntracks           = "Tracks",     // track selection used for vn estimate
  const char *nclusters         = "",           // clusters (not used) 
  const char *njets             = "Jets",       // jet selection for finding leading jet
  const char *nrho              = "Rho",        // name of nominal rho
  const char *lrho              = "LocalRho",   // name of local rho object
  Double_t   jetradius          = 0.3,          // jet radius (to remove leading jet)
  Double_t   jetptcut           = 1,            
  Double_t   jetareacut         = 0.557,
  const char* type              = "TPC",
  Int_t      leadhadtype        = 0,
  const char *name              = "AliAnalysisTaskLocalRhoDev",    // task name
  TString    fitOpts            = "WLQI",                       // options for tfitter
  UInt_t     fitType            = AliAnalysisTaskLocalRhoDev::kCombined,   // fitting strategy
  TArrayI    *centralities      = 0x0,                                  // centrality binning for qa
  UInt_t     runMode            = AliAnalysisTaskLocalRhoDev::kGrid        // run mode 
  )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskLocalRhoDev* jetTask = new AliAnalysisTaskLocalRhoDev(name, runMode);
  // inherited setters
  AliParticleContainer* partCont = jetTask->AddParticleContainer(ntracks);
  if(partCont) {
      partCont->SetName("Tracks");
      partCont->SetParticlePtCut(.15);
  }
  AliJetContainer* jetCont = jetTask->AddJetContainer(njets, type, jetradius);
  if(jetCont) {
      jetCont->SetName("Jets");
      jetCont->SetPercAreaCut(jetareacut);
      jetCont->SetRhoName(nrho);
      jetCont->ConnectParticleContainer(partCont);
  }
  // task specific setters
  jetTask->SetDebugMode(-1);
  jetTask->SetModulationFitType(fitType);
  jetTask->SetModulationFitOptions(fitOpts);
  jetTask->SetModulationFitMinMaxP(.01, 1);
  // if centralities haven't been specified use defaults
  if(!centralities) {
     Double_t c[] = {0., 10., 30., 50., 70., 90.};
     jetTask->SetCentralityClasses(new TArrayD(sizeof(c)/sizeof(c[0]), c));
  }
  // pass the expected run lists to the task. the total list is used for QA plots which are stored per run-number, the semi-good list is used to change the phi acceptance of jets and pico trakcs, and - if an alternatie is provided - switch to a 'small rho' task, which also runs on limited acceptance
  Int_t totalRuns[] = {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, /* up till here original good TPC list */169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309, /* original semi-good tpc list */169415, 169411, 169035, 168988, 168984, 168826, 168777, 168512, 168511, 168467, 168464, 168342, 168310, 168115, 168108, 168107, 167987, 167915, 167903, /*new runs, good according to RCT */ 169238, 169160, 169156, 169148, 169145, 169144 /* run swith missing OROC 8 but seem ok in QA */};
  jetTask->SetExpectedRuns(new TArrayI(sizeof(totalRuns)/sizeof(totalRuns[0]), totalRuns));

  Int_t semiGoodRuns[] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};
  jetTask->SetExpectedSemiGoodRuns(new TArrayI(sizeof(semiGoodRuns)/sizeof(semiGoodRuns[0]), semiGoodRuns));


  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );

  switch (runMode) {
      case AliAnalysisTaskLocalRhoDev::kLocal : {
          gStyle->SetOptFit(1);
          AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("good_fits_%s", name), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
          AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("bad_fits_%s", name),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s", AliAnalysisManager::GetCommonFileName()));
          mgr->ConnectOutput (jetTask, 2, coutput2);
          mgr->ConnectOutput (jetTask, 3, coutput3);
      } break;
      default: break;
  }
  return jetTask;
}

AliAnalysisTaskLocalRhoDev* AddTaskLocalRhoDev(
  const char *ntracks           = "Tracks",     // track selection used for vn estimate
  const char *nclusters         = "",           // clusters (not used) 
  const char *njets             = "Jets",       // jet selection for finding leading jet
  const char *nrho              = "Rho",        // name of nominal rho
  const char *lrho              = "LocalRho",   // name of local rho object
  Double_t   jetradius          = 0.3,          // jet radius (to remove leading jet)
  Double_t   jetptcut           = 1,            
  Double_t   jetareacut         = 0.557,
  const Int_t type              = 0,
  Int_t      leadhadtype        = 0,
  const char *name              = "AliAnalysisTaskLocalRhoDev",    // task name
  TString    fitOpts            = "WLQI",                       // options for tfitter
  UInt_t     fitType            = AliAnalysisTaskLocalRhoDev::kCombined,   // fitting strategy
  TArrayI    *centralities      = 0x0,                                  // centrality binning for qa
  UInt_t     runMode            = AliAnalysisTaskLocalRhoDev::kGrid        // run mode 
  )
{
  // overloaded addtask for backward compatibility
  UInt_t jetType = 0;
  jetType |= AliEmcalJetTask::kAKT;

  if (type==0)
    jetType |= AliEmcalJetTask::kFullJet; 
  else if (type==1) 
    jetType |= AliEmcalJetTask::kChargedJet; 
  else if (type==2) 
    jetType |= AliEmcalJetTask::kNeutralJet;

  if (jetradius==0.2) 
    jetType |= AliEmcalJetTask::kR020Jet; 
  else if (jetradius==0.3) 
    jetType |= AliEmcalJetTask::kR030Jet;
  else if (jetradius==0.4) 
    jetType |= AliEmcalJetTask::kR040Jet;
  else
    jetType |= AliEmcalJetTask::kRX1Jet;

  return AddTaskLocalRhoDev(ntracks, nclusters, njets, nrho, lrho, jetradius, jetptcut, jetareacut, jetType, leadhadtype, name, fitOpts, fitType, centralities, runMode);
}
