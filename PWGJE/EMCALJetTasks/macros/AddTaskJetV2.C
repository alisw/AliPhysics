/*
 * AddTask macro for class 
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands
 *
 * Note: this macro is pretty much a copy of AddTaskEmcalJetSample.C
 *
 */

AliAnalysisTaskJetV2* AddTaskJetV2(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t   jetradius          = 0.2,
  Double_t   jetptcut           = 1,
  Double_t   jetareacut         = 0.557,
  const char* type              = "TPC",
  Int_t      leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskJetV2",
  UInt_t     runMode            = AliAnalysisTaskJetV2::kGrid,
  Bool_t     fillQA             = kTRUE,
  TString    fitOpts            = "WLQI",
  UInt_t     fitType            = AliAnalysisTaskJetV2::kFourierSeries,
  TArrayD    *centralities      = 0x0,
  TRandom3   *randomizer        = 0x0,
  Double_t   trackptcut         = .15,
  Bool_t     LHC10h             = kFALSE
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

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (!strcmp(type, "TPC"))
    name += "_TPC";
  else if (!strcmp(type, "EMCAL"))
    name += "_EMCAL";
  else if (!strcmp(type, "USER")) 
    name += "_USER";

  // create instance of the object
  AliAnalysisTaskJetV2* jetTask = new AliAnalysisTaskJetV2(name, runMode);
  
  // create and connect data containers
  AliParticleContainer* partCont = jetTask->AddParticleContainer(ntracks);
  if(partCont) {
      partCont->SetName("Tracks");
      partCont->SetParticlePtCut(trackptcut);
  }
  TString tmp(nclusters);
  AliClusterContainer* clusterCont = 0x0;
  if(!tmp.IsNull()) {
      clusterCont = jetTask->AddClusterContainer(nclusters);
      jetTask->SetAnalysisType(AliAnalysisTaskJetV2::kFull);
  }
  AliJetContainer* jetCont = jetTask->AddJetContainer(njets, type, jetradius);
  if(jetCont) {
      jetCont->SetName("Jets");
      jetCont->SetPercAreaCut(jetareacut);
      jetCont->SetRhoName(nrho);
      if(partCont)      jetCont->ConnectParticleContainer(partCont);
      if(clusterCont)   jetCont->ConnectClusterContainer(clusterCont);
  }

  // task specific setters
  jetTask->SetFillQAHistograms(fillQA);
  jetTask->SetDebugMode(-1);
  jetTask->SetModulationFitType(fitType);
  jetTask->SetModulationFitOptions(fitOpts);
  jetTask->SetModulationFitMinMaxP(.01, 1);
  // if centralities haven't been specified use defaults
  if(!centralities) {
     Double_t c[] = {0., 10., 30., 50., 70., 90.};
     jetTask->SetCentralityClasses(new TArrayD(sizeof(c)/sizeof(c[0]), c));
  }
  // if a randomized hasn't specified use a safe default 
  if(!randomizer) jetTask->SetRandomSeed(new TRandom3(0));

  // pass the expected run lists to the task. the total list is used for QA plots which are stored per run-number, the semi-good list is used to change the phi acceptance of jets and pico trakcs, and - if an alternatie is provided - switch to a 'small rho' task, which also runs on limited acceptance
  Int_t totalRuns[] = {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, /* up till here original good TPC list */169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309, /* original semi-good tpc list */169415, 169411, 169035, 168988, 168984, 168826, 168777, 168512, 168511, 168467, 168464, 168342, 168310, 168115, 168108, 168107, 167987, 167915, 167903, /*new runs, good according to RCT */ 169238, 169160, 169156, 169148, 169145, 169144 /* run swith missing OROC 8 but seem ok in QA */};

  Int_t totalRuns10h[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137544, 137541, 137539, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135};

  // set the runnumbers for either 10h or 11h
  (LHC10h) ? jetTask->SetExpectedRuns(new TArrayI(sizeof(totalRuns10h)/sizeof(totalRuns10h[0]), totalRuns10h)) : jetTask->SetExpectedRuns(new TArrayI(sizeof(totalRuns)/sizeof(totalRuns[0]), totalRuns));

  Int_t semiGoodRuns[] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};

  // set the semi-good runnumbers for 10h or 11h (10h has no semi-good numbers, so pass a NULL pointer)
  (LHC10h) ? jetTask->SetExpectedSemiGoodRuns(0x0) : jetTask->SetExpectedSemiGoodRuns(new TArrayI(sizeof(semiGoodRuns)/sizeof(semiGoodRuns[0]), semiGoodRuns));

  // and if 10h, pass this info to the task so acceptance isn't changed
  if(LHC10h) jetTask->SetCollisionType(AliAnalysisTaskJetV2::kPbPb10h);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname+="_PWGJE";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );

  switch (runMode) {
      case AliAnalysisTaskJetV2::kLocal : {
          gStyle->SetOptFit(1);
          AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("good_fits_%s", name.Data()), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
          AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("bad_fits_%s", name.Data()),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s", AliAnalysisManager::GetCommonFileName()));
          mgr->ConnectOutput (jetTask, 2, coutput2);
          mgr->ConnectOutput (jetTask, 3, coutput3);
      } break;
      default: break;
  }
  return jetTask;
}
