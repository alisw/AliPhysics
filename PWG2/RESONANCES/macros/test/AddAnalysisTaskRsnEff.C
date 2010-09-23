//
// This macro add an analysis task for computing efficiency.
// It will have as output an AliCFContainer with several steps:
//
//  0) all resonances in MC which decay in the pair specified
//  1) all resonances in ESD (i.e.: reconstruction/acceptance effect)
//  2) subset of (1) whose daughters satisfy quality track cuts (covariance, chi square && nTPCclusters, DCA)
//  3) subset of (2) whose daughters satisty the PID cuts
//
Bool_t AddAnalysisTaskRsnEff(const char *dataLabel)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  // create task
  AliRsnAnalysisEffSE *task = new AliRsnAnalysisEffSE("phiEff");

  // pair definitions:
  // phi --> K+ K-
  AliRsnPairDef *pairPhi = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  task->AddPairDef(pairPhi);

  // axis definition
  //  0) transverse momentum (fixed bins of 100 MeV/c)
  //  1) rapidity (variable bins)
  //  2) invariant mass (fixed bins of 1 MeV/c^2)
  Double_t y[] = {-0.8, -0.7, -0.6, -0.5, 0.5, 0.6, 0.7, 0.8};
  Int_t    ny  = sizeof(y) / sizeof(y[0]);
  AliRsnValue *axisIM = new AliRsnValue("IM", AliRsnValue::kPairInvMass, 500, 0.9,  1.4);
  AliRsnValue *axisPt = new AliRsnValue("PT", AliRsnValue::kPairPt     , 100, 0.0, 10.0);
  AliRsnValue *axisY  = new AliRsnValue("Y" , AliRsnValue::kPairY      ,  ny, y);
  task->AddAxis(axisIM);
  task->AddAxis(axisPt);
  task->AddAxis(axisY);

  // define cuts for event selection:
  // this will determine the filling of bins in the "info" histograms
  // and should be computed as additional correction factor in efficiency
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  task->GetEventCuts()->AddCut(cutVertex);
  task->GetEventCuts()->SetCutScheme("cutVertex");

  //
  // *** STEP 0 - All resonances which decay in the specified pairs
  //
  // This step does not need any kind of definition, since
  // its requirement is automatically checked during execution,
  // but a cut manager needs to be defined at each step, even if empty.
  //
  AliRsnCutManager *mgr_step0 = new AliRsnCutManager("mc_step0", "");
  
  // add this as a step on MonteCarlo
  task->AddStepMC (mgr_step0);
  
  //
  // *** STEP 1 - All resonances which decay into tracked particles
  //
  // This step does not need any kind of definition, since
  // its requirement is automatically checked during execution,
  // but a cut manager needs to be defined at each step, even if empty.
  //
  AliRsnCutManager *mgr_step1 = new AliRsnCutManager("mc_step0", "");
  
  // add this as a step on reconstructed tracks
  task->AddStepESD(mgr_step1);
  
  //
  // *** STEP 2 - Track quality
  //
  
  // track cut -----------------------
  // --> global cuts for 2010 analysis
  AliRsnCutESD2010 *cuts2010 = new AliRsnCutESD2010("cutESD2010");
  // ----> set the flag for sim management
  cuts2010->SetMC(kTRUE);
  // ----> disable checking of PID, but include also ITS stand-alone
  cuts2010->SetCheckITS (kFALSE);
  cuts2010->SetCheckTPC (kFALSE);
  cuts2010->SetCheckTOF (kFALSE);
  cuts2010->SetUseGlobal(kTRUE);
  cuts2010->SetUseITSSA (kTRUE);
  // ----> set defaults for the rest
  cuts2010->InitializeToDefaults(kTRUE);
  
  // cut sets ---------------------------------
  // --> only common cuts for tracks are needed
  // --> standard 2010 cuts are applied only when working on ESD
  AliRsnCutSet *cutSetDaughterCommon_step2 = new AliRsnCutSet("commonDaughterCuts_step2", AliRsnCut::kDaughter);
  cutSetDaughterCommon_step2->AddCut(cuts2010);
  cutSetDaughterCommon_step2->SetCutScheme(cuts2010->GetName());
  AliRsnCutManager *mgr_step2   = new AliRsnCutManager("esd_step2", "");
  mgr_step2->SetCommonDaughterCuts(cutSetDaughterCommon_step2);
  
  // add this as a step on reconstructed tracks
  task->AddStepESD(mgr_step2);
  
  //
  // *** STEP 3 - PID
  //
  
  // track cut -----------------------
  // --> global cuts for 2010 analysis
  AliRsnCutESD2010 *cuts2010_pid = new AliRsnCutESD2010("cutESD2010_pid");
  // ----> set the flag for sim/data management
  cuts2010_pid->SetMC(kTRUE);
  // ----> enable PID checking
  cuts2010_pid->SetCheckITS (kTRUE);
  cuts2010_pid->SetCheckTPC (kTRUE);
  cuts2010_pid->SetCheckTOF (kTRUE);
  cuts2010_pid->SetUseGlobal(kTRUE);
  cuts2010_pid->SetUseITSSA (kTRUE);
  // ----> set defaults for the rest
  cuts2010_pid->InitializeToDefaults(kTRUE);
  
  // cut sets ---------------------------------
  // --> only common cuts for tracks are needed
  // --> standard 2010 cuts are applied only when working on ESD
  AliRsnCutSet *cutSetDaughterCommon_step3 = new AliRsnCutSet("commonDaughterCuts_step3", AliRsnCut::kDaughter);
  cutSetDaughterCommon_step3->AddCut(cuts2010_pid);
  cutSetDaughterCommon_step3->SetCutScheme(cuts2010_pid->GetName());
  AliRsnCutManager *mgr_step3   = new AliRsnCutManager("esd_step3", "");
  mgr_step3->SetCommonDaughterCuts(cutSetDaughterCommon_step3);

  // add this as a step on reconstructed tracks
  task->AddStepESD(mgr_step3);

  // add task to manager and setup containers
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // create paths for the output in the common file
  Char_t commonPath[500];
  sprintf(commonPath, "%s", AliAnalysisManager::GetCommonFileName());

  // create paths for the output in the common file
  AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("RsnInfoEff", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
  AliAnalysisDataContainer *outputHist = mgr->CreateContainer("RsnEff" , TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
  mgr->ConnectOutput(task, 1, outputInfo);
  mgr->ConnectOutput(task, 2, outputHist);

  return kTRUE;
}
