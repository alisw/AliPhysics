//
// This macro add an analysis task for computing efficiency.
// It will have as output an AliCFContainer with several steps:
//
//  0) all resonances in MC which decay in the pair specified
//  1) subset of (0) whose daughters are in acceptance
//  2) subset of (1) whose daughters satisfy quality track cuts (covariance, chi square && nTPCclusters)
//  3) subset of (2) whose daughters satisfy primary track cuts (nsigma to vertex, no kink daughters)
//  4) subset of (3) whose daughters satisty the BB TPC compatibility cut at 3 sigma
//
Bool_t AddRsnEfficiency(const char *dataLabel, const char *path)
{
  // load useful macros
  gROOT->LoadMacro(Form("%s/QualityCutsITS.C", path));
  gROOT->LoadMacro(Form("%s/QualityCutsTPC.C", path));
  
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  // create task
  AliRsnAnalysisEffSE *task1[2], *task2[2];
  task1[0] = new AliRsnAnalysisEffSE("RsnTaskEffNoSA_1");
  task1[1] = new AliRsnAnalysisEffSE("RsnTaskEffSA_1");
  task2[0] = new AliRsnAnalysisEffSE("RsnTaskEffNoSA_2");
  task2[1] = new AliRsnAnalysisEffSE("RsnTaskEffSA_2");
  task1[0]->SelectCollisionCandidates();
  task1[1]->SelectCollisionCandidates();
  task2[0]->SelectCollisionCandidates();
  task2[1]->SelectCollisionCandidates();

  // pair definition: 
  // kstar --> K+ pi- or vice versa
  AliRsnPairDef *pairKstar1 = new AliRsnPairDef(AliPID::kPion, '+', AliPID::kKaon, '-', 313, 0.896);
  AliRsnPairDef *pairKstar1 = new AliRsnPairDef(AliPID::kPion, '-', AliPID::kKaon, '+', 313, 0.896);
  
  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) rapidity
  // 3) multiplicity
  Double_t     mult[]   = {0., 6., 10., 15., 23., 100000000.0};
  Int_t        nmult    = sizeof(mult) / sizeof(mult[0]);
  AliRsnValue *axisIM   = new AliRsnValue("IM"  , AliRsnValue::kPairInvMass     , 0.9,  1.4, 0.001);
  AliRsnValue *axisPt   = new AliRsnValue("PT"  , AliRsnValue::kPairPt          , 0.0, 10.0, 0.100);
  AliRsnValue *axisY    = new AliRsnValue("Y"   , AliRsnValue::kPairY           ,-1.2,  1.2, 0.100);
  AliRsnValue *axisMult = new AliRsnValue("Mult", AliRsnValue::kEventMultESDCuts, nmult, mult);
  
  // initialize the support object: AliESDtrackCuts
  // configured using the standard values
  AliESDtrackCuts *cuts = new AliESDtrackCuts(QualityCutsTPC());
  axisMult->SetSupportObject(cuts);
  
  // define cuts for event selection:
  // this will determine the filling of bins in the "info" histograms
  // and should be computed as additional correction factor in efficiency
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  cutVertex->SetCheckPileUp(kTRUE);
  
  // define standard 2010 track quality/PID cuts:
  // - first  index: [0] = no PID, [1] = PID
  // - second index: [0] = no ITS, [1] = ITS
  AliRsnCutESD2010 *cuts2010[2][2];
  cuts2010[0][0] = new AliRsnCutESD2010("cutESD2010nopidNoSA");
  cuts2010[0][1] = new AliRsnCutESD2010("cutESD2010nopidSA");
  cuts2010[1][0] = new AliRsnCutESD2010("cutESD2010pidNoSA");
  cuts2010[1][1] = new AliRsnCutESD2010("cutESD2010pidSA");
  // define Bethe-Bloch parameters (only for MC, since this computes efficiency)
  Double_t bbPar[5] = {2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720};
  // since both indexes are 0/1, the boolean settings are done according to them, for clarity
  for (Int_t ipid = 0; ipid < 2; ipid++)
  {
    for (Int_t iits = 0; iits < 2; iits++)
    {
      // all work with MC here
      cuts2010[ipid][iits]->SetMC(kTRUE);
      
      // PID reference is kaons
      cuts2010[ipid][iits]->SetPID(AliPID::kKaon);
      
      // all use global tracks
      cuts2010[ipid][iits]->SetUseITSTPC(kTRUE);
      
      // other flags, depend on indexes
      cuts2010[ipid][iits]->SetUseITSSA((Bool_t)iits);
      cuts2010[ipid][iits]->SetCheckITS((Bool_t)ipid);
      cuts2010[ipid][iits]->SetCheckTPC((Bool_t)ipid);
      cuts2010[ipid][iits]->SetCheckTOF((Bool_t)ipid);
      
      // basic quality settings
      cuts2010[ipid][iits]->CopyCutsTPC(QualityCutsTPC());
      cuts2010[ipid][iits]->CopyCutsITS(QualityCutsITS());
      
      // (unused for No PID) set the ITS PID-related variables
      cuts2010[ipid][iits]->SetITSband(3.0);
      
      // (unused for No PID) set the TPC PID-related variables
      cuts2010[ipid][iits]->SetTPCrange(3.0, 5.0);
      cuts2010[ipid][iits]->SetTPCpLimit(0.35);
      cuts2010[ipid][iits]->GetESDpid()->GetTPCResponse().SetBetheBlochParameters(bbPar[0], bbPar[1], bbPar[2], bbPar[3], bbPar[4]);
      
      // (unused for No PID) set the TOF PID-related variables
      cuts2010[ipid][iits]->SetTOFrange(-3.0, 3.0);
    }
  }

  // define cut on dip angle:
  AliRsnCutValue *cutDip = new AliRsnCutValue("cutDip", AliRsnValue::kPairDipAngle, 0.02, 1.01);
  
  // define a common path for the output file
  Char_t commonPath[500];
  sprintf(commonPath, "%s", AliAnalysisManager::GetCommonFileName());
  
  // add all the steps
  // two-folded loop on the two tasks, where one contains the ITS-SA and the other doesn't
  for (Int_t itask = 0; itask < 2; itask++)
  {
    // add pair definition, to choose the checked resonance
    task1[itask]->AddPairDef(pairKstar1);
    task2[itask]->AddPairDef(pairKstar2);

    // add the output histogram axis
    task1[itask]->AddAxis(axisIM);
    task1[itask]->AddAxis(axisPt);
    task1[itask]->AddAxis(axisY);
    task1[itask]->AddAxis(axisMult);
    task2[itask]->AddAxis(axisIM);
    task2[itask]->AddAxis(axisPt);
    task2[itask]->AddAxis(axisY);
    task2[itask]->AddAxis(axisMult);
    
    // add the cut on primary vertex
    task1[itask]->GetEventCuts()->AddCut(cutVertex);
    task1[itask]->GetEventCuts()->SetCutScheme(cutVertex->GetName());
    task2[itask]->GetEventCuts()->AddCut(cutVertex);
    task2[itask]->GetEventCuts()->SetCutScheme(cutVertex->GetName());

    //
    // *** STEP 0 - All resonances which decay in the specified pair
    //
    // This step does not need any kind of definition, since
    // its requirement is automatically checked during execution,
    // but to avoid segfaults, it is better to initialize a cut manager.
    //
    AliRsnCutManager *mgr_step0 = new AliRsnCutManager("mc_step0", "");
    
    //
    // *** STEP 1 - All resonances which decay into tracked particles
    //
    // This step does not need any kind of definition, since
    // its requirement is automatically checked during execution,
    // but to avoid segfaults, it is better to initialize a cut manager.
    //
    AliRsnCutManager *mgr_step1 = new AliRsnCutManager("esd_step0", "");

    //
    // *** STEP 2 - Reconstruction & track quality
    //
    // Define a cut on track quality, disabling the PID cuts (first index = [0])
    //
    AliRsnCutManager *mgr_step2 = new AliRsnCutManager("esd_step2", "");
    AliRsnCutSet     *set_step2 = mgr_step2->GetCommonDaughterCuts();
    
    set_step2->AddCut(cuts2010[0][itask]);
    set_step2->SetCutScheme(cuts2010[0][itask]->GetName());
    
    //
    // *** STEP 3 - PID
    //
    // Define a cut on track quality, enabling the PID cuts (first index = [1])
    //
    AliRsnCutManager *mgr_step3 = new AliRsnCutManager("esd_step3", "");
    AliRsnCutSet     *set_step3 = mgr_step3->GetCommonDaughterCuts();
    
    set_step3->AddCut(cuts2010[1][itask]);
    set_step3->SetCutScheme(cuts2010[1][itask]->GetName());

    //
    // *** STEP 4 - Dip angle
    //
    // Add a cut on the pair dip angle
    //
    AliRsnCutManager *mgr_step4 = new AliRsnCutManager("esd_step4", "");
    AliRsnCutSet     *set_step4 = mgr_step4->GetMotherCuts();
    
    set_step4->AddCut(cutDip);
    set_step4->SetCutScheme(Form("%s", cutDip->GetName()));
    
    // add all steps to the task:
    // - first step computed on MC
    // - all other steps computed on reconstruction
    task1[itask]->AddStepMC (mgr_step0);
    task1[itask]->AddStepESD(mgr_step1);
    task1[itask]->AddStepESD(mgr_step2);
    task1[itask]->AddStepESD(mgr_step3);
    task1[itask]->AddStepESD(mgr_step4);
    task2[itask]->AddStepMC (mgr_step0);
    task2[itask]->AddStepESD(mgr_step1);
    task2[itask]->AddStepESD(mgr_step2);
    task2[itask]->AddStepESD(mgr_step3);
    task2[itask]->AddStepESD(mgr_step4);
    
    // add the task to the manager and connect to input
    mgr->AddTask(task1[itask]);
    mgr->AddTask(task2[itask]);
    mgr->ConnectInput(task1[itask], 0, mgr->GetCommonInputContainer());
    mgr->ConnectInput(task2[itask], 0, mgr->GetCommonInputContainer());
    
    // create paths for the output in the common file
    TString infoname1(task1[itask]->GetName()), infoname2(task2[itask]->GetName());
    TString histname1(task1[itask]->GetName()), histname2(task2[itask]->GetName());
    infoname1.ReplaceAll("TaskEff", "Info");
    histname1.ReplaceAll("TaskEff", "Hist");
    infoname2.ReplaceAll("TaskEff", "Info");
    histname2.ReplaceAll("TaskEff", "Hist");
    AliAnalysisDataContainer *outputInfo1 = mgr->CreateContainer(infoname1.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    AliAnalysisDataContainer *outputHist1 = mgr->CreateContainer(histname1.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    AliAnalysisDataContainer *outputInfo2 = mgr->CreateContainer(infoname2.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    AliAnalysisDataContainer *outputHist2 = mgr->CreateContainer(histname2.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    mgr->ConnectOutput(task1[itask], 1, outputInfo1);
    mgr->ConnectOutput(task1[itask], 2, outputHist1);
    mgr->ConnectOutput(task2[itask], 1, outputInfo2);
    mgr->ConnectOutput(task2[itask], 2, outputHist2);
  }

  return kTRUE;
}
