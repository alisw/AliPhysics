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
Bool_t AddRsnEfficiency(const char *dataLabel)
{
  // load useful macros
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/ConfigESDCutsITS.C");
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/ConfigESDCutsTPC.C");
  
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  // create task
  AliRsnAnalysisEffSE *task[2];
  task[0] = new AliRsnAnalysisEffSE("RsnTaskEffNoSA");
  task[1] = new AliRsnAnalysisEffSE("RsnTaskEffSA");

  // pair definition: 
  // phi --> K+ K-
  AliRsnPairDef *pairPhi = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  
  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) rapidity
  // 3) multiplicity
  Double_t mult[] = {0., 6., 10., 15., 23., 1E6};
  Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
  AliRsnValue *axisIM   = new AliRsnValue("IM"  , AliRsnValue::kPairInvMass     , 0.9,  1.4, 0.001);
  AliRsnValue *axisPt   = new AliRsnValue("PT"  , AliRsnValue::kPairPt          , 0.0, 10.0, 0.100);
  AliRsnValue *axisY    = new AliRsnValue("Y"   , AliRsnValue::kPairY           ,-1.2,  1.2, 0.100);
  AliRsnValue *axisMult = new AliRsnValue("Mult", AliRsnValue::kEventMultESDCuts, nmult, mult);
  
  // add the support cut to the value which computes the multiplicity
  AliESDtrackCuts *cuts = new AliESDtrackCuts;
  ConfigESDCutsTPC(cuts);
  axisMult->SetSupportObject(cuts);
  
  // define cuts for event selection:
  // this will determine the filling of bins in the "info" histograms
  // and should be computed as additional correction factor in efficiency
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 1, kFALSE);
  
  // define standard 2010 track quality/PID cuts:
  // - first  index: [0] = no PID, [1] = PID
  // - second index: [0] = no ITS, [1] = ITS
  AliRsnCutESD2010 *cuts2010[2][2];
  cuts2010[0][0] = new AliRsnCutESD2010("cutESD2010nopidNoSA");
  cuts2010[0][1] = new AliRsnCutESD2010("cutESD2010nopidSA");
  cuts2010[1][0] = new AliRsnCutESD2010("cutESD2010pidNoSA");
  cuts2010[1][1] = new AliRsnCutESD2010("cutESD2010pidSA");
  // since both indexes are 0/1, the boolean settings
  // are done according to them, for clarity
  for (Int_t ipid = 0; ipid < 2; ipid++)
  {
    for (Int_t iits = 0; iits < 2; iits++)
    {
      // all work with MC here
      cuts2010[ipid][iits]->SetMC(kTRUE);
      
      // all use global tracks
      cuts2010[ipid][iits]->SetUseITSTPC(kTRUE);
      
      // other flags, depend on indexes
      cuts2010[ipid][iits]->SetUseITSSA((Bool_t)iits);
      cuts2010[ipid][iits]->SetCheckITS((Bool_t)ipid);
      cuts2010[ipid][iits]->SetCheckTPC((Bool_t)ipid);
      cuts2010[ipid][iits]->SetCheckTOF((Bool_t)ipid);
      
      // basic quality settings
      ConfigESDCutsITS(cuts2010[ipid][iits]->GetCutsITS());
      ConfigESDCutsTPC(cuts2010[ipid][iits]->GetCutsTPC());
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
    task[itask]->AddPairDef(pairPhi);

    // add the output histogram axis
    task[itask]->AddAxis(axisIM);
    task[itask]->AddAxis(axisPt);
    task[itask]->AddAxis(axisY);
    task[itask]->AddAxis(axisMult);
    
    // add the cut only when working on ESD, not on MC only
    task[itask]->GetEventCuts()->AddCut(cutVertex);
    task[itask]->GetEventCuts()->SetCutScheme(cutVertex->GetName());

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
    task[itask]->AddStepMC (mgr_step0);
    task[itask]->AddStepESD(mgr_step1);
    task[itask]->AddStepESD(mgr_step2);
    task[itask]->AddStepESD(mgr_step3);
    task[itask]->AddStepESD(mgr_step4);
    
    // add the task to the manager and connect to input
    mgr->AddTask(task[itask]);
    mgr->ConnectInput(task[itask], 0, mgr->GetCommonInputContainer());
    
    // create paths for the output in the common file
    TString infoname(task[itask]->GetName());
    TString histname(task[itask]->GetName());
    infoname.ReplaceAll("TaskEff", "Info");
    histname.ReplaceAll("TaskEff", "Hist");
    AliAnalysisDataContainer *outputInfo = mgr->CreateContainer(infoname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    AliAnalysisDataContainer *outputHist = mgr->CreateContainer(histname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
    mgr->ConnectOutput(task[itask], 1, outputInfo);
    mgr->ConnectOutput(task[itask], 2, outputHist);
  }

  return kTRUE;
}
