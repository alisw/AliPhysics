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
Bool_t AddRsnAnalysisTaskEffPhi
(
   const char *evtopts,
   const char *options, 
   const char *path = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train"
)
{
   // ==================================================================================================================
   // == OPTIONS =======================================================================================================
   // ==================================================================================================================
   
   // Instead or getting confused with plenty of arguments in the macro (with default values),
   // we use a unique string of options with a set of conventional strings to set up the job:
   // -- "MC"/"DATA" --> what kind of sample
   // -- "ITS"/"TPC" --> what tracks to use (ITS standalone and/or TPC+ITS)
   // -- "xxxPID"    --> add the PID cuts for the detector xxx.
   //
   // In this point, these options are converted into boolean variables.
   
   TString opt(options);
   opt.ToUpper();
   opt.ReplaceAll(" ", "");
   
   Bool_t addITS = opt.Contains("ITS");
   Bool_t addTPC = opt.Contains("TPC");
   Bool_t useITS = opt.Contains("ITSPID");
   Bool_t useTPC = opt.Contains("TPCPID");
   Bool_t useTOF = opt.Contains("TOFPID");
   
   // correct options when needed
   if (!addITS) useITS = kFALSE;
   if (!addTPC) useTPC = useTOF = kFALSE;
   
   // ==================================================================================================================
   // == DEFINITIONS ===================================================================================================
   // ==================================================================================================================
   
   // We put here the definitions of all objects which are needed in the following, in order to have then
   // a more readable code in the points where these objects are added to the analysis manager.
   
   // pair definition:
   // phi --> K+ K-
   AliRsnPairDef *pairPhi = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
   
   // ==================================================================================================================
   // == CUTS AND AXES =================================================================================================
   // ==================================================================================================================
   
   //
   // Track quality for ITS standalone:
   // this cut is used to select tracks of good quality, irrespective of the PID.
   // When adding status flags, the second argument tells if each considered flag
   // must be active or not in the track status, since the ITS-SA tracks need that
   // some of them are OFF (e.g.: kTPCin)
   //
   AliRsnCutTrackQuality *cutQualityITS = new AliRsnCutTrackQuality("cutQualityITS");
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSin    , kTRUE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kTPCin    , kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSrefit , kTRUE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kTPCrefit , kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSpureSA, kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSpid   , kTRUE);
   cutQualityITS->SetPtRange(0.15, 1E+20);
   cutQualityITS->SetEtaRange(-0.8, 0.8);
   cutQualityITS->SetDCARPtFormula("0.0595+0.0182/pt^1.55");
   cutQualityITS->SetDCAZmax(2.0);
   cutQualityITS->SetSPDminNClusters(1);
   cutQualityITS->SetITSminNClusters(4);
   cutQualityITS->SetITSmaxChi2(2.0);
   cutQualityITS->SetTPCminNClusters(0);
   cutQualityITS->SetTPCmaxChi2(1E+10);
   cutQualityITS->SetRejectKinkDaughters();
      
   //
   // Track quality for TPC+ITS:
   // works exactly like the one above, but has settings for selecting TPC+ITS tracks
   // in this case, the flags required are all necessary, so here the procedure is simpler
   //
   AliRsnCutTrackQuality *cutQualityTPC = new AliRsnCutTrackQuality("cutQualityTPC");
   cutQualityTPC->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   cutQualityTPC->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   cutQualityTPC->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);
   cutQualityTPC->SetPtRange(0.15, 1E+20);
   cutQualityTPC->SetEtaRange(-0.8, 0.8);
   cutQualityTPC->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   cutQualityTPC->SetDCAZmax(2.0);
   cutQualityTPC->SetSPDminNClusters(1);
   cutQualityTPC->SetITSminNClusters(0);
   cutQualityTPC->SetITSmaxChi2(1E+20);
   cutQualityTPC->SetTPCminNClusters(70);
   cutQualityTPC->SetTPCmaxChi2(4.0);
   cutQualityTPC->SetRejectKinkDaughters();
   
   //
   // ITS PID
   // In this implementation, it is a 3sigma cut around the Bethe-Bloch value.
   //
   // NOTE:
   // --> The initialization of the BB is different between data and MC.
   // --> The cut is the same for all momenta.
   //
   AliRsnCutPIDITS *cutPIDITSkaon = new AliRsnCutPIDITS("cutPIDITSkaon", AliPID::kKaon, -3.0, 3.0);
   cutPIDITSkaon->SetMC(kTRUE);
   cutPIDITSkaon->SetMomentumRange(0.0, 1E20);
   
   //
   // TPC PID
   // In this implementation, there are two instances:
   // - below 350 MeV --> 5 sigma cut
   // - above 350 MeV --> 3 sigma cut
   //
   // NOTE:
   // --> The initialization of the BB is different between data and MC.
   //
   AliRsnCutPIDTPC *cutPIDTPCkaonLow  = new AliRsnCutPIDTPC("cutPIDTPCkaonLow" , AliPID::kKaon, -5.0, 5.0);
   AliRsnCutPIDTPC *cutPIDTPCkaonHigh = new AliRsnCutPIDTPC("cutPIDTPCkaonHigh", AliPID::kKaon, -3.0, 3.0);
   
   // assign the momentum range and tell to reject tracks outside it
   cutPIDTPCkaonLow ->SetMomentumRange(0.00, 0.35);
   cutPIDTPCkaonHigh->SetMomentumRange(0.35, 1E20);
   cutPIDTPCkaonLow ->SetRejectOutside(kTRUE);
   cutPIDTPCkaonHigh->SetRejectOutside(kTRUE);
   
   // BB parameterization depends on data sample (MC, data)
   // the momentum range is passed and tracks outside it are rejected
   Double_t bbPar[5];
   bbPar[0] = 2.15898 / 50.0;
   bbPar[1] = 1.75295E1;
   bbPar[2] = 3.40030E-9;
   bbPar[3] = 1.96178;
   bbPar[4] = 3.91720;
   cutPIDTPCkaonLow ->SetBBParam(bbPar);
   cutPIDTPCkaonHigh->SetBBParam(bbPar);
   
   //
   // TOF PID
   // In this implementation it is a 3sigma cout aroung expected kaon time.
   //
   // NOTE:
   // --> It is important to choose if this cut must reject tracks not matched in TOF.
   //     Usually, if TPC pid is used, we can accept them, otherwise we must reject.
   //     (here we assume TPC is used)
   //
   AliRsnCutPIDTOF *cutPIDTOFkaon = new AliRsnCutPIDTOF("cutPIDTOFkaon", AliPID::kKaon, -3.0, 3.0);
   cutPIDTOFkaon->SetRejectUnmatched(!useTPC);
   
   //
   // Rapidity cut
   // Only thing to consider is that it needs a support object to define mass
   //
   AliRsnCutValue *cutRapidity = new AliRsnCutValue("cutY", AliRsnValue::kPairY, -0.5, 0.5);
   cutRapidity->GetValueObj()->SetSupportObject(pairPhi);
   
   //
   // Axes
   //
   // NOTE:
   // --> multiplicity has variable bins, defined by array below
   //
   
   Double_t mult[] = { 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13.,  14.,  15.,  16.,  17.,  18.,  19., 
                      20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 500.};
   Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
   
   AliRsnValue *axisIM      = new AliRsnValue("IM"  , AliRsnValue::kPairInvMass   ,  0.9, 1.4, 0.001);
   AliRsnValue *axisRes     = new AliRsnValue("Res" , AliRsnValue::kPairInvMassRes, -0.5, 0.5, 0.001);
   AliRsnValue *axisPt      = new AliRsnValue("PT"  , AliRsnValue::kPairPt        ,  0.0, 5.0, 0.1  );
   AliRsnValue *axisY       = new AliRsnValue("Y"   , AliRsnValue::kPairY         , -1.1, 1.1, 0.1  );
   AliRsnValue *axisMultESD = new AliRsnValue("MESD", AliRsnValue::kEventMultESDCuts, nmult, mult);
   AliRsnValue *axisMultSPD = new AliRsnValue("MSPD", AliRsnValue::kEventMultSPD    , nmult, mult);
   AliRsnValue *axisMultMC  = new AliRsnValue("MMC" , AliRsnValue::kEventMultMC     , nmult, mult);
   
   // ==================================================================================================================
   // == PRELIMINARY OPERATIONS ========================================================================================
   // ==================================================================================================================

   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   // create task
   AliRsnAnalysisTaskEffPair *task = new AliRsnAnalysisTaskEffPair(Form("RsnEff_%s", options));
   task->SelectCollisionCandidates();

   // add pair definition, to choose the checked resonance
   task->AddDef(pairPhi);

   // add the output histogram axis
   task->AddAxis(axisIM);
   task->AddAxis(axisRes);
   task->AddAxis(axisPt);
   task->AddAxis(axisY);
   task->AddAxis(axisMultSPD);
   task->AddAxis(axisMultMC);

   // ==================================================================================================================
   // == EVENT CUTS ====================================================================================================
   // ==================================================================================================================

   gROOT->LoadMacro(Form("%s/AddRsnEventComputations.C", path));
   AddRsnEventComputations(kTRUE, evtopts);

   // ==================================================================================================================
   // == STEPS =========================================================================================================
   // ==================================================================================================================
   
   Char_t qualityITS[255], qualityTPC[255];
   Char_t pidITS[255], pidTPC[255], pidTOF[255];
   Char_t schemeITS[255], schemeTPC[255], scheme[255];
   
   sprintf(qualityITS, "%s"     , cutQualityITS->GetName());
   sprintf(qualityTPC, "%s"     , cutQualityTPC->GetName());
   sprintf(pidITS    , "%s"     , cutPIDITSkaon->GetName());
   sprintf(pidTPC    , "(%s|%s)", cutPIDTPCkaonHigh->GetName(), cutPIDTPCkaonLow->GetName());
   sprintf(pidTOF    , "%s"     , cutPIDTOFkaon->GetName());
   sprintf(schemeITS , "");
   sprintf(schemeTPC , "");
   sprintf(scheme    , "");

   //
   // *** STEP 0 - All resonances which decay in the specified pair
   //
   // This step does not need any kind of definition, since
   // its requirement is automatically checked during execution,
   // but to avoid segfaults, it is needed to initialize a cut manager.
   //
   AliRsnCutManager *mgr_step0 = new AliRsnCutManager("mc_step0", "");

   //
   // *** STEP 1 - Track quality
   //
   // All resonances whose daughters were reconstructed
   // and pass quality track cuts will enter this step
   //
   AliRsnCutManager *mgr_step1 = new AliRsnCutManager("rec_step1", "");
   AliRsnCutSet     *set_step1 = mgr_step1->GetCommonDaughterCuts();

   if (addTPC && addITS) {
      set_step1->AddCut(cutQualityTPC);
      set_step1->AddCut(cutQualityITS);
      set_step1->SetCutScheme(Form("%s|%s", qualityTPC, qualityITS));
   } else if (addTPC) {
      set_step1->AddCut(cutQualityTPC);
      set_step1->SetCutScheme(qualityTPC);
   } else if (addITS) {
      set_step1->AddCut(cutQualityITS);
      set_step1->SetCutScheme(qualityITS);
   } else {
      ::Error("Need to ad at least one between ITS and TPC tracks");
      return kFALSE;
   }
   ::Info("AddRsnAnalysisTaskEffPhi", "Cut scheme for step #1: %s", set_step1->GetCutScheme().Data());

   //
   // *** STEP 2 - PID
   //
   // Add all TPC cuts, according to options
   //
   AliRsnCutManager *mgr_step2 = new AliRsnCutManager("esd_step2", "");
   AliRsnCutSet     *set_step2 = mgr_step2->GetCommonDaughterCuts();

   if (addITS && useITS) {
      sprintf(schemeITS, "%s & %s", qualityITS, pidITS);
      set_step2->AddCut(cutPIDITSkaon);
   }
   if (addTPC) {
      if (useTPC && useTOF) {
         set_step2->AddCut(cutPIDTPCkaonLow);
         set_step2->AddCut(cutPIDTPCkaonHigh);
         set_step2->AddCut(cutPIDTOFkaon);
         sprintf(schemeTPC, "%s & %s", pidTPC, pidTOF);
      } else if (useTPC) {
         set_step2->AddCut(cutPIDTPCkaonLow);
         set_step2->AddCut(cutPIDTPCkaonHigh);
         sprintf(schemeTPC, "%s & %s", pidTPC);
      } else if (useTOF) {
         set_step2->AddCut(cutPIDTOFkaon);
         sprintf(schemeTPC, "%s & %s", pidTOF);
      }
   }
      
   // final scheme depends on what of the above were added
   // in case both ITS-SA and TPC tracks are added, the scheme
   // is the OR of the cuts for the first and those for the second
   // category of tracks
   if (strlen(schemeITS) > 0 && strlen(schemeTPC) > 0) {
      sprintf(scheme, "(%s) | (%s)", schemeITS, schemeTPC);
   } else if (strlen(schemeITS) > 0) {
      sprintf(scheme, "%s", schemeITS);
   } else if (strlen(schemeTPC) > 0) {
      sprintf(scheme, "%s", schemeTPC);
   } else {
      ::Error("Scheme is empty!");
      return kFALSE;
   }
   // check scheme
   set_step2->SetCutScheme(scheme);
   ::Info("AddRsnAnalysisTaskEffPhi", "Cut scheme for step #1: %s", set_step2->GetCutScheme().Data());

   // add all steps to the task:
   // - first step computed on MC
   // - all other steps computed on reconstruction
   task->AddStepMC(mgr_step0);
   task->AddStepRec(mgr_step1);
   task->AddStepRec(mgr_step2);

   // add the task to the manager and connect to input
   mgr->AddTask(task);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

   // create paths for the output in the common file
   TString infoname(task->GetName());
   TString histname(task->GetName());
   infoname.Append("_Info");
   histname.Append("_Hist");
   AliAnalysisDataContainer *outputInfo = mgr->CreateContainer(infoname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
   AliAnalysisDataContainer *outputHist = mgr->CreateContainer(histname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
   mgr->ConnectOutput(task, 1, outputInfo);
   mgr->ConnectOutput(task, 2, outputHist);

   return kTRUE;
}
