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
   
   // for cuts and axes, load the support macro
   //gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PWG2/RESONANCES");
   //gROOT->LoadMacro(Form("%s/CPhiCutsAndAxes.C", path));
   
   AliRsnCutTrackQuality *cutQualityITS     = TrackQualityITS();
   AliRsnCutTrackQuality *cutQualityTPC     = TrackQualityTPC();
   AliRsnCutPIDITS       *cutPIDITSkaon     = PIDITS(kTRUE, AliRsnDaughter::kKaon, 0.0, 1E20, 3.0, "cutPIDITSkaon");
   AliRsnCutPIDTPC       *cutPIDTPCkaonLow  = PIDTPC(kTRUE, AliRsnDaughter::kKaon, 0.0, 0.350, 5.0, "cutPIDTPCkaonLow");
   AliRsnCutPIDTPC       *cutPIDTPCkaonHigh = PIDTPC(kTRUE, AliRsnDaughter::kKaon, 0.350, 1E20, 3.0, "cutPIDTPCkaonHigh");
   AliRsnCutPIDTOF       *cutPIDTOFkaon     = PIDTOF(AliRsnDaughter::kKaon, 3.0, !useTPC, "cutPIDTOFkaon");
   
   AliRsnValue *axisIM      = AxisIM();
   AliRsnValue *axisPt      = AxisPt();
   AliRsnValue *axisRes     = AxisRes();
   AliRsnValue *axisMultSPD = AxisMultSPD();
   AliRsnValue *axisMultMC  = AxisMultMC();
   AliRsnValue *axisY       = AxisY();
   
   // ==================================================================================================================
   // == PRELIMINARY OPERATIONS ========================================================================================
   // ==================================================================================================================

   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   // create task
   AliRsnAnalysisEffSE *task = new AliRsnAnalysisEffSE(Form("RsnEff_%s", options));
   task->SelectCollisionCandidates();

   // add pair definition, to choose the checked resonance
   task->AddPairDef(pairPhi);

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
   task->AddStepESD(mgr_step1);
   task->AddStepESD(mgr_step2);

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
