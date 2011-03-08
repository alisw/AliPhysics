//
// *** Configuration script for phi->KK analysis with 2010 runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t RsnConfigPhi
(
   Bool_t      isMC,
   const char *options,
   const char *path     = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train",
   const char *taskName = "RSNtask"
)
{
   // ==================================================================================================================
   // == PRELIMINARY OPERATIONS ========================================================================================
   // ==================================================================================================================

   // retrieve task from manager, using its name
   AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
   AliRsnAnalysisTask *task = (AliRsnAnalysisTask*)mgr->GetTask(taskName);
   if (!task) {
      Error("RsnConfigMonitor", "Task not found");
      return kFALSE;
   }
   
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
   
   // pair definitions --> decay channels:
   // in our case, unlike-charged KK pairs for the signal, and like-charged ones for background
   AliRsnPairDef *pairDef[3];
   pairDef[0] = new AliRsnPairDef(AliRsnDaughter::kKaon, '+', AliRsnDaughter::kKaon, '-', 333, 1.019455); // unlike
   pairDef[1] = new AliRsnPairDef(AliRsnDaughter::kKaon, '+', AliRsnDaughter::kKaon, '+', 333, 1.019455); // like ++
   pairDef[2] = new AliRsnPairDef(AliRsnDaughter::kKaon, '-', AliRsnDaughter::kKaon, '-', 333, 1.019455); // like --

   // computation objects:
   // these are the objects which drive the computations, whatever it is (histos or tree filling)
   // and all tracks passed to them will be given the mass according to the reference pair definition
   // we create two unlike-sign pair computators, one for all pairs and another for true pairs (useful in MC)
   AliRsnPairFunctions *pair[4];
   pair[0] = new AliRsnPairFunctions(Form("%s_kaonP_kaonM_phi", opt.Data()), pairDef[0]); // unlike - true
   pair[1] = new AliRsnPairFunctions(Form("%s_kaonP_kaonM_all", opt.Data()), pairDef[0]); // unlike - all
   pair[2] = new AliRsnPairFunctions(Form("%s_kaonP_kaonM_all", opt.Data()), pairDef[1]); // like ++
   pair[3] = new AliRsnPairFunctions(Form("%s_kaonP_kaonM_all", opt.Data()), pairDef[2]); // like --

   // set additional option for true pairs (slot [0])
   pair[0]->SetOnlyTrue(kTRUE);
   pair[0]->SetCheckDecay(kTRUE);
   
   // for cuts and axes, load the support macro
   gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PWG2/RESONANCES");
   gROOT->LoadMacro(Form("%s/CPhiCutsAndAxes.C++", path));
   
   // ==================================================================================================================
   // == SINGLE DAUGHTER CUTS ==========================================================================================
   // ==================================================================================================================
   
   // we need to define cuts for track quality and for PID
   // which one will be used, depend on the analysis type
   
   // use functions in the loaded macro to create cuts
   AliRsnCutTrackQuality *cutQualityITS     = TrackQualityITS();
   AliRsnCutTrackQuality *cutQualityTPC     = TrackQualityTPC();
   AliRsnCutPIDITS       *cutPIDITSkaon     = PIDITS(isMC, AliRsnDaughter::kKaon, 0.0, 1E20, 3.0, "cutPIDITSkaon");
   AliRsnCutPIDTPC       *cutPIDTPCkaonLow  = PIDTPC(isMC, AliRsnDaughter::kKaon, 0.0, 0.350, 5.0, "cutPIDTPCkaonLow");
   AliRsnCutPIDTPC       *cutPIDTPCkaonHigh = PIDTPC(isMC, AliRsnDaughter::kKaon, 0.350, 1E20, 3.0, "cutPIDTPCkaonHigh");
   AliRsnCutPIDTOF       *cutPIDTOFkaon     = PIDTOF(AliRsnDaughter::kKaon, 3.0, !useTPC, "cutPIDTOFkaon");
   
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
   
   // choose cuts to add depending on used tracks
   for (Int_t i = 0; i < 4; i++) {
      // assign cuts to common manager for daughters
      AliRsnCutSet *cutSet = pair[i]->GetCutManager()->GetCommonDaughterCuts();
      // add cuts and define schemes depending on options
      if (addITS) {
         // if adding ITS standalone,
         // by default we use the quality check
         // and if required we add the PID (option "ITSPID")
         cutSet->AddCut(cutQualityITS);
         sprintf(schemeITS, "%s", qualityITS);
         if (useITS) {
            sprintf(schemeITS, "%s & %s", qualityITS, pidITS);
            cutSet->AddCut(cutPIDITSkaon);
         }
      }
      if (addTPC) {
         // if adding TPC+ITS tracks,
         // by default we use the quality check
         // and then we can use no PID, or both, 
         // or only one among TPC and TOF
         // the TPC PID cut is always an 'OR' of 
         // the two ones, in order to cover the full momentum range
         cutSet->AddCut(cutQualityTPC);
         sprintf(schemeTPC, "%s", qualityTPC);
         if (useTPC && useTOF) {
            cutSet->AddCut(cutPIDTPCkaonLow);
            cutSet->AddCut(cutPIDTPCkaonHigh);
            cutSet->AddCut(cutPIDTOFkaon);
            sprintf(schemeTPC, "%s & %s & %s", qualityTPC, pidTPC, pidTOF);
         } else if (useTPC) {
            cutSet->AddCut(cutPIDTPCkaonLow);
            cutSet->AddCut(cutPIDTPCkaonHigh);
            sprintf(schemeTPC, "%s & %s", qualityTPC, pidTPC);
         } else if (useTOF) {
            cutSet->AddCut(cutPIDTOFkaon);
            sprintf(schemeTPC, "%s & %s", qualityTPC, pidTOF);
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
      cutSet->SetCutScheme(scheme);
      ::Info("RsnConfigPhi", "Scheme for daughter cuts: %s", cutSet->GetCutScheme().Data());
   }
      
   // ==================================================================================================================
   // == PAIR CUTS =====================================================================================================
   // ==================================================================================================================
   
   // we define here our rapidity window
   AliRsnCutValue *cutRapidity = RapidityRange(pairDef[0], -0.5, 0.5);

   // in this case, we add the cut to the specific cut sets of all pairs
   // and we must then loop over all pairs, to add cuts to the related sets
   for (Int_t ipair = 0; ipair < 4; ipair++) {
      pair[ipair]->GetCutManager()->GetMotherCuts()->AddCut(cutRapidity);
      pair[ipair]->GetCutManager()->GetMotherCuts()->SetCutScheme(cutRapidity->GetName());
      ::Info("RsnConfigPhi", "Scheme for pair cuts: %s", pair[ipair]->GetCutManager()->GetMotherCuts()->GetCutScheme().Data());
   }
   
   // ==================================================================================================================
   // == FUNCTIONS FOR HISTOGRAMS ======================================================================================
   // ==================================================================================================================
   
   // we define a histogram with inv mass and other support binning values
   // and, when possible, a resolution on inv. mass
   
   AliRsnValue *axisIM      = AxisIM();
   AliRsnValue *axisPt      = AxisPt();
   AliRsnValue *axisRes     = AxisRes();
   AliRsnValue *axisMultSPD = AxisMultSPD();
   AliRsnValue *axisMultMC  = AxisMultMC();

   // create function for inv. mass and add axes
   AliRsnFunction *fcnIM = new AliRsnFunction;
   if (!fcnIM->AddAxis(axisIM)     ) return kFALSE;
   if (!fcnIM->AddAxis(axisPt)     ) return kFALSE;
   if (!fcnIM->AddAxis(axisMultSPD)) return kFALSE;
   fcnIM->UseSparse();
   
   // create function for inv. mass and add axes
   AliRsnFunction *fcnIMRes = new AliRsnFunction;
   if (!fcnIMRes->AddAxis(axisIM)     ) return kFALSE;
   if (!fcnIMRes->AddAxis(axisPt)     ) return kFALSE;
   if (!fcnIMRes->AddAxis(axisMultSPD)) return kFALSE;
   if (!fcnIMRes->AddAxis(axisRes)    ) return kFALSE;
   
   // add functions to pairs
   pair[0]->AddFunction(fcnIMRes);
   for (Int_t ipair = 1; ipair < 4; ipair++) pair[ipair]->AddFunction(fcnIM);
   
   // ==================================================================================================================
   // == CONCLUSION ====================================================================================================
   // ==================================================================================================================
   
   // add all created AliRsnPair objects to the AliRsnAnalysisManager in the task
   task->GetAnalysisManager()->Add(pair[1]);
   task->GetAnalysisManager()->Add(pair[2]);
   task->GetAnalysisManager()->Add(pair[3]);
   if (isMC) task->GetAnalysisManager()->Add(pair[0]);
   
   return kTRUE;
}
