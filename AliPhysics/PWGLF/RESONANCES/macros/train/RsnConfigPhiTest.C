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
   AliRsnAnalysisTask *task,
   Bool_t              isMC,
   Bool_t              isMix,
   const char         *options
)
{
   void myError  (const char *msg) {::Error  ("RsnConfigPhi", msg);}
   void myWarning(const char *msg) {::Warning("RsnConfigPhi", msg);}
   void myInfo   (const char *msg) {::Info   ("RsnConfigPhi", msg);}

   if (!task) myError("NULL task");
   
   // ==================================================================================================================
   // == OPTIONS =======================================================================================================
   // ==================================================================================================================
   
   // Instead or getting confused with plenty of arguments in the macro (with default values),
   // we use a unique string of options with a set of conventional strings to set up the job:
   // -- "MC"/"DATA" --> what kind of sample
   // -- "ITS"/"TPC" --> what tracks to use (ITS standalone and/or TPC+ITS)
   // -- "xxxPID"    --> add the PID cuts for the detector xxx.
   // -- "MULT"      --> add axes for multiplicity
   //
   // In this point, these options are converted into boolean variables.
   
   TString opt(options);
   opt.ToUpper();
   opt.ReplaceAll(" ", "");
   
   Bool_t addITS  = opt.Contains("ITS");
   Bool_t addTPC  = opt.Contains("TPC");
   Bool_t useITS  = opt.Contains("ITSPID");
   Bool_t useTPC  = opt.Contains("TPCPID");
   Bool_t useTOF  = opt.Contains("TOFPID");
   Bool_t useMult = opt.Contains("MULT");
   
   // correct options when needed
   if (!addITS) useITS = kFALSE;
   if (!addTPC) useTPC = useTOF = kFALSE;
   
   // ==================================================================================================================
   // == DEFINITIONS ===================================================================================================
   // ==================================================================================================================
   
   // daughter definitions
   AliRsnDaughterDef *defKaonP = new AliRsnDaughterDef(AliRsnDaughter::kKaon, '+');
   AliRsnDaughterDef *defKaonM = new AliRsnDaughterDef(AliRsnDaughter::kKaon, '-');
   
   // pair definitions --> decay channels:
   // in our case, unlike-charged KK pairs for the signal, and like-charged ones for background
   AliRsnPairDef *pairDef[3];
   pairDef[0] = new AliRsnPairDef(defKaonP, defKaonM, 333, 1.019455); // unlike
   pairDef[1] = new AliRsnPairDef(defKaonP, defKaonP, 333, 1.019455); // like ++
   pairDef[2] = new AliRsnPairDef(defKaonM, defKaonM, 333, 1.019455); // like --

   // computation objects:
   // these are the objects which drive the computations, whatever it is (histos or tree filling)
   // and all tracks passed to them will be given the mass according to the reference pair definition
   // we create two unlike-sign pair computators, one for all pairs and another for true pairs (useful in MC)
   AliRsnLoopPair *pair[4];
   pair[0] = new AliRsnLoopPair(Form("%s_kaonP_kaonM_phi", opt.Data()), 0, 0, pairDef[0]); // unlike - true
   pair[1] = new AliRsnLoopPair(Form("%s_kaonP_kaonM_all", opt.Data()), 0, 0, pairDef[0]); // unlike - all
   pair[2] = new AliRsnLoopPair(Form("%s_kaonP_kaonP_all", opt.Data()), 0, 0, pairDef[1]); // like ++
   pair[3] = new AliRsnLoopPair(Form("%s_kaonM_kaonM_all", opt.Data()), 0, 0, pairDef[2]); // like --
   
   // loop on events
   AliRsnLoopEvent *event = new AliRsnLoopEvent(Form("%s_evt", opt.Data()));

   // set additional option for true pairs (slot [0])
   pair[0]->SetOnlyTrue(kTRUE);
   pair[0]->SetCheckDecay(kTRUE);
   
   // assign the ID of the entry lists to be used by each pair to get selected daughters
   // in our case, the AliRsnInputHandler contains only one list for selecting kaons
   for (Int_t i = 0; i < 4; i++) pair[i]->SetListID(0, 0);
   
   // ==================================================================================================================
   // == COMPUTED VALUES ===============================================================================================
   // ==================================================================================================================
   
   // All values which should be computed are defined here and passed to the computation objects,
   // since they define all that is computed bye each one, and, in case one output is a histogram
   // they define the binning and range for that value
   //
   // NOTE:
   // --> multiplicity bins have variable size
   
   Double_t mult[] = { 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13.,  14.,  15.,  16.,  17.,  18.,  19., 
                      20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 500.};
   Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
   
   AliRsnValue *axisIM      = new AliRsnValue("IM"  , AliRsnValue::kPairInvMass   ,  0.9, 1.4, 0.001);
   AliRsnValue *axisRes     = new AliRsnValue("RES" , AliRsnValue::kPairInvMassRes, -0.5, 0.5, 0.001);
   AliRsnValue *axisPt      = new AliRsnValue("PT"  , AliRsnValue::kPairPt        ,  0.0, 5.0, 0.1  );
   AliRsnValue *axisY       = new AliRsnValue("Y"   , AliRsnValue::kPairY         , -1.1, 1.1, 0.1  );
   AliRsnValue *axisMultESD = new AliRsnValue("MESD", AliRsnValue::kEventMultESDCuts, nmult, mult);
   AliRsnValue *axisMultSPD = new AliRsnValue("MSPD", AliRsnValue::kEventMultSPD    , nmult, mult);
   AliRsnValue *axisMultMC  = new AliRsnValue("MMC" , AliRsnValue::kEventMultMC     , nmult, mult);

   // add values to pairs
   // NOTE: in previous package version, they were added to functions
   for (Int_t i = 0; i < 4; i++) {
      if (i == 0) pair[i]->AddValue(axisRes);
      pair[i]->AddValue(axisIM);
      pair[i]->AddValue(axisPt);
      pair[i]->AddValue(axisY);
      if (useMult) {
         ::Info("RsnConfigPhi", "Adding multiplicity computations");
         pair[i]->AddValue(axisMultESD);
         pair[i]->AddValue(axisMultSPD);
         if (isMC) pair[i]->AddValue(axisMultMC);
         event->AddValue(axisMultESD);
         event->AddValue(axisMultSPD);
         if (isMC) event->AddValue(axisMultMC);
      }
   }
      
   // ==================================================================================================================
   // == PAIR CUTS =====================================================================================================
   // ==================================================================================================================
   
   // Rapidity cut
   // Only thing to consider is that it needs a support object to define mass
   AliRsnCutValue *cutRapidity = new AliRsnCutValue("cutY", AliRsnValue::kPairY, -0.5, 0.5);
   cutRapidity->GetValueObj()->SetSupportObject(pairDef[0]);

   // in this case, we add the cut to the specific cut sets of all pairs
   // and we must then loop over all pairs, to add cuts to the related sets
   for (Int_t ipair = 0; ipair < 4; ipair++) {
      pair[ipair]->GetPairCuts()->AddCut(cutRapidity);
      pair[ipair]->GetPairCuts()->SetCutScheme(cutRapidity->GetName());
      ::Info("RsnConfigPhi", "Scheme for pair cuts: %s", pair[ipair]->GetPairCuts()->GetCutScheme().Data());
   }
   
   // ==================================================================================================================
   // == OUTPUTS =======================================================================================================
   // ==================================================================================================================
   
   // now we define the outputs
   // in the new version of package, we use same syntax as TNtuple:
   // for each output object we define a name, and a list of variables
   // which must contain a combination of names of the axes added to the pair
   // and third argument is an enum to decide what kind of output we want

   AliRsnOutput *ntpMult = new AliRsnOutput("ntp1", "MSPD:MESD", AliRsnOutput::kNtuple);
   AliRsnOutput *ntpIMPt = new AliRsnOutput("ntp2", "IM:PT"    , AliRsnOutput::kNtuple);
   AliRsnOutput *fcnIMPt = new AliRsnOutput("out1", "IM:PT"    , AliRsnOutput::kHistoDefault);
   AliRsnOutput *fcnIM   = new AliRsnOutput("out2", "IM"       , AliRsnOutput::kHistoDefault);
   AliRsnOutput *fcnPt   = new AliRsnOutput("out3", "PT"       , AliRsnOutput::kHistoDefault);
   AliRsnOutput *fcnRes  = new AliRsnOutput("out4", "IM:PT:RES", AliRsnOutput::kHistoSparse);
   
   // add outputs to pairs
   for (Int_t ipair = 0; ipair < 4; ipair++) {
      if (!ipair) pair[ipair]->AddOutput(fcnRes);
      pair[ipair]->AddOutput(ntpIMPt);
      //pair[ipair]->AddOutput(fcnIMPt);
      //pair[ipair]->AddOutput(fcnIM);
      //pair[ipair]->AddOutput(fcnPt);
   }
   if (useMult) event->AddOutput(ntpMult);
   
   // ==================================================================================================================
   // == CONCLUSION ====================================================================================================
   // ==================================================================================================================
   
   // add all created objects to the task
   for (Int_t i = 0; i < 4; i++) {
      if (i == 0 && !isMC) continue;
      if (isMix && i != 1) continue;
      pair[i]->SetMixed(isMix);
      task->Add(pair[i]);
   }
   if (useMult) task->Add(event);
   
   return kTRUE;
}
