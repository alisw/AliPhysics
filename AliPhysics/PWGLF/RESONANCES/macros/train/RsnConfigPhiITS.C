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
Bool_t RsnConfigPhiITS
(
   AliRsnAnalysisTask *task,
   Bool_t              isMC,
   Bool_t              isMix,
   Bool_t              useCentrality
)
{
   void myError  (const char *msg) {::Error  ("RsnConfigPhiITS", msg);}
   void myWarning(const char *msg) {::Warning("RsnConfigPhiITS", msg);}
   void myInfo   (const char *msg) {::Info   ("RsnConfigPhiITS", msg);}

   if (!task) myError("NULL task");
   
   const char *suffix = "itsstd";
      
   // ==================================================================================================================
   // == DEFINITIONS ===================================================================================================
   // ==================================================================================================================
   
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
   AliRsnLoopPair *pair[4];
   pair[0] = new AliRsnLoopPair(Form("%s_kaonP_kaonM_phi", suffix), 0, 0, pairDef[0]); // unlike - true
   pair[1] = new AliRsnLoopPair(Form("%s_kaonP_kaonM_all", suffix), 0, 0, pairDef[0]); // unlike - all
   pair[2] = new AliRsnLoopPair(Form("%s_kaonP_kaonP_all", suffix), 0, 0, pairDef[1]); // like ++
   pair[3] = new AliRsnLoopPair(Form("%s_kaonM_kaonM_all", suffix), 0, 0, pairDef[2]); // like --

   // set additional option for true pairs (slot [0])
   pair[0]->SetOnlyTrue(kTRUE);
   pair[0]->SetCheckDecay(kTRUE);
   
   // assign the ID of the entry lists to be used by each pair to get selected daughters
   // in our case, the AliRsnInputHandler contains only one list for selecting kaons
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *multi = dynamic_cast<AliMultiInputEventHandler*>(mgr->GetInputEventHandler());
   if (!multi) {
      myError("Needed a multi input handler!");
      return kFALSE;
   }
   TObjArray *array = multi->InputEventHandlers();
   AliRsnInputHandler *rsn = (AliRsnInputHandler*)array->FindObject("rsnInputHandler");
   if (!rsn) {
      myError("Needed an RSN event handler");
      return kFALSE;
   }
   AliRsnDaughterSelector *sel = rsn->GetSelector();
   Int_t id = sel->GetID("kaonTPC", kTRUE);
   if (id < 0) {
      myError("Kaons are not added in the selector");
      return kFALSE;
   }
   myInfo(Form("Selected list is in position #%d", id));
   for (Int_t i = 0; i < 4; i++) {
      pair[i]->SetListID(0, id);
      pair[i]->SetListID(1, id);
   }
   
   // ----------------------------------------------------------------------------------------------
   // -- EVENT CUTS --------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // in the function for events, we don't cut on centrality or multiplicity, 
   // since it becomes an axis of the output histogram

   // primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   // we switch on the check for pileup
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   cutVertex->SetCheckPileUp(kTRUE);
      
   // primary vertex is always used
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(cutVertex->GetName());
   
   // set cuts for each loop
   for (Int_t i = 0; i < 4; i++) {
      pair[i]->SetEventCuts(eventCuts);
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
   }
   
   // ==================================================================================================================
   // == COMPUTED VALUES & OUTPUTS =====================================================================================
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
   AliRsnValue *axisMultESD = new AliRsnValue("MESD", AliRsnValue::kEventMultESDCuts, nmult, mult);
   AliRsnValue *axisMultSPD = new AliRsnValue("MSPD", AliRsnValue::kEventMultSPD    , nmult, mult);
   AliRsnValue *axisMultMC  = new AliRsnValue("MMC" , AliRsnValue::kEventMultMC     , nmult, mult);
   AliRsnValue *axisCentV0  = new AliRsnValue("CNT" , AliRsnValue::kEventCentralityV0 , 0.0, 100.0, 10.0);

   // create outputs
   AliRsnListOutput *out[2];
   AliRsnListOutput *out[0] = new AliRsnListOutput("phi", AliRsnListOutput::kHistoSparse);
   AliRsnListOutput *out[1] = new AliRsnListOutput("all", AliRsnListOutput::kHistoSparse);
   
   // add values to outputs
   out[0]->AddValue(axisRes);
   for (Int_t i = 0; i < 2; i++) {
      out[i]->AddValue(axisIM);
      out[i]->AddValue(axisPt);
      if (useCentrality) {
         out[i]->AddValue(axisCentV0);
      } else {
         out[i]->AddValue(axisMultESD);
         out[i]->AddValue(axisMultSPD);
         if (isMC) out[i]->AddValue(axisMultMC);
      }
   }
   
   // add outputs to pairs
   pair[0]->AddOutput(out[0]);
   for (Int_t ipair = 1; ipair < 4; ipair++) {
      pair[ipair]->AddOutput(out[1]);
   }
   
   // ==================================================================================================================
   // == CONCLUSION ====================================================================================================
   // ==================================================================================================================
   
   if (isMC && !isMix) task->Add(pair[0]);
   task->Add(pair[1]);
   if (!isMix) {
      task->Add(pair[2]);
      task->Add(pair[3]);
   }
   
   return kTRUE;
}
