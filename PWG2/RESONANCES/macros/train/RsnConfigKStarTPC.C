//
// *** Configuration script for KStar->Kpi analysis with 2010 runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t RsnConfigKStarTPC
(
   AliRsnAnalysisTask *task,
   Bool_t              isMC,
   Bool_t              useCentrality,
   AliRsnCutSet       *eventCuts
)
{
   if (!task) ::Error("RsnConfigKStarTPC", "NULL task");
   
   // we define here a suffix to differentiate names of different setups for the same resonance
   // and we define also the name of the list of tracks we want to select for the analysis
   // (if will fail if no lists with this name were added to the RsnInputHandler)
   const char *suffix     = "tpc";
   const char *listName1  = "kaonTPC";
   const char *listName2  = "pionTPC";
   Bool_t      useCharged = kTRUE;
   Int_t       listID1    = -1;
   Int_t       listID2    = -1;
   
   // find the index of the corresponding list in the RsnInputHandler
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *multi = dynamic_cast<AliMultiInputEventHandler*>(mgr->GetInputEventHandler());
   if (multi) {
      TObjArray *array = multi->InputEventHandlers();
      AliRsnInputHandler *rsn = (AliRsnInputHandler*)array->FindObject("rsnInputHandler");
      if (rsn) {
         AliRsnDaughterSelector *sel = rsn->GetSelector();
         listID1 = sel->GetID(listName1, useCharged);
         listID2 = sel->GetID(listName2, useCharged);
      }
   }
   if (listID1 >= 0 && listID2 >= 0)
      ::Info("RsnConfigKStarTPC.C", "Required lists '%s' and '%s' stay in positions %d and %d", listName1, listName2, listID1, listID2);
   else {
      ::Error("RsnConfigKStarTPC.C", "Required lists '%s' and '%s' absent in handler!", listName1, listName2);
      return kFALSE;
   }
            
   // ----------------------------------------------------------------------------------------------
   // -- DEFINITIONS -------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // PAIR DEFINITIONS:
   // this contains the definition of particle species and charge for both daughters of a resonance,
   // which are used for the following purposes:
   // --> species is used to assign the mass to the daughter (e.g. for building invariant mass)
   // --> charge is used to select what tracks to use when doing the computation loops
   // When a user wants to compute a like-sign background, he must define also a pair definition
   // for each like-sign: in case of charged track decays, we need one for ++ and one for --
   // Last two arguments are necessary only in some cases (but it is not bad to well initialize them):
   // --> PDG code of resonance, which is used for selecting true pairs, when needed
   // --> nominal resonance mass, which is used for computing quantities like Y or Mt
   AliRsnPairDef *kstar_kaonP_pionM = new AliRsnPairDef(AliRsnDaughter::kKaon, '+', AliRsnDaughter::kPion, '-', 313, 0.896);
   AliRsnPairDef *kstar_kaonM_pionP = new AliRsnPairDef(AliRsnDaughter::kKaon, '-', AliRsnDaughter::kPion, '+', 313, 0.896);
   AliRsnPairDef *kstar_kaonP_pionP = new AliRsnPairDef(AliRsnDaughter::kKaon, '+', AliRsnDaughter::kPion, '+', 313, 0.896);
   AliRsnPairDef *kstar_kaonM_pionM = new AliRsnPairDef(AliRsnDaughter::kKaon, '-', AliRsnDaughter::kPion, '-', 313, 0.896);

   // PAIR LOOPS:
   // these are the objects which drive the computations and fill the output histograms
   // each one requires to be initialized with an AliRsnPairDef object, which provided masses,
   // last argument tells if the pair is for mixing or not (this can be also set afterwards, anyway)
   const Int_t     nPairs = 8;
   Bool_t          addPair[nPairs] = {1, 1, 1, 1, 1, 1, 1, 1};
   AliRsnLoopPair *kstarLoop[nPairs];
   kstarLoop[0] = new AliRsnLoopPair(Form("%s_kstar_kaonP_pionM"     , suffix), kstar_kaonP_pionM, kFALSE);
   kstarLoop[1] = new AliRsnLoopPair(Form("%s_kstar_kaonP_pionM_true", suffix), kstar_kaonP_pionM, kFALSE);
   kstarLoop[2] = new AliRsnLoopPair(Form("%s_kstar_kaonP_pionM_mix" , suffix), kstar_kaonP_pionM, kTRUE );
   kstarLoop[3] = new AliRsnLoopPair(Form("%s_kstar_kaonM_pionP"     , suffix), kstar_kaonM_pionP, kFALSE);
   kstarLoop[4] = new AliRsnLoopPair(Form("%s_kstar_kaonM_pionP_true", suffix), kstar_kaonM_pionP, kFALSE);
   kstarLoop[5] = new AliRsnLoopPair(Form("%s_kstar_kaonM_pionP_mix" , suffix), kstar_kaonM_pionP, kTRUE );
   kstarLoop[6] = new AliRsnLoopPair(Form("%s_kstar_kaonP_pionP"     , suffix), kstar_kaonP_pionP, kFALSE);
   kstarLoop[7] = new AliRsnLoopPair(Form("%s_kstar_kaonM_pionM"     , suffix), kstar_kaonM_pionM, kFALSE);

   // set additional option for true pairs
   // 1) we select only pairs coming from the same mother, which must have the right PDG code (from pairDef)
   // 2) we select only pairs decaying according to the right channel (from pairDef species+charge definitions)
   kstarLoop[1]->SetOnlyTrue(kTRUE);
   kstarLoop[1]->SetCheckDecay(kTRUE);
   kstarLoop[4]->SetOnlyTrue(kTRUE);
   kstarLoop[4]->SetCheckDecay(kTRUE);
   
   // don't add true pairs if not MC
   if (!isMC) addPair[1] = addPair[4] = 0;
   
   // ----------------------------------------------------------------------------------------------
   // -- PAIR CUTS ---------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // for pairs we define a rapidity windows, defined through a cut
   // --> NOTE: it needs a support AliRsnPairDef from which it takes the mass
   AliRsnValueStd *valRapidity = new AliRsnValueStd("valY", AliRsnValueStd::kPairY);
   AliRsnCutValue *cutRapidity = new AliRsnCutValue("kstar_cutY", -0.5, 0.5, isMC);
   valRapidity->SetSupportObject(kstar_kaonP_pionM);
   cutRapidity->SetValueObj(valRapidity);
   
   // add the cut to a cut set (will be simple, there's only one)
   AliRsnCutSet *pairCuts = new AliRsnCutSet("kstar_pairCuts", AliRsnTarget::kMother);
   pairCuts->AddCut(cutRapidity);
   pairCuts->SetCutScheme(cutRapidity->GetName());
   
   // ----------------------------------------------------------------------------------------------
   // -- COMPUTED VALUES & OUTPUTS -----------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // All values which should be computed are defined here and passed to the computation objects,
   // since they define all that is computed bye each one, and, in case one output is a histogram
   // they define the binning and range for that value
   //
   // NOTE:
   // --> multiplicity bins have variable size
   
   Double_t mult[] = { 0.,  1.,  2.,  3.,  4.,  5.,   6.,   7.,   8.,   9.,  10., 11., 12., 13., 
                      14., 15., 16., 17., 18., 19.,  20.,  21.,  22.,  23.,  24.,  25., 30., 35., 
                      40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 500.};
   Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
   
   AliRsnValueStd *axisIM      = new AliRsnValueStd("kIM" , AliRsnValueStd::kPairInvMass       ,   0.5,   1.5, 0.0025);
   AliRsnValueStd *axisRes     = new AliRsnValueStd("RES" , AliRsnValueStd::kPairInvMassRes    ,  -0.5,   0.5, 0.001);
   AliRsnValueStd *axisPt      = new AliRsnValueStd("PT"  , AliRsnValueStd::kPairPt            ,   0.0,   5.0, 0.1  );
   AliRsnValueStd *axisCentV0  = new AliRsnValueStd("CNT" , AliRsnValueStd::kEventCentralityV0 ,   0.0, 100.0, 5.0  );
   AliRsnValueStd *axisMultESD = new AliRsnValueStd("MESD", AliRsnValueStd::kEventMultESDCuts  , nmult, mult);
   AliRsnValueStd *axisMultSPD = new AliRsnValueStd("MSPD", AliRsnValueStd::kEventMultSPD      , nmult, mult);
   AliRsnValueStd *axisMultTRK = new AliRsnValueStd("MTRK", AliRsnValueStd::kEventMult         , nmult, mult);
   AliRsnValueStd *axisMultMC  = new AliRsnValueStd("MMC" , AliRsnValueStd::kEventMultMC       , nmult, mult);

   // create outputs:
   // we define one for true pairs, where we add resolution, and another without it, for all others
   // it seems that it is much advantageous to use sparse histograms when adding more than 2 axes
   AliRsnListOutput *out[2];
   out[0] = new AliRsnListOutput("res"  , AliRsnListOutput::kHistoSparse);
   out[1] = new AliRsnListOutput("nores", AliRsnListOutput::kHistoSparse);
   
   // add values to outputs:
   // if centrality is required, we add it only, otherwise we add all multiplicities
   // other axes (invmass, pt) are always added
   for (Int_t i = 0; i < 2; i++) {
      out[i]->AddValue(axisIM);
      out[i]->AddValue(axisPt);
      if (useCentrality) {
         ::Info("RsnConfigKStarTPC.C", "Adding centrality axis");
         out[i]->AddValue(axisCentV0);
      } else {
         ::Info("RsnConfigKStarTPC.C", "Adding multiplicity axes");
         //out[i]->AddValue(axisMultESD);
         //out[i]->AddValue(axisMultSPD);
         out[i]->AddValue(axisMultTRK);
         if (isMC) out[i]->AddValue(axisMultMC);
      }
   }
   // resolution only in the first
   out[0]->AddValue(axisRes);
   
   // ----------------------------------------------------------------------------------------------
   // -- ADD SETTINGS TO LOOPS AND LOOPS TO TASK ---------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   for (Int_t ip = 0; ip < nPairs; ip++) {
      // skip pairs not to be added
      if (!addPair[ip]) continue;
      // assign list IDs
      kstarLoop[ip]->SetListID(0, listID1);
      kstarLoop[ip]->SetListID(1, listID2);
      // assign event cuts
      kstarLoop[ip]->SetEventCuts(eventCuts);
      // assign pair cuts
      kstarLoop[ip]->SetPairCuts(pairCuts);
      // assign outputs
      if (ip != 1)
         kstarLoop[ip]->AddOutput(out[1]);
      else
         kstarLoop[ip]->AddOutput(out[0]);
      // add to task
      task->Add(kstarLoop[ip]);
   }
   
   return kTRUE;
}
