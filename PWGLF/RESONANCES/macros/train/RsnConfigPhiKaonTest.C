//
// Test config macro for RSN package.
// It configures:
// 1) a monitor for all tracks passing quality cuts
// 2) a monitor for all tracks passing quality + PID cuts
// 3) an unlike-sign invariant-mass + pt distribution for K+K- pairs
//
Bool_t RsnConfigPhiKaonTest
(
   AliRsnAnalysisTask *task,
   Bool_t              isMC
)
{
   if (!task) {
      ::Error("RsnConfigPhiKaonTest.C", "NULL task");
      return kFALSE;
   }
   
   const char *suffix = "test";
      
   // ----------------------------------------------------------------------------------------------
   // -- DEFINITIONS -------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // daughter definition for monitor loops
   // since it is intended to loop over all 'track' like objects (in the sense that we exclude V0s and cascades),
   // we initialize it using the constructor that requires an AliRsnDaughter::EType and a charge, but since
   // we want to loop over both charges, we set it to anything which is not '+' '-' or '0', which are tokens for
   // selecting only positive, only negative or only neutral
   AliRsnDaughterDef *tracks = new AliRsnDaughterDef(AliRsnDaughter::kTrack, 0);
   
   // definition of pair decay tree for phi resonance
   // here we *must* specify a particle species and a charge, in order to check the decay tree
   // last arguments are the PDG code and nominal mass of the resonance, which are needed when
   // one wants to select true pairs only and/or he wants to compute rapidity or Mt 
   AliRsnPairDef *pairDef = new AliRsnPairDef(AliRsnDaughter::kKaon, '+', AliRsnDaughter::kKaon, '-', 333, 1.019455);
   
   // definition of loop objects:
   // (a) 1 monitor for all tracks passing quality cuts
   // (b) 1 monitor for all tracks passing quality+PID cuts
   // (c) 1 pair filled with all tracks passing same cuts as (b)
   // (d) 1 pair like (c) but for mixing
   // (e) 1 pair like (c) but with true pairs only
   // NOTE: (c) and (d) are instantiated with same settings, they will be made
   //       different after some settings done in second moment
   AliRsnLoopDaughter *loopQuality = new AliRsnLoopDaughter(Form("%s_mon_quality", suffix), 0, tracks);
   AliRsnLoopDaughter *loopPID     = new AliRsnLoopDaughter(Form("%s_mon_pid"    , suffix), 0, tracks);
   AliRsnLoopPair     *loopPhi     = new AliRsnLoopPair    (Form("%s_unlike"     , suffix), pairDef);
   AliRsnLoopPair     *loopPhiMix  = new AliRsnLoopPair    (Form("%s_unlike"     , suffix), pairDef);
   AliRsnLoopPair     *loopPhiTrue = new AliRsnLoopPair    (Form("%s_trues"      , suffix), pairDef);
   
   // set additional option for true pairs (slot [0])
   loopPhiTrue->SetOnlyTrue(kTRUE);
   loopPhiTrue->SetCheckDecay(kTRUE);
   
   // set mixing options
   loopPhi    ->SetMixed(kFALSE);
   loopPhiMix ->SetMixed(kTRUE);
   loopPhiTrue->SetMixed(kFALSE);
   
   // assign the ID of the entry lists to be used by each pair to get selected daughters
   // in our case, the AliRsnInputHandler contains only one list for selecting kaons
   Int_t idQuality, idPID;
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
   idQuality = sel->GetID("qualityTPC", kTRUE);
   idPID = sel->GetID("kaonTPC", kTRUE);
   if (idQuality < 0 || idPID < 0) {
      myError("List problems");
      return kFALSE;
   }
   loopQuality->SetListID(idQuality);
   loopPID    ->SetListID(idPID);
   loopPhi    ->SetListID(0, idPID);
   loopPhi    ->SetListID(1, idPID);
   loopPhiMix ->SetListID(0, idPID);
   loopPhiMix ->SetListID(1, idPID);
   loopPhiTrue->SetListID(0, idPID);
   loopPhiTrue->SetListID(1, idPID);
   
   // ----------------------------------------------------------------------------------------------
   // -- EVENT CUTS --------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------

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
   
   // add the event cuts to all loops
   loopQuality->SetEventCuts(eventCuts);
   loopPID    ->SetEventCuts(eventCuts);
   loopPhi    ->SetEventCuts(eventCuts);
   loopPhi    ->SetEventCuts(eventCuts);
   loopPhiMix ->SetEventCuts(eventCuts);
   loopPhiMix ->SetEventCuts(eventCuts);
   loopPhiTrue->SetEventCuts(eventCuts);
   loopPhiTrue->SetEventCuts(eventCuts);
   
   // ----------------------------------------------------------------------------------------------
   // -- PAIR CUTS ---------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // for pairs we define a rapidity windows, defined through a cut
   // --> NOTE: it needs a support AliRsnPairDef from which it takes the mass
   AliRsnValueStd *valRapidity = new AliRsnValueStd("valY", AliRsnValueStd::kPairY);
   AliRsnCutValue *cutRapidity = new AliRsnCutValue("cutY", -0.5, 0.5, isMC);
   valRapidity->SetSupportObject(pairDef);
   cutRapidity->SetValueObj(valRapidity);
   
   // cut set
   AliRsnCutSet *pairCuts = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   pairCuts->AddCut(cutRapidity);
   pairCuts->SetCutScheme(cutRapidity->GetName());
   
   // add cut to pair loops only
   loopPhi    ->SetPairCuts(pairCuts);
   loopPhi    ->SetPairCuts(pairCuts);
   loopPhiMix ->SetPairCuts(pairCuts);
   loopPhiMix ->SetPairCuts(pairCuts);
   loopPhiTrue->SetPairCuts(pairCuts);
   loopPhiTrue->SetPairCuts(pairCuts);
   
   // ----------------------------------------------------------------------------------------------
   // -- COMPUTED VALUES & OUTPUTS -----------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   AliRsnValueStd *axisIM     = new AliRsnValueStd("IM"  , AliRsnValueStd::kPairInvMass   , 0.9,   1.4, 0.001);
   AliRsnValueStd *axisPt     = new AliRsnValueStd("PT"  , AliRsnValueStd::kPairPt        , 0.0,   5.0, 0.1  );
   AliRsnValueStd *axisMomTPC = new AliRsnValueStd("pTPC", AliRsnValueStd::kTrackPtpc     , 0.0,   5.0, 0.01 );
   AliRsnValueStd *axisSigTPC = new AliRsnValueStd("sTPC", AliRsnValueStd::kTrackTPCsignal, 0.0, 500.0, 2.0  );
   
   // output for monitors:
   // 2D histogram with TPC signal vs TPC momentum
   AliRsnListOutput *outMonitor = new AliRsnListOutput("mon", AliRsnListOutput::kHistoDefault);
   outMonitor->AddValue(axisMomTPC);
   outMonitor->AddValue(axisSigTPC);
   
   // output for pairs:
   // 2D histogram with inv.mass vs pt
   AliRsnListOutput *outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoDefault);
   outPair->AddValue(axisIM);
   outPair->AddValue(axisPt);
   
   // add outputs to loops
   loopQuality->AddOutput(outMonitor);
   loopPID    ->AddOutput(outMonitor);
   loopPhi    ->AddOutput(outPair);
   loopPhiMix ->AddOutput(outPair);
   loopPhiTrue->AddOutput(outPair);
   
   // ----------------------------------------------------------------------------------------------
   // -- CONCLUSION --------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   task->Add(loopQuality);
   task->Add(loopPID    );
   task->Add(loopPhi    );
   task->Add(loopPhiMix );
   task->Add(loopPhiTrue);
   
   return kTRUE;
}
