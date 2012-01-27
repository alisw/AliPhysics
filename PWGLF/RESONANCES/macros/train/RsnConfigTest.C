Bool_t RsnConfigTest
(
   AliRsnAnalysisTask *task,
   Bool_t              isMC
)
{
   // find the index of the corresponding list in the RsnInputHandler
   const char *listNameQuality = "qualityTPC";
   const char *listNamePID     = "kaonTPC";
   Int_t       qualityID       = -1;
   Int_t       pidID           = -1;
   AliAnalysisManager        *mgr   = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *multi = dynamic_cast<AliMultiInputEventHandler*>(mgr->GetInputEventHandler());
   if (multi) {
      TObjArray *array = multi->InputEventHandlers();
      TObjArrayIter next(array);
      TObject *obj;
      while ( (obj = next()) ) {
         if (obj->InheritsFrom(AliRsnInputHandler::Class())) {
            AliRsnInputHandler *rsn = (AliRsnInputHandler*)obj;
            AliRsnDaughterSelector *sel = rsn->GetSelector();
            qualityID = sel->GetID(listNameQuality, kTRUE);
            pidID = sel->GetID(listNamePID, kTRUE);
         }
      }
   }
   if (qualityID < 0) {
      ::Error("RsnConfigTest", "Selector does not contain list for quality only");
      return kFALSE;
   }
   if (pidID < 0) {
      ::Error("RsnConfigTest", "Selector does not contain list for quality+PID");
      return kFALSE;
   }
   ::Info("RsnConfigTest", "ID for cut set named '%10s' = %d", listNameQuality, qualityID);
   ::Info("RsnConfigTest", "ID for cut set named '%10s' = %d", listNamePID, pidID);
   
   // add pair computation
   //AddPairLoop(task, isMC, pidID, pidID, "test");
   
   // add monitor computation
   AddMonitorLoop(task, isMC, qualityID, pidID, "test");
   return kTRUE;
}

AliRsnCutSet* EventCuts() 
{
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
   
   return eventCuts;
}

AliRsnCutSet* PairCuts(AliRsnPairDef *support)
{
   // for pairs we define a rapidity windows, defined through a cut
   // --> NOTE: it needs a support AliRsnPairDef from which it takes the mass
   AliRsnValueStd *valRapidity = new AliRsnValueStd("valY", AliRsnValueStd::kPairY);
   AliRsnCutValue *cutRapidity = new AliRsnCutValue("cutY", -0.5, 0.5, kFALSE);
   valRapidity->SetSupportObject(support);
   cutRapidity->SetValueObj(valRapidity);
   
   // cut set
   AliRsnCutSet *pairCuts = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   pairCuts->AddCut(cutRapidity);
   pairCuts->SetCutScheme(cutRapidity->GetName());
}

void AddPairOutput(AliRsnLoopPair *pair)
{
   // axes
   AliRsnValueStd *axisIM = new AliRsnValueStd("IM", AliRsnValueStd::kPairInvMass, 0.9, 1.4, 0.001);
   AliRsnValueStd *axisPt = new AliRsnValueStd("PT", AliRsnValueStd::kPairPt     , 0.0, 5.0, 0.1  );
     
   // output: 2D histogram of inv. mass vs. pt
   AliRsnListOutput *outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoDefault);
   outPair->AddValue(axisIM);
   outPair->AddValue(axisPt);
   
   // add outputs to loop
   pair->AddOutput(outPair);
}

void AddMonitorOutput(AliRsnLoopDaughter *mon)
{
   // axes
   AliRsnValueStd *axisMomTPC = new AliRsnValueStd("pTPC", AliRsnValueStd::kTrackPtpc     , 0.0,   5.0, 0.01 );
   AliRsnValueStd *axisSigTPC = new AliRsnValueStd("sTPC", AliRsnValueStd::kTrackTPCsignal, 0.0, 500.0, 2.0  );
   
   // output: 2D histogram of TPC signal vs. TPC momentum
   AliRsnListOutput *outMonitor = new AliRsnListOutput("mon", AliRsnListOutput::kHistoDefault);
   outMonitor->AddValue(axisMomTPC);
   outMonitor->AddValue(axisSigTPC);
   
   // add outputs to loop
   mon->AddOutput(outMonitor);
}

void AddPairLoop(AliRsnAnalysisTask *task, Bool_t isMC, Int_t listID1, Int_t listID2, const char *suffix = "")
{
   // pair definition
   AliRsnPairDef  *pairDef    = new AliRsnPairDef(AliRsnDaughter::kKaon, '+', AliRsnDaughter::kKaon, '-', 333, 1.019455);
   
   // loop object creation
   AliRsnLoopPair *loopPhi    = new AliRsnLoopPair(Form("%s_unlike", suffix), pairDef, kFALSE);
   AliRsnLoopPair *loopPhiMix = new AliRsnLoopPair(Form("%s_unlike", suffix), pairDef, kTRUE );
   
   // assign ID of selector lists to be used (KK --> the same)
   loopPhi->SetListID(0, listID1);
   loopPhi->SetListID(1, listID2);
   loopPhiMix->SetListID(0, listID1);
   loopPhiMix->SetListID(1, listID2);
   
   // assign cuts for events
   AliRsnCutSet *eventCuts = EventCuts();
   loopPhi->SetEventCuts(eventCuts);
   loopPhiMix->SetEventCuts(eventCuts);
   
   // assign cuts for pairs
   AliRsnCutSet *pairCuts = PairCuts(pairDef);
   loopPhi->SetPairCuts(eventCuts);
   loopPhiMix->SetPairCuts(eventCuts);
   
   // add outputs
   AddPairOutput(loopPhi);
   AddPairOutput(loopPhiMix);
   
   // add loops to task
   task->Add(loopPhi);
   task->Add(loopPhiMix);
}

void AddMonitorLoop(AliRsnAnalysisTask *task, Bool_t isMC, Int_t listID1, Int_t listID2, const char *suffix = "")
{
   // monitor definition
   AliRsnDaughterDef *tracks = new AliRsnDaughterDef(AliRsnDaughter::kTrack /*'+' or '-'*/);
   
   // loop object
   AliRsnLoopDaughter *loopMon1 = new AliRsnLoopDaughter(Form("%s_mon1", suffix), listID1, tracks);
   AliRsnLoopDaughter *loopMon2 = new AliRsnLoopDaughter(Form("%s_mon2", suffix), listID2, tracks);
   
   // add cuts on events
   AliRsnCutSet *eventCuts = EventCuts();
   loopMon1->SetEventCuts(eventCuts);
   loopMon2->SetEventCuts(eventCuts);
   
   // add monitors
   AddMonitorOutput(loopMon1);
   AddMonitorOutput(loopMon2);
   
   // add loop to task
   task->Add(loopMon1);
   task->Add(loopMon2);
}
