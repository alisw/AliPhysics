AliRsnCutSet *AddRsnCommonEventCuts(AliAnalysisTaskSE *task=0) {

//   return 0;

   Bool_t valid;
   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Double_t eventCutVertex = AliRsnTrainManager::GetGlobalDbl("RsnEventCutPrimaryVertex",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);

   // primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   // we switch on the check for pileup
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", eventCutVertex, 0, kFALSE);
   if (collisionType==0) cutVertex->SetCheckPileUp(kTRUE);

   // primary vertex is always used
   AliRsnCutSet *commonEventCuts = new AliRsnCutSet("commonEventCuts", AliRsnTarget::kEvent);
   commonEventCuts->AddCut(cutVertex);
   commonEventCuts->SetCutScheme(cutVertex->GetName());

   // if task is mini
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsn = (AliRsnMiniAnalysisTask *)task;
      taskRsn->SetEventCuts(commonEventCuts);
   }
   return commonEventCuts;
}
