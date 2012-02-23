AliRsnCutSet *AddRsnCommonPairCuts(AliAnalysisTaskSE *task=0,Bool_t isPP=kTRUE) {
   AliRsnCutSet *cutsPair = 0;
   //   TODO
//   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
//   //    cutY->SetRangeD(-0.5, 0.5);
//   cutY->SetRangeD(-0.9, 0.9);
//
//   cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
//   cutsPair->AddCut(cutY);
//   cutsPair->SetCutScheme(cutY->GetName());

   return cutsPair;
}
