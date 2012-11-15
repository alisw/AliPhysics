void RsnTrainSettingsExtra(Double_t primaryVertex,
                           Int_t useCommonQualityCut,
                           Int_t numMix,
                           Int_t filterBit,
                           Int_t useRapidity,
                           Int_t useAOD49Patch,
                           Int_t useMixDiffMult,
                           Int_t useMixDiffVz,
                           Int_t useMixDiffAngle) {

   AliRsnTrainManager::SetGlobalDbl("RsnEventCutPrimaryVertex",primaryVertex);
   AliRsnTrainManager::SetGlobalInt("RsnCommonQualityCut",useCommonQualityCut);
   AliRsnTrainManager::SetGlobalInt("RsnPhysSelFilterBit",filterBit);
   AliRsnTrainManager::SetGlobalInt("RsnUseRapidity",useRapidity);
   AliRsnTrainManager::SetGlobalInt("RsnUseAOD049Patch",useAOD49Patch);

   // for now we will use only on/off (0/1), maybe we can use number also, but let's see if neede
   AliRsnTrainManager::SetGlobalInt("RsnNumMix",numMix);
   AliRsnTrainManager::SetGlobalInt("RsnMixDiffMult",useMixDiffMult);
   AliRsnTrainManager::SetGlobalInt("RsnMixDiffVz",useMixDiffVz);
   AliRsnTrainManager::SetGlobalInt("RsnMixDiffAngle",useMixDiffAngle);


   return;
}
