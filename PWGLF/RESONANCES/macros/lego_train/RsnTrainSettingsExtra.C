void RsnTrainSettingsExtra(Double_t primaryVertex,
                           Int_t useCommonQualityCut,
                           Int_t numMix,
                           Int_t filterBit,
                           Int_t useRapidity,
                           Int_t useAOD49Patch,
                           Double_t useMixDiffMult,
                           Double_t useMixDiffVz,
                           Double_t useMixDiffAngle,
                           Int_t printRefresh=-1,
                           Int_t useMixLike=0) {

   AliRsnTrainManager::SetGlobalDbl("RsnEventCutPrimaryVertex",primaryVertex);
   AliRsnTrainManager::SetGlobalInt("RsnCommonQualityCut",useCommonQualityCut);
   AliRsnTrainManager::SetGlobalInt("RsnPhysSelFilterBit",filterBit);
   AliRsnTrainManager::SetGlobalInt("RsnUseRapidity",useRapidity);
   AliRsnTrainManager::SetGlobalInt("RsnUseAOD049Patch",useAOD49Patch);

   // for now we will use only on/off (0/1), maybe we can use number also, but let's see if neede
   AliRsnTrainManager::SetGlobalInt("RsnNumMix",numMix);
   AliRsnTrainManager::SetGlobalDbl("RsnMixDiffMult",useMixDiffMult);
   AliRsnTrainManager::SetGlobalDbl("RsnMixDiffVz",useMixDiffVz);
   AliRsnTrainManager::SetGlobalDbl("RsnMixDiffAngle",useMixDiffAngle);

   AliRsnTrainManager::SetGlobalInt("RsnMixPrintRefresh",printRefresh);
   AliRsnTrainManager::SetGlobalInt("RsnMixLike",useMixLike);

   return;
}
