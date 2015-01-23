void RsnTrainCommonSettings(TString type,TString rsnPart,TString extraMacro="",TString extraMacroArgs="") {

   Bool_t valid;
   AliRsnTrainManager::GetGlobalStr("LegoTrainPath",valid);
   if (!valid) {
      TString legoTrainPath = "$ALICE_PHYSICS/PWGLF/RESONANCES/macros/lego_train";
      AliRsnTrainManager::SetGlobalStr("LegoTrainPath",legoTrainPath.Data());
   }

   // removing Option part fo Rsn particle
   if (rsnPart.Contains(":")) rsnPart.Remove(rsnPart.Index(":"),rsnPart.Length());
   AliRsnTrainManager::SetGlobalStr("RsnParticle",rsnPart.Data());

   // CollisionType (pp=0,PbPb=1,pPb=2)
   if (type.Contains("pp")) AliRsnTrainManager::SetGlobalInt("IsCollisionType",0);
   else if (type.Contains("PbPb")) AliRsnTrainManager::SetGlobalInt("IsCollisionType",1);
   else if (type.Contains("pPb")) AliRsnTrainManager::SetGlobalInt("IsCollisionType",2);

   // data type
   if (type.Contains("ESD")) AliRsnTrainManager::SetGlobalInt("IsESD",1);
   else AliRsnTrainManager::SetGlobalInt("IsESD",0);

   // flag if we are using MC
   if (type.Contains("MC")) AliRsnTrainManager::SetGlobalInt("IsMC",1);
   else AliRsnTrainManager::SetGlobalInt("IsMC",0);

   // flag if we want to use event Mixing
   if (type.Contains("MIX")) AliRsnTrainManager::SetGlobalInt("IsMixing",1);
   else AliRsnTrainManager::SetGlobalInt("IsMixing",0);

   // Use Rsn Mini
   if (type.Contains("MINI")) AliRsnTrainManager::SetGlobalInt("IsRsnMini",1);
   else AliRsnTrainManager::SetGlobalInt("IsRsnMini",0);


   // current RSN base defaults (Will be changed in future)
   if (!extraMacro.IsNull()) {
      extraMacro.ReplaceAll(".C","");
      Printf("Running Extra Macro %s(%s)",extraMacro.Data(),extraMacroArgs.Data());
      gROOT->ProcessLine(TString::Format("%s(%s)",extraMacro.Data(),extraMacroArgs.Data()).Data());
   }
   AliRsnTrainManager::SetGlobalInt("RsnQA",0,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnNumMix",5,kFALSE);
   AliRsnTrainManager::SetGlobalDbl("RsnEventCutPrimaryVertex",10.0,kFALSE);
   AliRsnTrainManager::SetGlobalStr("RsnLegoTrainCommonCutOption","mon",kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnPhysSelFilterBit",-1,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnCommonQualityCut",-1,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnUseRapidity",0,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnOutputFull",1,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnUseMCMomentum",0,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnUseMCMonitoring",0,kFALSE);
   AliRsnTrainManager::SetGlobalInt("RsnUseAOD049Patch",0,kFALSE);

   AliRsnTrainManager::SetGlobalDbl("RsnMixDiffMult",10.0,kFALSE);
   AliRsnTrainManager::SetGlobalDbl("RsnMixDiffVz",1.0,kFALSE);
   AliRsnTrainManager::SetGlobalDbl("RsnMixDiffAngle",-1.0,kFALSE);

   // expert options (don't change)
   AliRsnTrainManager::SetGlobalInt("RsnMixPrintRefresh",-1,kFALSE);

}
