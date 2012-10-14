void RsnTrainCommonSettings(TString type) {

//   TString legoTrainPath = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";
//   AliRsnTrainManager::SetGlobalStr("LegoTrainPath",legoTrainPath.Data());

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
   AliRsnTrainManager::SetGlobalInt("IsMixing",0);
   // number of mixing
   AliRsnTrainManager::SetGlobalInt("NumMix",5);

   // Use Rsn Mini
   if (type.Contains("MINI")) AliRsnTrainManager::SetGlobalInt("IsRsnMini",1);
   else AliRsnTrainManager::SetGlobalInt("IsRsnMini",0);


   // current RSN base defaults (Will be changed in future)

   AliRsnTrainManager::SetGlobalDbl("RsnEventCutPrimaryVertex",10.0);
   AliRsnTrainManager::SetGlobalStr("RsnLegoTrainCommonCutOption","mon");
   AliRsnTrainManager::SetGlobalInt("RsnPhysSelFilterBit",-1);
   AliRsnTrainManager::SetGlobalInt("RsnCommonQualityCut",-1);
   AliRsnTrainManager::SetGlobalInt("RsnUseRapidity",0);
   AliRsnTrainManager::SetGlobalInt("RsnOutputFull",1);
   AliRsnTrainManager::SetGlobalInt("RsnUseMCMomentum",0);
   AliRsnTrainManager::SetGlobalInt("RsnUseMCMonitoring",0);

   // expert options (don't change)
   AliRsnTrainManager::SetGlobalInt("RsnMixPrintRefresh",-1);

}
