void MakeADQAParamEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB") {
//========================================================================
//
// Steering macro for AD reconstruction parameters
//
// Author: Michal Broz
//
//========================================================================

  const char* macroname = "MakeADQAParamEntry.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(outputCDB);
  cdb->SetRun(0);

  AliADQAParam * ADQAParam = new AliADQAParam();
  ADQAParam->SetMaxPedDiff(5);
  ADQAParam->SetMaxPedWidth(1.5);
  ADQAParam->SetSatMed(0.7);
  ADQAParam->SetSatHigh(0.9);
  ADQAParam->SetSatHuge(0.99);
  
  ADQAParam->SetAsynchronBB(5.0);
  ADQAParam->SetAsynchronBG(10.0);
  
  ADQAParam->SetMaxBBVariation(0.2);
  ADQAParam->SetMaxBGVariation(0.2);
  
  ADQAParam->SetMaxNoFlagRate(0.2);
  
  ADQAParam->SetChargeChannelZoomMin(0);
  ADQAParam->SetChargeChannelZoomMax(50);
  
  ADQAParam->SetNTdcTimeBins(3062); 
  ADQAParam->SetTdcTimeMin(0.976562); 
  ADQAParam->SetTdcTimeMax(300.0); 
  
  ADQAParam->SetNChargeChannelBins(10000);
  ADQAParam->SetChargeChannelMin(0);
  ADQAParam->SetChargeChannelMax(10000);
  
  ADQAParam->SetNChargeSideBins(7000);
  ADQAParam->SetChargeSideMin(10);
  ADQAParam->SetChargeSideMax(70000);

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Michal Broz");
  md->SetComment("Reconstruction parameters for AD");
  //md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  //md->SetBeamPeriod(0);
  AliCDBId id("AD/Calib/QAParam", 0, AliCDBRunRange::Infinity());
  cdb->Put(ADQAParam, id, md);

  return;
}
