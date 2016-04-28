void MakeADQAParamEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB") {
//========================================================================
//
// Steering macro for AD QA parameters
//
// Author: Michal Broz
//
//========================================================================

  const char* macroname = "MakeADQAParamEntry.C";

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(outputCDB);
  cdb->SetRun(0);

  AliADQAParam * ADQAParam = new AliADQAParam();
  ADQAParam->SetMaxPedDiff(2);
  ADQAParam->SetMaxPedWidth(1.5);
  ADQAParam->SetSatMed(0.2);
  ADQAParam->SetSatHigh(0.5);
  ADQAParam->SetSatHuge(0.7);
  
  ADQAParam->SetAsynchronBB(5.0);
  ADQAParam->SetAsynchronBG(10.0);
  
  ADQAParam->SetMaxBBVariation(0.2);
  ADQAParam->SetMaxBGVariation(0.2);
  
  ADQAParam->SetMaxNoFlagRate(0.2);
  ADQAParam->SetMaxNoTimeRate(0.05);
  
  ADQAParam->SetChargeChannelZoomMin(0);
  ADQAParam->SetChargeChannelZoomMax(50);
  
  ADQAParam->SetNTdcTimeBins(3062); 
  ADQAParam->SetTdcTimeMin(0.976562); 
  ADQAParam->SetTdcTimeMax(300.0); 
  
  ADQAParam->SetNTdcWidthBins(159); 
  ADQAParam->SetTdcWidthMin(6.25); 
  ADQAParam->SetTdcWidthMax(1000.0); 
  
  ADQAParam->SetNChargeChannelBins(10000);
  ADQAParam->SetChargeChannelMin(0);
  ADQAParam->SetChargeChannelMax(10000);
  
  ADQAParam->SetNChargeSideBins(7000);
  ADQAParam->SetChargeSideMin(10);
  ADQAParam->SetChargeSideMax(70000);

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Michal Broz");
  md->SetComment("QA parameters for AD");
  //md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  //md->SetBeamPeriod(0);
  AliCDBId id("AD/Calib/QAParam", 0, AliCDBRunRange::Infinity());
  cdb->Put(ADQAParam, id, md);

  return;
}
