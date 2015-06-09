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
  ADQAParam->SetMaxPedDiff(1);
  ADQAParam->SetMaxPedWidth(1.5);
  ADQAParam->SetSatMed(0.1);
  ADQAParam->SetSatHigh(0.2);
  ADQAParam->SetSatHuge(0.5);
  
  ADQAParam->SetNTdcTimeBins(3062); 
  ADQAParam->SetTdcTimeMin(0.976562); 
  ADQAParam->SetTdcTimeMax(300.0); 
  
  ADQAParam->SetNChargeChannelBins(10000);
  ADQAParam->SetChargeChannelMin(1);
  ADQAParam->SetChargeChannelMax(10001);
  
  ADQAParam->SetNChargeSideBins(8000);
  ADQAParam->SetChargeSideMin(10);
  ADQAParam->SetChargeSideMax(80000);

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
