void MakeADRecoParamEntry(AliRecoParam::EventSpecie_t defaultEventSpecie=AliRecoParam::kLowMult, const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB") {
//========================================================================
//
// Steering macro for AD reconstruction parameters
//
// Author: Brigitte Cheynis
//
//========================================================================

  const char* macroname = "MakeADRecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(outputCDB);
  cdb->SetRun(0);
  
  TObjArray *recoParamArray = new TObjArray();

  {
    AliADRecoParam * ADRecoParam = new AliADRecoParam;
    ADRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    recoParamArray->AddLast(ADRecoParam);
  }
  {
    AliADRecoParam * ADRecoParam = new AliADRecoParam;
    // the following two settings are needed to high lumi runs in 2011
    ADRecoParam->SetStartClock(9);
    ADRecoParam->SetEndClock(11);
    ADRecoParam->SetNPreClocks(1);
    ADRecoParam->SetNPostClocks(1);
    ADRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    recoParamArray->AddLast(ADRecoParam);
  }
  {
    AliADRecoParam * ADRecoParam = new AliADRecoParam;
    ADRecoParam->SetStartClock(9);
    ADRecoParam->SetEndClock(11);
    ADRecoParam->SetNPreClocks(6);
    ADRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    recoParamArray->AddLast(ADRecoParam);
  }

  // Set the defaultEventSpecie
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i =0; i < recoParamArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam *param = (AliDetectorRecoParam *)recoParamArray->UncheckedAt(i);
    if (!param) continue;
    if (defaultEventSpecie & param->GetEventSpecie()) {
      param->SetAsDefault();
      defaultIsSet = kTRUE;
    }
  }

  if (!defaultIsSet) {
    Error(macroname,"The default reconstruction parameters are not set! Exiting...");
    return;
  }

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Michal Broz");
  md->SetComment("Reconstruction parameters for AD");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("AD/Calib/RecoParam", 0, AliCDBRunRange::Infinity());
  cdb->Put(recoParamArray, id, md);

  return;
}
