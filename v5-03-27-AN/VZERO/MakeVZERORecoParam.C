void MakeVZERORecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult) {
//========================================================================
//
// Steering macro for VZERO reconstruction parameters
//
// Author: Brigitte Cheynis
//
//========================================================================

  const char* macroname = "MakeVZERORecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  
  TObjArray *recoParamArray = new TObjArray();

  {
    AliVZERORecoParam * vzeroRecoParam = new AliVZERORecoParam;
    vzeroRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    recoParamArray->AddLast(vzeroRecoParam);
  }
  {
    AliVZERORecoParam * vzeroRecoParam = new AliVZERORecoParam;
    // the following two settings are needed to high lumi runs in 2011
    vzeroRecoParam->SetStartClock(9);
    vzeroRecoParam->SetEndClock(11);
    vzeroRecoParam->SetNPreClocks(1);
    vzeroRecoParam->SetNPostClocks(1);
    vzeroRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    recoParamArray->AddLast(vzeroRecoParam);
  }
  {
    AliVZERORecoParam * vzeroRecoParam = new AliVZERORecoParam;
    vzeroRecoParam->SetStartClock(9);
    vzeroRecoParam->SetEndClock(11);
    vzeroRecoParam->SetNPreClocks(6);
    vzeroRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    recoParamArray->AddLast(vzeroRecoParam);
  }

  // Set the default
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i =0; i < recoParamArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam *param = (AliDetectorRecoParam *)recoParamArray->UncheckedAt(i);
    if (!param) continue;
    if (default & param->GetEventSpecie()) {
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
  md->SetResponsible("Brigitte Cheynis");
  md->SetComment("Reconstruction parameters for VZERO");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("VZERO/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}
