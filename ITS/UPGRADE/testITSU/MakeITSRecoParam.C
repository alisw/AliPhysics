void MakeITSRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult, const char* cdbURI="local://") {
//========================================================================
//
// Steering macro for ITS reconstruction parameters
//
// Author: A.Dainese
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================

  const char* macroname = "MakeITSRecoParam.C";
  //
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");
  //
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);

  int nLr = 7;
  
  TObjArray *recoParamArray = new TObjArray();
  //
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetCosmicTestParam();
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    itsRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    itsRecoParam->SetTitle("Cosmic");
    recoParamArray->AddLast(itsRecoParam);
  }
  //
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetLowFluxParam();
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    itsRecoParam->SetTitle("LowMult");
    recoParamArray->AddLast(itsRecoParam);
  }
  {
    AliITSURecoParam * itsRecoParam = AliITSURecoParam::GetHighFluxParam();
    //
    itsRecoParam->SetNLayers(nLr);
    //
    //******************************************************************
    itsRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    itsRecoParam->SetTitle("HighMult");
    recoParamArray->AddLast(itsRecoParam);
  }
  //
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
  md->SetResponsible("Andrea Dainese");
  md->SetComment("Reconstruction parameters ITS.");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);
  //
  return;
}

