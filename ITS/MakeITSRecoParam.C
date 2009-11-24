void MakeITSRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult) {
//========================================================================
//
// Steering macro for ITS reconstruction parameters
//
// Author: A.Dainese
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================


  const char* macroname = "MakeITSRecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  
  TObjArray *recoParamArray = new TObjArray();

  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
    itsRecoParam->SetFactorSAWindowSizes(20);
    itsRecoParam->SetClusterErrorsParam(2);
    itsRecoParam->SetFindV0s(kFALSE);
    itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
    itsRecoParam->SetUseAmplitudeInfo(kFALSE);
    // In case we want to switch off a layer
    //  itsRecoParam->SetLayerToSkip(<N>);
    //  itsRecoParam->SetLayerToSkip(4);
    //  itsRecoParam->SetLayerToSkip(5);
    //    itsRecoParam->SetLayerToSkip(2);
    //    itsRecoParam->SetLayerToSkip(3);
    //itsRecoParam->SetSAOnePointTracks();
    itsRecoParam->SetClusterMisalError(0.1); // [cm]
    itsRecoParam->SetSAUseAllClusters();
    itsRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    recoParamArray->AddLast(itsRecoParam);
  }
  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetLowFluxParam();
    itsRecoParam->SetClusterErrorsParam(2);
    itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    // Plane Efficiency evaluation with tracklets Method
    itsRecoParam->SetIPlanePlaneEff(-1);
    itsRecoParam->SetComputePlaneEff(kTRUE,kFALSE);
    itsRecoParam->SetUseTrackletsPlaneEff(kTRUE);
    itsRecoParam->SetTrackleterPhiWindowL2(0.07);
    itsRecoParam->SetTrackleterZetaWindowL2(0.4);
    itsRecoParam->SetTrackleterPhiWindowL1(0.10);
    itsRecoParam->SetTrackleterZetaWindowL1(0.6);
    itsRecoParam->SetUpdateOncePerEventPlaneEff(kTRUE);
    itsRecoParam->SetMinContVtxPlaneEff(3);
    itsRecoParam->SetOptTrackletsPlaneEff(kTRUE);
    recoParamArray->AddLast(itsRecoParam);
  }
  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetHighFluxParam();
    itsRecoParam->SetClusterErrorsParam(2);
    itsRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    recoParamArray->AddLast(itsRecoParam);
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
  md->SetResponsible("Andrea Dainese");
  md->SetComment("Reconstruction parameters ITS");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}

