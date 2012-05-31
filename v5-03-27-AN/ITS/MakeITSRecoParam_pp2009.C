void MakeITSRecoParam_pp2009(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult, const char* cdbURI) {
//========================================================================
//
// Steering macro for ITS reconstruction parameters
//
// Author: A.Dainese
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================


  const char* macroname = "MakeITSRecoParam_pp2009.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);
  
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
    itsRecoParam->SetMinNPointsSA(2);
    itsRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    recoParamArray->AddLast(itsRecoParam);
  }
  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetLowFluxParam();
    itsRecoParam->SetClusterErrorsParam(2);
    //****** FIRST PHYSICS 2009 *********************

    /* //  First Collisions 23.11.2009 (same as COSMICS 2009) 
    // find independently ITS SA tracks 
    itsRecoParam->SetSAUseAllClusters();
    itsRecoParam->SetOuterStartLayerSA(AliITSgeomTGeo::GetNLayers()-2);

    // to maximize efficiency
    itsRecoParam->SetAllowProlongationWithEmptyRoad();
    
    // larger seach windows for SA (in case of large misalignments)
    itsRecoParam->SetNLoopsSA(33);
    itsRecoParam->SetFactorSAWindowSizes(20);
    
    // additional error due to misal (B off)
    itsRecoParam->SetClusterMisalErrorY(1.0,1.0,1.0,1.0,1.0,1.0); // [cm]
    itsRecoParam->SetClusterMisalErrorZ(1.0,1.0,1.0,1.0,1.0,1.0); // [cm]
    // additional error due to misal (B on)
    itsRecoParam->SetClusterMisalErrorYBOn(0.0,0.0,0.1,0.1,0.1,0.1); // [cm]
    itsRecoParam->SetClusterMisalErrorZBOn(0.1,0.1,0.1,0.1,0.1,0.1); // [cm]
    */

    //----  Collisions since 03.12.2009 
    // find independently ITS SA tracks 
    itsRecoParam->SetSAUseAllClusters();
    itsRecoParam->SetOuterStartLayerSA(2);

    // to maximize efficiency (change to kTRUE on 15.12.2009)
    //itsRecoParam->SetAllowProlongationWithEmptyRoad(kFALSE);
    itsRecoParam->SetAllowProlongationWithEmptyRoad(kTRUE);
    
    // larger seach windows for SA (in case of large misalignments)
    itsRecoParam->SetFactorSAWindowSizes(2);
    
    // additional error due to misal (B off)
    itsRecoParam->SetClusterMisalErrorY(0.01,0.01,0.1,0.1,0.1,0.1); // [cm]
    itsRecoParam->SetClusterMisalErrorZ(0.01,0.01,0.1,0.1,0.1,0.1); // [cm]
    // additional error due to misal (B on)
    itsRecoParam->SetClusterMisalErrorYBOn(0.01,0.01,0.1,0.1,0.1,0.1); // [cm]
    itsRecoParam->SetClusterMisalErrorZBOn(0.01,0.01,0.1,0.1,0.1,0.1); // [cm]
    //----

    // SDD configuration 
    itsRecoParam->SetUseSDDCorrectionMaps(kFALSE);
    itsRecoParam->SetUseSDDClusterSizeSelection(kTRUE);
    itsRecoParam->SetMinClusterChargeSDD(30.);

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
    // itsRecoParam->SetOptTrackletsPlaneEff(kTRUE); // activate it for MC (very important) !
    // Removal of tracklets reconstructed in the SPD overlaps 
    itsRecoParam->SetTrackleterRemoveClustersFromOverlaps(kTRUE);
  
    //******************************************************************

    itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
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
  md->SetComment("Reconstruction parameters ITS. Use large misal errors for cosmics and lowflux, 18 Nov 2009");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}


