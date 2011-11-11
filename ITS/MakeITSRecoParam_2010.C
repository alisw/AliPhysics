void MakeITSRecoParam_2010(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult, const char* cdbURI="local://") {
//========================================================================
//
// Steering macro for ITS reconstruction parameters
//
// Author: A.Dainese
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================


  const char* macroname = "MakeITSRecoParam_2010.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);
  
  TObjArray *recoParamArray = new TObjArray();

  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
    // find independently ITS SA tracks 
    itsRecoParam->SetSAUseAllClusters();
    itsRecoParam->SetMinNPointsSA(2);
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

    itsRecoParam->SetClusterErrorsParam(2);
    itsRecoParam->SetFindV0s(kFALSE);
    itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
    itsRecoParam->SetUseAmplitudeInfo(kFALSE);
    itsRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    itsRecoParam->SetTitle("Cosmic");
    recoParamArray->AddLast(itsRecoParam);
  }
  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetLowFluxParam();
    itsRecoParam->SetClusterErrorsParam(2);

    // find independently ITS SA tracks 
    itsRecoParam->SetSAUseAllClusters();
    itsRecoParam->SetOuterStartLayerSA(2);

    itsRecoParam->SetAllowProlongationWithEmptyRoad(kTRUE);
    
    // larger seach windows for SA (in case of large misalignments)
    itsRecoParam->SetFactorSAWindowSizes(2);
    
    // Misalignment syst errors decided at ITS meeting 25.03.2010
    // Errors in Z reduced on 11.10.2010 for SPD and SDD
    // additional error due to misal (B off)
    itsRecoParam->SetClusterMisalErrorY(0.0010,0.0010,0.0300,0.0300,0.0020,0.0020); // [cm]
    itsRecoParam->SetClusterMisalErrorZ(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000); // [cm]
    // additional error due to misal (B on)
    itsRecoParam->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020); // [cm]
    itsRecoParam->SetClusterMisalErrorZBOn(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000); // [cm]
    //----

    // SDD configuration 
    itsRecoParam->SetUseSDDCorrectionMaps(kTRUE); // changed 30.04.2010
    itsRecoParam->SetUseSDDClusterSizeSelection(kTRUE);
    itsRecoParam->SetMinClusterChargeSDD(30.);
    itsRecoParam->SetUseUnfoldingInClusterFinderSDD(kFALSE);

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
    // 
    itsRecoParam->SetTrackleterRemoveClustersFromOverlaps(kTRUE);
    itsRecoParam->SetTrackleterPhiWindow(0.08);
    itsRecoParam->SetTrackleterThetaWindow(0.025);
    itsRecoParam->SetTrackleterScaleDThetaBySin2T(kTRUE);
    //
    // V0 finder (A. Marin)
    itsRecoParam->GetESDV0Params()->SetMaxPidProbPionForb(0.9);

    //******************************************************************

    itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    itsRecoParam->SetTitle("LowMult");
    recoParamArray->AddLast(itsRecoParam);
  }
  {
    AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetHighFluxParam();
    itsRecoParam->SetClusterErrorsParam(2);

    // find independently ITS SA tracks for nContrSPD<50
    itsRecoParam->SetSAUseAllClusters();
    itsRecoParam->SetMaxSPDcontrForSAToUseAllClusters(50);

    itsRecoParam->SetImproveWithVertex(kTRUE);
    // Misalignment syst errors decided at ITS meeting 25.03.2010
    // additional error due to misal (B off)
    itsRecoParam->SetClusterMisalErrorY(0.0010,0.0010,0.0300,0.0300,0.0020,0.0020); // [cm]
    itsRecoParam->SetClusterMisalErrorZ(0.0100,0.0100,0.0100,0.0100,0.0500,0.0500); // [cm]
    // additional error due to misal (B on)
    itsRecoParam->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020); // [cm]
    itsRecoParam->SetClusterMisalErrorZBOn(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000); // [cm]
    //----

    //Vertexer Z
    itsRecoParam->SetVertexerZ();


    // tracklets
    itsRecoParam->SetTrackleterPhiWindowL2(0.07);
    itsRecoParam->SetTrackleterZetaWindowL2(0.4);
    itsRecoParam->SetTrackleterPhiWindowL1(0.10);
    itsRecoParam->SetTrackleterZetaWindowL1(0.6);
    //
    itsRecoParam->SetTrackleterPhiWindow(0.06);
    itsRecoParam->SetTrackleterThetaWindow(0.025);
    itsRecoParam->SetTrackleterScaleDThetaBySin2T(kTRUE);
    //
    // Removal of tracklets reconstructed in the SPD overlaps 
    itsRecoParam->SetTrackleterRemoveClustersFromOverlaps(kTRUE);

    // SDD configuration 
    itsRecoParam->SetUseSDDCorrectionMaps(kTRUE); 
    itsRecoParam->SetUseSDDClusterSizeSelection(kTRUE);
    itsRecoParam->SetMinClusterChargeSDD(30.);
    itsRecoParam->SetUseUnfoldingInClusterFinderSDD(kFALSE);

    itsRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    itsRecoParam->SetTitle("HighMult");
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
  md->SetComment("Reconstruction parameters ITS.");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}

