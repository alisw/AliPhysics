/// \file MakeTPCRecoParam.C
///
/// Steering macro for TPC reconstruction parameters
///
/// Contact: marian.ivanov@cern.ch
///
/// \author M.Ivanov


void MakeTPCRecoParam_Run2Dist(const char* storage="local://OCDB", AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult) {

  const char* macroname = "MakeTPCRecoParam.C";

  double sysCovMatColl[5] = {0.02,0.3, 2.5e-4, 2.5e-4, 5.0e-3}; // systematic error for low/high flux
  double sysCovMatCosm[5] = {0.05,0.5, 5.0e-4, 5.0e-4, 2.0e-2}; // systematic error for cosmics and calb (old one from Marian)
  // vector for charging-up masking
  TVectorF zreg(1),zregsigInv(1);
  zreg[0] = -5.;
  zregsigInv[0] = 1./2.0; // 2cm sigma in Z

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage(storage);
  //
  //
  TObjArray *recoParamArray = new TObjArray();  
  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
    tpcRecoParam->SetEventSpecie(AliRecoParam::kLowMult);

    //    Double_t sysError[5]={0.01,0.02, 0.01/150., 0.02/150.,0.01/(150*150.)};
    tpcRecoParam->SetSystematicError(sysCovMatColl);
    tpcRecoParam->SetMinMaxCutAbs(2.5);
    tpcRecoParam->SetMinLeftRightCutAbs(5.);
    tpcRecoParam->SetMinUpDownCutAbs(5.);
    //
    tpcRecoParam->SetMinMaxCutSigma(2.5);
    tpcRecoParam->SetMinLeftRightCutSigma(5.);
    tpcRecoParam->SetMinUpDownCutSigma(5.);
    //
    tpcRecoParam->SetUseHLTClusters(4);
    tpcRecoParam->SetUseOulierClusterFilter(1);
    tpcRecoParam->SetUseOnePadCluster(1);
    tpcRecoParam->SetUseHLTOnePadCluster(0);
    tpcRecoParam->SetCtgRange(10);
 
    // we use in data Ion tail and X-talk corrections
    tpcRecoParam->SetUseIonTailCorrection(1);
    tpcRecoParam->SetCrosstalkCorrection(1.);
    tpcRecoParam->SetCrosstalkCorrectionMissingCharge(1.);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(0);

    tpcRecoParam->SetUseGainCorrectionTime(1);

    // Settings needed for SP distortion corrections
    tpcRecoParam->SetUseCorrectionMap(kTRUE);
    tpcRecoParam->SetAccountDistortions(1);
    // scaling
    tpcRecoParam->SetUseLumiType(AliLumiTools::kLumiCTP);
    tpcRecoParam->SetCorrMapTimeDepMethod(AliTPCRecoParam::kCorrMapGlobalScalingLumi);
    // add difference to reference distortions as syst. error
    tpcRecoParam->SetUseDistortionFractionAsErrorY(-1);
    tpcRecoParam->SetUseDistortionFractionAsErrorZ(-1);
    //
    // add difference to reference dispersion as syst. error
    tpcRecoParam->SetUseDistDispFractionAsErrorY(0.6);
    tpcRecoParam->SetUseDistDispFractionAsErrorZ(0.3);
    //
    // eliminate pads with too large distortions XYZ or dispersions
    tpcRecoParam->SetBadPadMaxDistXYZD(7.0,10.,999, 1.0);
    
    // SP corrections are coupled to these VDrift settings
    tpcRecoParam->SetUseDriftCorrectionTime(1);
    tpcRecoParam->SetUseDriftCorrectionGY(1);

    // disable obsolete corrections
    tpcRecoParam->SetUseComposedCorrection(kFALSE);
    tpcRecoParam->SetUseExBCorrection(kFALSE);
    tpcRecoParam->SetUseAlignmentTime(kFALSE); 
    tpcRecoParam->SetUseFieldCorrection(0);
    tpcRecoParam->SetUseRPHICorrection(0);
    tpcRecoParam->SetUseRadialCorrection(0);
    tpcRecoParam->SetUseQuadrantAlignment(0);
    tpcRecoParam->SetUseSectorAlignment(0);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(kFALSE);
    tpcRecoParam->SetUseAlignmentTime(kFALSE);
    //
    // charging-up zone errors
    tpcRecoParam->SetSystErrClInnerRegZ( new TVectorF(zreg) );
    tpcRecoParam->SetSystErrClInnerRegZSigInv( new TVectorF(zregsigInv) );
    double *clInner = (double*)tpcRecoParam->GetSystematicErrorClusterInner();
    clInner[0] = 1.0; // 1.5 cm err at maximum
    clInner[1] = 5.0; // dumped with 5 cm
    //
    recoParamArray->AddLast(tpcRecoParam);
  }
  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetHighFluxParam();
    tpcRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
 
    //    Double_t sysError[5]={0.01,0.02, 0.01/150., 0.02/150.,0.01/(150*150.)};
    tpcRecoParam->SetSystematicError(sysCovMatColl);
    //
    tpcRecoParam->SetUseHLTClusters(4);
    tpcRecoParam->SetUseOulierClusterFilter(1);
    tpcRecoParam->SetUseOnePadCluster(1);
    tpcRecoParam->SetUseHLTOnePadCluster(0);
    tpcRecoParam->SetCtgRange(10);
 
    // we use in data Ion tail and X-talk corrections
    tpcRecoParam->SetUseIonTailCorrection(1);
    tpcRecoParam->SetCrosstalkCorrection(1.);
    tpcRecoParam->SetCrosstalkCorrectionMissingCharge(1.);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(0);

    tpcRecoParam->SetUseGainCorrectionTime(1);

    // Settings needed for SP distortion corrections
    tpcRecoParam->SetUseCorrectionMap(kTRUE);
    tpcRecoParam->SetAccountDistortions(1);
    // scaling
    tpcRecoParam->SetUseLumiType(AliLumiTools::kLumiCTP);
    tpcRecoParam->SetCorrMapTimeDepMethod(AliTPCRecoParam::kCorrMapGlobalScalingLumi);
    // add difference to reference distortions as syst. error
    tpcRecoParam->SetUseDistortionFractionAsErrorY(-1);
    tpcRecoParam->SetUseDistortionFractionAsErrorZ(-1);
    //
    // add difference to reference dispersion as syst. error
    tpcRecoParam->SetUseDistDispFractionAsErrorY(0.6);
    tpcRecoParam->SetUseDistDispFractionAsErrorZ(0.3);
    //
    // eliminate pads with too large distortions XYZ or dispersions
    tpcRecoParam->SetBadPadMaxDistXYZD(7.0,10.,999, 1.0);
    
    // SP corrections are coupled to these VDrift settings
    tpcRecoParam->SetUseDriftCorrectionTime(1);
    tpcRecoParam->SetUseDriftCorrectionGY(1);

    // disable obsolete corrections
    tpcRecoParam->SetUseComposedCorrection(kFALSE);
    tpcRecoParam->SetUseExBCorrection(kFALSE);
    tpcRecoParam->SetUseAlignmentTime(kFALSE); 
    tpcRecoParam->SetUseFieldCorrection(0);
    tpcRecoParam->SetUseRPHICorrection(0);
    tpcRecoParam->SetUseRadialCorrection(0);
    tpcRecoParam->SetUseQuadrantAlignment(0);
    tpcRecoParam->SetUseSectorAlignment(0);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(kFALSE);
    tpcRecoParam->SetUseAlignmentTime(kFALSE);
    //
    // charging-up zone errors
    tpcRecoParam->SetSystErrClInnerRegZ( new TVectorF(zreg) );
    tpcRecoParam->SetSystErrClInnerRegZSigInv( new TVectorF(zregsigInv) );
    double *clInner = (double*)tpcRecoParam->GetSystematicErrorClusterInner();
    clInner[0] = 1.0; // 1.5 cm err at maximum
    clInner[1] = 5.0; // dumped with 5 cm
    //    
    recoParamArray->AddLast(tpcRecoParam);
  }
  
  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLaserTestParam(kFALSE);
    tpcRecoParam->SetEventSpecie(AliRecoParam::kCalib);
    //
    tpcRecoParam->SetUseHLTClusters(4);
    tpcRecoParam->SetUseOulierClusterFilter(1);
    tpcRecoParam->SetUseOnePadCluster(1);
    tpcRecoParam->SetUseHLTOnePadCluster(1);
    tpcRecoParam->SetCtgRange(10);
 
    // we use in data Ion tail and X-talk corrections
    tpcRecoParam->SetUseIonTailCorrection(1);
    tpcRecoParam->SetCrosstalkCorrection(1.);
    tpcRecoParam->SetCrosstalkCorrectionMissingCharge(1.);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(0);
    //
    tpcRecoParam->SetSystematicError(sysCovMatCosm);
    //
    recoParamArray->AddLast(tpcRecoParam);
  }

  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kFALSE);
    tpcRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    //
    tpcRecoParam->SetTimeInterval(60,1000);
    tpcRecoParam->SetUseHLTClusters(4);
    tpcRecoParam->SetUseOulierClusterFilter(1);
    tpcRecoParam->SetUseOnePadCluster(1);
    tpcRecoParam->SetUseHLTOnePadCluster(0);
    tpcRecoParam->SetCtgRange(10);
 
    // in cosmics we use ComposedCorrection ? 
    tpcRecoParam->SetUseComposedCorrection(kTRUE);
    // disable obsolete corrections
    tpcRecoParam->SetUseExBCorrection(kFALSE);
    tpcRecoParam->SetUseAlignmentTime(kFALSE); 
    tpcRecoParam->SetUseFieldCorrection(0);
    tpcRecoParam->SetUseRPHICorrection(0);
    tpcRecoParam->SetUseRadialCorrection(0);
    tpcRecoParam->SetUseQuadrantAlignment(0);
    tpcRecoParam->SetUseSectorAlignment(0);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(kFALSE);
    tpcRecoParam->SetUseAlignmentTime(kFALSE);

    // we use in data Ion tail and X-talk corrections
    tpcRecoParam->SetUseIonTailCorrection(1);
    tpcRecoParam->SetCrosstalkCorrection(1.);
    tpcRecoParam->SetCrosstalkCorrectionMissingCharge(1.);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(0);

    tpcRecoParam->SetTimeInterval(60,940);
    //    Double_t sysError[5]={0.3,1, 0.3/150., 1./150.,0.3/(150*150.)};
    //    tpcRecoParam->SetSystematicError(sysError);
    tpcRecoParam->SetSystematicError(sysCovMatCosm);
    //
    tpcRecoParam->SetMinMaxCutAbs(2.5);
    tpcRecoParam->SetMinLeftRightCutAbs(5.);
    tpcRecoParam->SetMinUpDownCutAbs(5.);
    //
    tpcRecoParam->SetMinMaxCutSigma(2.5);
    tpcRecoParam->SetMinLeftRightCutSigma(5.);
    tpcRecoParam->SetMinUpDownCutSigma(5.);
    recoParamArray->AddLast(tpcRecoParam);
    //
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
  md->SetResponsible("Marian Ivanov");
  md->SetComment("Reconstruction parameters TPC, as used for Run2 data with distortions");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("TPC/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);
  return;
}

