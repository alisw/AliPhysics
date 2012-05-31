void MakePHOSRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult)
{
  //========================================================================
  //
  // Steering macro for PHOS reconstruction parameters
  //
  // Author: Yuri Kharlov
  // 01.11.2009
  //========================================================================
  // $Id$ */

  const char* macroname = "MakePHOSRecoParam.C";
  const Int_t firstRun = 0;

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  
  TObjArray *recoParamArray = new TObjArray();
  
  {
    // Reconstruction parameters for cosmic run 2009
    AliPHOSRecoParam * phosRecoParam = AliPHOSRecoParam::GetDefaultParameters();
    phosRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    phosRecoParam->SetEMCSubtractPedestals(kFALSE);
    phosRecoParam->SetEMCMinE(0.025);
    phosRecoParam->SetEMCClusteringThreshold(0.05);
    phosRecoParam->SetEMCFitterVersion("v0");
    phosRecoParam->SetEMCUnfolding(kFALSE);
    phosRecoParam->SetEMCEnergyCorrectionOn(kFALSE);
    phosRecoParam->SetName("Cosmic2009");
    recoParamArray->AddLast(phosRecoParam);
  }
  {
    // Reconstruction parameters for the first pp run 2009
    AliPHOSRecoParam * phosRecoParam = AliPHOSRecoParam::GetDefaultParameters();
    phosRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    phosRecoParam->SetEMCSubtractPedestals(kFALSE);
    phosRecoParam->SetEMCRawDigitThreshold(2);
    phosRecoParam->SetEMCMinE(0.012);
    phosRecoParam->SetEMCClusteringThreshold(0.22);
    phosRecoParam->SetEMCFitterVersion("v0");
    phosRecoParam->SetEMCUnfolding(kTRUE);
    phosRecoParam->SetEMCSampleQualityCut(10.);
    phosRecoParam->SetName("LowMult2009");
    recoParamArray->AddLast(phosRecoParam);
  }
  {
    // Reconstruction parameters for the first PbPb run 2010
    AliPHOSRecoParam * phosRecoParam = AliPHOSRecoParam::GetDefaultParameters();
    phosRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    phosRecoParam->SetEMCSubtractPedestals(kFALSE);
    phosRecoParam->SetEMCRawDigitThreshold(2);
    phosRecoParam->SetEMCMinE(0.015);
    phosRecoParam->SetEMCClusteringThreshold(0.40);
    phosRecoParam->SetEMCFitterVersion("v0");
    phosRecoParam->SetEMCUnfolding(kTRUE);
    phosRecoParam->SetName("HighMult2009");
    recoParamArray->AddLast(phosRecoParam);
  }
  { // Reconstruction parameters for "calibration" events
    AliPHOSRecoParam * phosRecoParam = AliPHOSRecoParam::GetDefaultParameters();
    phosRecoParam->SetEventSpecie(AliRecoParam::kCalib);
    phosRecoParam->SetEMCSubtractPedestals(kFALSE);
    phosRecoParam->SetEMCMinE(0.025);
    phosRecoParam->SetEMCClusteringThreshold(0.05);
    phosRecoParam->SetEMCFitterVersion("v0");
    phosRecoParam->SetEMCUnfolding(kFALSE);
    phosRecoParam->SetEMCEnergyCorrectionOn(kFALSE);
    phosRecoParam->SetName("Calib2009");
    recoParamArray->AddLast(phosRecoParam);
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
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("Reconstruction parameters PHOS for cosmic, lowMult and highMult");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("PHOS/Calib/RecoParam",firstRun,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}
