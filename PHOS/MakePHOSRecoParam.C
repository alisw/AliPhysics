void MakePHOSRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kDefault)
{
  //========================================================================
  //
  // Steering macro for PHOS reconstruction parameters
  //
  // Author: Yuri Kharlov
  // 13.08.2008
  //========================================================================


  const char* macroname = "MakePHOSRecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  
  TObjArray *recoParamArray = new TObjArray();

  {
    AliPHOSRecoParam * phosRecoParam = AliPHOSRecoParam::GetDefaultParameters();
    phosRecoParam->SetEventSpecie(AliRecoParam::kDefault);
    phosRecoParam->SetEMCSubtractPedestals(kFALSE);
    phosRecoParam->SetEMCRawDigitThreshold(2);
    phosRecoParam->SetEMCMinE(0.012);
    phosRecoParam->SetEMCClusteringThreshold(0.20);
    phosRecoParam->SetEMCDecoderVersion("v1");
    phosRecoParam->SetEMCSampleQualityCut(10.);
    recoParamArray->AddLast(phosRecoParam);
  }
  {
    AliPHOSRecoParam * phosRecoParam = AliPHOSRecoParam::GetDefaultParameters();
    recoParamArray->AddLast(phosRecoParam);
    phosRecoParam->SetEMCClusteringThreshold(0.02008);
    phosRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
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
  md->SetComment("Reconstruction parameters PHOS");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("PHOS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);

  return;
}

