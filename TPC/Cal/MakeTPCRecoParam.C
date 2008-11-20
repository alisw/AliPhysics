void MakeTPCRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult) {
//========================================================================
//
// Steering macro for TPC reconstruction parameters
//
// Author: M.Ivanov
// Contact: marian.ivanov@cern.ch
//
//========================================================================


  const char* macroname = "MakeTPCRecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  //
  //
  TObjArray *recoParamArray = new TObjArray();  
  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
    Double_t sysError[5]={0.01,0.02, 0.01/150., 0.02/150.,0.01/(150*150.)};
    tpcRecoParam->SetSystematicError(sysError);
    tpcRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    tpcRecoParam->SetMinMaxCutAbs(2.5);
    tpcRecoParam->SetMinLeftRightCutAbs(5.);
    tpcRecoParam->SetMinUpDownCutAbs(5.);
    //
    tpcRecoParam->SetMinMaxCutSigma(2.5);
    tpcRecoParam->SetMinLeftRightCutSigma(5.);
    tpcRecoParam->SetMinUpDownCutSigma(5.);
    
    recoParamArray->AddLast(tpcRecoParam);
  }
  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetHighFluxParam();
    tpcRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    Double_t sysError[5]={0.01,0.02, 0.01/150., 0.02/150.,0.01/(150*150.)};
    tpcRecoParam->SetSystematicError(sysError);
    recoParamArray->AddLast(tpcRecoParam);
  }
  
  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLaserTestParam(kFALSE);
    tpcRecoParam->SetEventSpecie(AliRecoParam::kCalib);
    recoParamArray->AddLast(tpcRecoParam);
  }

  {
    AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kFALSE);
    tpcRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    tpcRecoParam->SetTimeInterval(60,940);
    Double_t sysError[5]={0.3,1, 0.3/150., 1./150.,0.3/(150*150.)};
    tpcRecoParam->SetSystematicError(sysError);
    //
    tpcRecoParam->SetMinMaxCutAbs(2.5);
    tpcRecoParam->SetMinLeftRightCutAbs(5.);
    tpcRecoParam->SetMinUpDownCutAbs(5.);
    //
    tpcRecoParam->SetMinMaxCutSigma(2.5);
    tpcRecoParam->SetMinLeftRightCutSigma(5.);
    tpcRecoParam->SetMinUpDownCutSigma(5.);
    recoParamArray->AddLast(tpcRecoParam);
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
  md->SetComment("Reconstruction parameters TPC");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("TPC/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);
  return;
}

