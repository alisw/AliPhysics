void MakePMDRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult){
  // Create PMD Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  TObjArray *recoParamArray = new TObjArray();

  {
    AliPMDRecoParam* pmdRecoParam = AliPMDRecoParam::GetPbPbParam();
    pmdRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    recoParamArray->AddLast(pmdRecoParam);
  }

  {
    AliPMDRecoParam* pmdRecoParam = AliPMDRecoParam::GetPPParam();
    pmdRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    recoParamArray->AddLast(pmdRecoParam);
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

  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Basanta Nandi");
  md->SetComment("Reconstruction parameters PMD");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("PMD/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  man->GetDefaultStorage()->Put(recoParamArray,id, md);
}


