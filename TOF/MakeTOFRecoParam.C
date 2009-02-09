void MakeTOFRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  TObjArray *recoParamArray = new TObjArray();

  {
	  AliTOFRecoParam* tofRecoParam = AliTOFRecoParam::GetPbPbparam();
	  tofRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
	  recoParamArray->AddLast(tofRecoParam);
  }

  {
       	  AliTOFRecoParam* tofRecoParam = AliTOFRecoParam::GetPPparam();
	  tofRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
	  recoParamArray->AddLast(tofRecoParam);
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
  md->SetResponsible("Chiara Zampolli");
  md->SetComment("Reconstruction parameters TOF");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("TOF/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  man->GetDefaultStorage()->Put(recoParamArray,id, md);
}


