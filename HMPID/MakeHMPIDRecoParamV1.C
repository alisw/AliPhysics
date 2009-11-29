void MakeHMPIDRecoParamV1(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult)
{
  //
  // Create HMPID Reco Parameters in OCDB
  //  
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  TObjArray *recoParamArray = new TObjArray();
  
  {
	  AliHMPIDRecoParamV1* hmpRecoParam = AliHMPIDRecoParamV1::GetLowFluxParam();
	  hmpRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
          hmpRecoParam->SetAsDefault();
	  recoParamArray->AddLast(hmpRecoParam);
          hmpRecoParam->PrintParameters();
  }
  {
	  AliHMPIDRecoParamV1* hmpRecoParam = AliHMPIDRecoParamV1::GetHighFluxParam();
	  hmpRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
	  recoParamArray->AddLast(hmpRecoParam);
          hmpRecoParam->PrintParameters();

  }
  {
	  AliHMPIDRecoParamV1* hmpRecoParam = AliHMPIDRecoParamV1::GetCosmicParam();
	  hmpRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
	  recoParamArray->AddLast(hmpRecoParam);
          hmpRecoParam->PrintParameters();
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
    Error("The default reconstruction parameters are not set! Exiting...");
    return;
  }
  
 
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Domenico DiBari");
  md->SetComment("Reconstruction parameters of HMPID, v1");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("HMPID/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  man->GetDefaultStorage()->Put(recoParamArray,id, md);
  
  return;
}
