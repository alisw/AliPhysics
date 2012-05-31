void MakeZDCRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult){
//========================================================================
//
// Steering macro to create and store in OCDB
//       ZDC reconstruction parameters
//
// Contact: chiara.oppedisano@to.infn.it
//
//========================================================================

  AliCDBManager* cdb = AliCDBManager::Instance();
  //if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  TObjArray *recoParamArray = new TObjArray();
  
 
  AliZDCRecoParampp* zdcppRecoParam = AliZDCRecoParampp::GetLowFluxParam();
  zdcppRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
  recoParamArray->AddLast(zdcppRecoParam);

  AliZDCRecoParamPbPb* zdcAARecoParam = AliZDCRecoParamPbPb::GetHighFluxParam();
  zdcAARecoParam->SetEventSpecie(AliRecoParam::kHighMult);
  recoParamArray->AddLast(zdcAARecoParam);
  
  // Set the default
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i =0; i<recoParamArray->GetEntriesFast(); i++) {
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

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Chiara Oppedisano");
  md->SetComment("Reconstruction parameters for ZDC");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("ZDC/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recoParamArray,id, md);
 
}
