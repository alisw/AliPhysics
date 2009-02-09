void MakeZDCRecoParam(Int_t type=0){
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
  
  AliZDCRecoParam *zdcRecoParam = 0;
  //
  switch(type) {
  case 0:
    zdcRecoParam = (AliZDCRecoParampp*) AliZDCRecoParampp::GetppRecoParam();
    zdcRecoParam->SetBeamEnergy(5000.); 
    break;
  case 1:
    zdcRecoParam = (AliZDCRecoParamPbPb*) AliZDCRecoParamPbPb::GetPbPbRecoParam();
    break;
  case default:
    printf("Event type not implemented\n");
    return;
    break;
  }
  //
  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Chiara Oppedisano");
  md->SetComment("Reconstruction parameters for ZDC");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  md->SetObjectClassName("AliZDCRecoParam");
  AliCDBId id("ZDC/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  AliCDBManager::Instance()->GetDefaultStorage()->Put(zdcRecoParam,id, md);

}
