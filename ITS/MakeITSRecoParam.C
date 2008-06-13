void MakeITSRecoParam(Int_t type=1) {
//========================================================================
//
// Steering macro for ITS reconstruction parameters
//
// Author: A.Dainese
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================


  const char* macroname = "MakeITSRecoParam.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://OCDB");
  
  AliITSRecoParam *itsRecoParam = 0;
  switch(type) {
  case 0:
    itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
    break;
  case 1:
    itsRecoParam = AliITSRecoParam::GetLowFluxParam();
    break;
  case 2:
    itsRecoParam = AliITSRecoParam::GetHighFluxParam();
    break;
  case default:
    printf("Wrong event type\n");
    return;
    break;
  }
  itsRecoParam->SetClusterErrorsParam(2);
  //itsRecoParam->SetClusterMisalError(1.0); // [cm]
  /*
  itsRecoParam->SetFindV0s(kTRUE);
  itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
  itsRecoParam->SetLayerToSkip(0);
  itsRecoParam->SetLayerToSkip(1);
  itsRecoParam->SetLayerToSkip(2);
  itsRecoParam->SetLayerToSkip(3);
  itsRecoParam->SetLayerToSkip(4);
  itsRecoParam->SetLayerToSkip(5);
  */

  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Andrea Dainese");
  md->SetComment("Reconstruction parameters ITS");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  AliCDBManager::Instance()->GetDefaultStorage()->Put(itsRecoParam,id, md);


  return;
}
