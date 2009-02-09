//____________________________________________________
void AliTRDmakeRecoParam()
{
  AliCDBMetaData *metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexandru Bercuci");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-21-01"); //root version
  metaData->SetComment("Ideal reconstruction parameters for low, high and cosmic runs");
  
  AliCDBId id("TRD/Calib/RecoParam", 0, AliCDBRunRange::Infinity()); 
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) {
    return;
  }
  gStorLoc->Put(CreateRecoParamObject(), id, metaData); 

  return;
}


//____________________________________________________
TObjArray* CreateRecoParamObject()
{
  TObjArray *recos = new TObjArray(4);

  AliTRDrecoParam *rec = 0x0;
  recos->AddLast(rec = AliTRDrecoParam::GetLowFluxParam());
  rec->SetAsDefault();
  rec->SetNameTitle("LOW", "TRD Low Flux Reco Param");
  // further settings for low flux reco param
  // reco->SetThisAndThat()

  recos->AddLast(rec = AliTRDrecoParam::GetLowFluxParam());
  rec->SetNameTitle("HLT", "TRD HLT Reco Param");
  rec->SetChi2Y(.1);
  rec->SetChi2Z(5.);

  recos->AddLast(rec = AliTRDrecoParam::GetHighFluxParam());
  rec->SetNameTitle("HIGH", "TRD High Flux Reco Param");

  recos->AddLast(rec = AliTRDrecoParam::GetCosmicTestParam());
  rec->SetNameTitle("COSMIC", "TRD Cosmic Reco Param");

  return recos;
}
