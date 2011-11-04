//____________________________________________________
void AliTRDmakeTrkDB(const Char_t *file)
{
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc) return;
  
  // Attach clusters likelihoods
  AliCDBMetaData *metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexandru Bercuci");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-30-02"); //root version
  metaData->SetComment(
    "Likelihoods for Attach Cluster.\n"
    " Tunned on p-p run 159580.");
  AliCDBId id("TRD/Calib/TrkAttach", 151536, AliCDBRunRange::Infinity());
  AliTRDCalTrkAttach attach;
  if(!attach.LoadReferences(file)) return;
//   attach.SetNsgmDy(Int_t ns0, Int_t ns1);
//   attach.SetLikeMinRelDecrease(Float_t p0, Float_t p1);
//   attach.SetRClikeLimit(Float_t rc);
  attach.SetScaleCov(5.);
  gStorLoc->Put(&attach, id, metaData); 

  return;
}

