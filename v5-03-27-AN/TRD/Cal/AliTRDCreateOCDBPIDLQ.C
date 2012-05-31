void AliTRDCreateOCDBPIDLQ(const char *fn){

  TObjArray *content = new TObjArray;
  TFile *in = TFile::Open(fn);
  TKey *key = NULL;
  TObject *tmp = NULL;
  TIter iter(in->GetListOfKeys());
  while((key = (TKey *)iter())){
    tmp = key->ReadObj();
    printf("Putting %s into the OCDB\n", tmp->GetName());
    content->Add(tmp);
  }

  AliCDBMetaData *metaData= new AliCDBMetaData(); 
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Markus Fasel");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-21-01"); //root version
  metaData->SetComment("TRD PID Reference Histos for the 1D Likelihood method");
  
  AliCDBId id("TRD/Calib/PIDLQ1D", 0, AliCDBRunRange::Infinity()); 
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *gStorLoc = man->GetStorage("local:///u/mfasel/OCDB");//$ALICE_ROOT/OCDB");
  if (!gStorLoc) {
    return;
  }
  gStorLoc->Put(content, id, metaData); 
  in->Close();

  return;
}