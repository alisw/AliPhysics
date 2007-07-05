void UpdateCDBIdealGeom(){
  // produce the ideal geometry with the current AliRoot and store it in the
  // CDB

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  AliCDBId id("GRP/Geometry/Data",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *md= new AliCDBMetaData();

  // Get root version
  const char* rootv = gROOT->GetVersion();

  // Get AliRoot version from file to set it in the metadata of the entry
  TFile *fv= TFile::Open("$ALICE_ROOT/CVS/Repository?filetype=raw","READ");
  Int_t size = fv->GetSize();
  char *buf = new Char_t[size];
  memset(buf, '\0', size);
  fv->Seek(0);
  char* alirootv;
  if ( fv->ReadBuffer(buf, size) ) {
    Printf("Error reading AliRoot version from file to buffer!");
    alirootv = "";
  }
  if(buf=="AliRoot"){
    alirootv="HEAD";
  }else{
    alirootv = buf;
    md->SetAliRootVersion(alirootv);
    md->SetComment(Form("Geometry produced with root version %s and AliRoot version %s",rootv,alirootv));
  }
  
  gAlice->Init();
  
  if(!gGeoManager){
    Printf("Unable to produce a valid geometry to be put in the CDB!");
    return;
  }
  
  Printf("Storing in CDB geometry produced with root version %s and AliRoot version %s",rootv,alirootv);
  man->Put(gGeoManager,id,md);

}
