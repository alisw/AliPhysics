void UpdateCDBVertexDiamond() {
  // produce the trigger descriptorwith the current AliRoot and store it in the
  // CDB
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  AliCDBId id("GRP/Calib/MeanVertex",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *metadata= new AliCDBMetaData();

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
    metadata->SetResponsible("Tapan Nayak");
    metadata->SetAliRootVersion(alirootv);
    metadata->SetComment(Form("Default trigger description produced with root version %s and AliRoot version %s",rootv,alirootv));
  }

  Printf("Storing in CDB the default trigger description produced with root version %s and AliRoot version %s",rootv,alirootv);

  Double_t position[3] = {0.0,0.0,0.0};
  Double_t sigma[3] = {0.0,0.0,0.0};
  AliESDVertex *vertex = new AliESDVertex(position,sigma,"Default");
  vertex->Print();

  man->Put(vertex,id,metadata);
}

