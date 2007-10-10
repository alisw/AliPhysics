void UpdateCDBVertexDiamond(Double_t xmed = 0., Double_t ymed = 0., Double_t sigx = 0.005, Double_t sigy = 0.005, Double_t sigz = 5.3) {
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
    metadata->SetResponsible("prino@to.infn.it");
    metadata->SetComment("Default mean vertex position");
    metadata->SetAliRootVersion(alirootv);
    metadata->SetComment(Form("Default trigger description produced with root version %s and AliRoot version %s",rootv,alirootv));
  }

  Printf("Storing in CDB the default trigger description produced with root version %s and AliRoot version %s",rootv,alirootv);

  Double_t resolx=35./10000.;
  Double_t resoly=35./10000.;
  Double_t sigma[3],position[3];
  position[0]=xmed;
  position[1]=ymed;
  position[2]=0.;
  sigma[0]=TMath::Sqrt(sigx*sigx+resolx*resolx);
  sigma[1]=TMath::Sqrt(sigy*sigy+resoly*resoly);
  sigma[2]=sigz;

  AliESDVertex *vertex = new AliESDVertex(position,sigma,"vtxmean");
  vertex->PrintStatus();

  man->Put(vertex,id,metadata);
}

