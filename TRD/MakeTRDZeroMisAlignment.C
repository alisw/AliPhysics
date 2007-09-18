void MakeTRDZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TRD
  //
  const char* macroname = "MakeTRDZeroMisAlignment.C";
  TClonesArray *array = new TClonesArray("AliAlignObjParams",1000);
  TClonesArray &alobj = *array;
   
  Int_t sActive[18]={0,0,1,1,1,0,1,0,0,0,0,1,1,0,1,1,0,0};
  Double_t dx=0.,dy=0.,dz=0.,rx=0.,ry=0.,rz=0.;

  Int_t j=0;
  UShort_t volid;
  const char *symname;

  // create the supermodules' alignment objects
  for (Int_t iSect; iSect<18; iSect++) {
    TString sm_symname(Form("TRD/sm%02d",iSect));
    if( (TString(gSystem->Getenv("PARTGEOM")) == TString("kTRUE")) && !sActive[iSect] ) continue;
    new(alobj[j++]) AliAlignObjParams(sm_symname.Data(),0,dx,dy,dz,rx,ry,rz,kTRUE);
  }
  
   // create the chambers' alignment objects
  Int_t chId;
  for (Int_t iLayer = AliGeomManager::kTRD1; iLayer <= AliGeomManager::kTRD6; iLayer++) {
    chId=-1;
    for (Int_t iSect = 0; iSect < 18; iSect++){
      for (Int_t iCh = 0; iCh < 5; iCh++) {
	chId++;
	volid = AliGeomManager::LayerToVolUID(iLayer,chId);
      symname = AliGeomManager::SymName(volid);
	if( (TString(gSystem->Getenv("PARTGEOM")) == TString("kTRUE")) && !sActive[iSect] ) continue;
      new(alobj[j++]) AliAlignObjParams(symname,volid,dx,dy,dz,rx,ry,rz,kTRUE);
    }
  }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "TRDzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"TRDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Zero misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}


