void MakeTRDFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TRD
  // Expects to read objects for FRAME
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjParams",1000);
  const char* macroname = "MakeTRDFullMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;
  TString Storage;
  
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  }else{
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
  }    
		  
  // load FRAME full misalignment objects (if needed, the macro
  // for FRAME has to be ran in advance) and apply them to geometry
  AliCDBPath fpath("GRP","Align","Data");
  AliCDBEntry *eFrame;
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    Info(macroname,"Loading FRAME alignment objects from CDB storage %s",
      Storage.Data());
    eFrame = storage->Get(fpath.GetPath(),cdb->GetRun());
  }else{
    eFrame = cdb->Get(fpath.GetPath());
  }
  if(!eFrame) Fatal(macroname,"Could not get the specified CDB entry!");
  TClonesArray* arFrame = (TClonesArray*) eFrame->GetObject();
  arFrame->Sort();
  Int_t nvols = arFrame->GetEntriesFast();
  Bool_t flag = kTRUE;
  for(Int_t j=0; j<nvols; j++)
  {
    AliAlignObj* alobj = (AliAlignObj*) arFrame->UncheckedAt(j);
    if (alobj->ApplyToGeometry() == kFALSE) flag = kFALSE;
  }
  if(!flag) Fatal(macroname,"Error in the application of FRAME objects");

   
  // sigmas for the supermodules
  Double_t smdx=0.3; // 3 mm
  Double_t smdy=0.3; // 3 mm
  Double_t smdz=0.3; // 3 mm
  Double_t smrx=0.4/1000/TMath::Pi()*180; // 0.4 mrad
  Double_t smry=2.0/1000/TMath::Pi()*180; // 2 mrad
  Double_t smrz=0.4/1000/TMath::Pi()*180; // 0.4 mrad


  // sigmas for the chambers
  Double_t chdx=0.1; // 1 mm
  Double_t chdy=0.1; // 1 mm
  Double_t chdz=0.1; // 1 mm
  Double_t chrx=1.0/1000/TMath::Pi()*180; // 1 mrad
  Double_t chry=1.0/1000/TMath::Pi()*180; // 1 mrad
  Double_t chrz=0.7/1000/TMath::Pi()*180; // 0.7 mrad

  Int_t sActive[18]={1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1};
  Double_t dx,dy,dz,rx,ry,rz;

  Int_t j=0;
  TRandom *ran = new TRandom(4357);
  UShort_t volid;
  const char *symname;

  // create the supermodules' alignment objects
  for (int iSect; iSect<18; iSect++) {
    TString sm_symname(Form("TRD/sm%02d",iSect));
    ran->Rannor(dx,rx);
    ran->Rannor(dy,ry);
    ran->Rannor(dz,rz);
    dx*=smdx;
    dy*=smdy;
    dz*=smdz;
    rx*=smrx;
    ry*=smry;
    rz*=smrz;
    if( (TString(gSystem->Getenv("REALSETUP")) == TString("kTRUE")) && !sActive[iSect] ) continue;
    new((*array)[j++]) AliAlignObjParams(sm_symname.Data(),0,dx,dy,dz,rx,ry,rz,kFALSE);
  }
  // apply supermodules' alignment objects
  Int_t smCounter=0;
  for(Int_t iSect=0; iSect<18; iSect++){
    if( (TString(gSystem->Getenv("REALSETUP")) == TString("kTRUE")) && !sActive[iSect] ) continue;
    AliAlignObjParams* smobj =
      (AliAlignObjParams*)array->UncheckedAt(smCounter++);
    if(!smobj->ApplyToGeometry()){
      Fatal(macroname,Form("application of full misalignment object for sector %d failed!",iSect));
      return;
    }
  }

  // create the chambers' alignment objects
  ran = new TRandom(4357);
  Int_t chId;
  for (Int_t iLayer = AliGeomManager::kTRD1; iLayer <= AliGeomManager::kTRD6; iLayer++) {
    chId=-1;
    for (Int_t iSect = 0; iSect < 18; iSect++){
      for (Int_t iCh = 0; iCh < 5; iCh++) {
      ran->Rannor(dx,rx);
      ran->Rannor(dy,ry);
      ran->Rannor(dz,rz);
      dx*=chdx;
      dy*=chdy;
      dz*=chdz;
      rx*=chrx;
      ry*=chry;
      rz*=chrz;
      chId++;
      if ((iSect==13 || iSect==14 || iSect==15) && iCh==2) continue;
      volid = AliGeomManager::LayerToVolUID(iLayer,chId);
      if( (TString(gSystem->Getenv("REALSETUP")) == TString("kTRUE")) && !sActive[iSect] ) continue;
      symname = AliGeomManager::SymName(volid);
      new((*array)[j++]) AliAlignObjParams(symname,volid,dx,dy,dz,rx,ry,rz,kFALSE);
    }
  }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "TRDfullMisalignment.root";
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
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Full misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();
}

