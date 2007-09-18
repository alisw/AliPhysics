void MakeTOFFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TOF
  // Expects to read objects for FRAME
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjParams",2000);
  const char* macroname = "MakeTOFFullMisAlignment.C";
  
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
  // for FRAME has to be run in advance) and apply them to geometry
  AliCDBPath fpath("GRP","Align","Data");
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    Info(macroname,"Loading FRAME alignment objects from CDB storage %s",
	Storage.Data());
    AliCDBEntry *eFrame = storage->Get(fpath.GetPath(),cdb->GetRun());
  }else{
    AliCDBEntry *eFrame = cdb->Get(fpath.GetPath());
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
  if(!flag) Fatal(macroname,"Error in the application of FRAME alignment objects");

  //Produce objects for TOF supermodules
  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy vol id 

  Int_t nSMTOF = 18;
  Int_t sActive[18]={0,1,1,0,0,0,1,1,0,1,1,1,1,0,0,1,1,1};
  Int_t j=0;
  Double_t smdx, smdy, smdz=0., dpsi=0., dtheta, dphi=0.;
  TRandom *rnd   = new TRandom(2345);
  Double_t sigmatr = 0.4; // max shift in cm w.r.t. local ideal RS
  Double_t sigmarot = 0.06; // max rot in deg w.r.t. local ideal RS (~ 1 mrad)
  
  for(Int_t iSect=0; iSect<nSMTOF; iSect++) {
    TString symname(Form("TOF/sm%02d",iSect));
    smdx = rnd->Gaus(0.,sigmatr);
    smdy = rnd->Gaus(0.,sigmatr);
    dtheta = rnd->Gaus(0.,sigmarot);
    if( (TString(gSystem->Getenv("PARTGEOM")) == TString("kTRUE")) && !sActive[iSect] ) continue;
    new((*array)[j++]) AliAlignObjParams(symname.Data(), dvoluid, smdx, smdy, smdz, dpsi, dtheta, dphi, kFALSE);
  }
  // Apply objects for TOF supermodules 
  Int_t smCounter=0;
  for(Int_t iSect=0; iSect<nSMTOF; iSect++){
    if( (TString(gSystem->Getenv("PARTGEOM")) == TString("kTRUE")) && !sActive[iSect] ) continue;
    AliAlignObjParams* smobj = (AliAlignObjParams*)array->UncheckedAt(smCounter++);
    Info(macroname,Form("Applying object for sector %d ",iSect));
    if(!smobj->ApplyToGeometry()){
      Fatal(macroname,Form("application of full misalignment object for sector %d failed!",iSect));
      return;
    }
  }

  //Produce objects for TOF strips (same sigmas as for residual misalignment)
  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;
  Int_t strId=-1;

  Double_t sdx=0., sdy=0., sdz=0., sdpsi=0., sdtheta=0., sdphi=0.;
  TRandom *rnds   = new TRandom(4357);
  sigmatr = 0.1; // max shift in cm w.r.t. local ideal RS

  Int_t nstrA=15;
  Int_t nstrB=19;
  Int_t nstrC=19;
  Int_t nSectors=18;
  Int_t nStrips=nstrA+2*nstrB+2*nstrC;

  for (Int_t iSect = 0; iSect < nSectors; iSect++) {
    for (Int_t istr = 1; istr <= nStrips; istr++) {
    sdy = rnds->Gaus(0.,sigmatr);
    sdz = rnds->Gaus(0.,sigmatr);
      strId++;
      if( (TString(gSystem->Getenv("PARTGEOM")) == TString("kTRUE")) && !sActive[iSect] ) continue;
      new((*array)[j++]) AliAlignObjParams(AliGeomManager::SymName(idTOF,strId),AliGeomManager::LayerToVolUID(idTOF,strId), sdx, sdy, sdz, sdpsi, sdtheta, sdphi, kFALSE);
    }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save in file
    const char* filename = "TOFfullMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects in file %s", filename);
    f.cd();
    f.WriteObject(array,"TOFAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Silvia Arcelli");
    md->SetComment("Full misalignment for TOF");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TOF/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}


