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
  
  if( gSystem->Getenv("TOCDB") == TString("kTRUE") ){
    TString Storage = gSystem->Getenv("STORAGE");
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
  Info(macroname,"Loading FRAME alignment objects from CDB storage %s",
      Storage.Data());
  AliCDBPath fpath("GRP","Align","Data");
  if( gSystem->Getenv("TOCDB") == TString("kTRUE") ){
    AliCDBEntry *eFrame = storage->Get(fpath.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    TClonesArray* arFrame = (TClonesArray*) eFrame->GetObject();
    arFrame->Sort();
    Int_t nvols = arFrame->GetEntriesFast();
    Bool_t flag = kTRUE;
    for(Int_t j=0; j<nvols; j++)
      {
	AliAlignObj* alobj = (AliAlignObj*) arFrame->UncheckedAt(j);
	if (alobj->ApplyToGeometry() == kFALSE) flag = kFALSE;
      }
    if(!flag)
      Fatal(macroname,"Error in the application of FRAME alignment objects");
  }else{
    AliCDBEntry *eFrame = cdb->Get(fpath.GetPath());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    TClonesArray* arFrame = (TClonesArray*) eFrame->GetObject();
    arFrame->Sort();
    Int_t nvols = arFrame->GetEntriesFast();
    Bool_t flag = kTRUE;
    for(Int_t j=0; j<nvols; j++)
      {
	AliAlignObj* alobj = (AliAlignObj*) arFrame->UncheckedAt(j);
	if (alobj->ApplyToGeometry() == kFALSE) flag = kFALSE;
      }
    if(!flag)
      Fatal(macroname,"Error in the application of FRAME alignment objects");
  }

  //Produce objects for TOF supermodules
  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy vol id 

  Int_t nSMTOF = 18;
  Int_t j=0;
  Double_t smdx, smdy, smdz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(2345);
  Double_t sigmatr = 0.4; // max shift in cm w.r.t. local ideal RS
  Double_t sigmarot = 0.06; // max rot in deg w.r.t. local ideal RS (~ 1 mrad)
  
  for(Int_t i=0; i<nSMTOF; i++) {
    TString symname(Form("TOF/sm%02d",i));
    smdx = rnd->Gaus(0.,sigmatr);
    smdy = rnd->Gaus(0.,sigmatr);
    smdz =0;
    dpsi = 0.;
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = 0.;
    new((*array)[j++]) AliAlignObjParams(symname.Data(), dvoluid, smdx, smdy, smdz, dpsi, dtheta, dphi, kFALSE);
  }
  // Apply objects for TOF supermodules 
  for(Int_t i=0; i<nSMTOF; i++){
    AliAlignObjParams* smobj = (AliAlignObjParams*)array->UncheckedAt(i);
    if(!smobj->ApplyToGeometry()){
      cout<<"application of object "<<i<<" failed!"<<endl;
      return;
    }
  }

  //Produce objects for TOF strips (same sigmas as for residual misalignment)
  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;

  Double_t sdx=0.; 
  Double_t sdy=0.; 
  Double_t sdz=0.;
  Double_t sdpsi, sdtheta, sdphi;
  TRandom *rnds   = new TRandom(4357);
  sigmatr = 0.1; // max shift in cm w.r.t. local ideal RS

  for(i=0; i<AliGeomManager::LayerSize(idTOF); i++) {
    sdx = 0;
    sdy = rnds->Gaus(0.,sigmatr);
    sdz = rnds->Gaus(0.,sigmatr);
    sdpsi = 0.;
    sdtheta = 0.;
    sdphi = 0.;
    new((*array)[j++]) AliAlignObjParams(AliGeomManager::SymName(idTOF,i), AliGeomManager::LayerToVolUID(idTOF,i), sdx, sdy, sdz, sdpsi, sdtheta, sdphi, kFALSE);
  }

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
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


