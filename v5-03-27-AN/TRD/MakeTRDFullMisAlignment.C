void MakeTRDFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TRD
  // Expects to read objects for FRAME
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjParams",1000);
  const char* macroname = "MakeTRDFullMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
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

  // Sigmas for the chambers
  Double_t smdx    = 0.3; // 3 mm
  Double_t smdy    = 0.3; // 3 mm
  Double_t smdz    = 0.3; // 3 mm
  Double_t smrx    = 0.4 / 1000.0 / TMath::Pi()*180; // 0.4 mrad
  Double_t smry    = 2.0 / 1000.0 / TMath::Pi()*180; // 2.0 mrad
  Double_t smrz    = 0.4 / 1000.0 / TMath::Pi()*180; // 0.4 mrad
  // Truncation for the chambers
  Double_t cutSmdx = 3.0 * smdx;
  Double_t cutSmdy = 3.0 * smdy;
  Double_t cutSmdz = 3.0 * smdz;

  // Sigmas for the chambers
  Double_t chdx    = 0.05;  // 0.5 mm
  Double_t chdy    = 0.1;   // 1.0 mm
  Double_t chdz    = 0.007; // 70  microns
  Double_t chrx    = 0.0005 / 1000.0 / TMath::Pi()*180; // 0 mrad
  Double_t chry    = 0.0005 / 1000.0 / TMath::Pi()*180; // 0 mrad
  Double_t chrz    = 0.3    / 1000.0 / TMath::Pi()*180; // 0.3 mrad
  // Truncation for the chambers
  Double_t cutChdx = 1.0  * chdx;
  Double_t cutChdy = 1.0  * chdy;
  Double_t cutChdz = 0.14 * chdz;

  Int_t sActive[18]={1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1};
  Double_t dx,dy,dz,rx,ry,rz;

  Int_t j=0;
  UShort_t volid;
  const char *symname;

  // create the supermodules' alignment objects
  for (int iSect; iSect<18; iSect++) {
    TString sm_symname(Form("TRD/sm%02d",iSect));
    dx = AliMathBase::TruncatedGaus(0.0,smdx,cutSmdx); 
    dy = AliMathBase::TruncatedGaus(0.0,smdy,cutSmdy); 
    dz = AliMathBase::TruncatedGaus(0.0,smdz,cutSmdz); 
    rx = gRandom->Rndm() * 2.0*smrx - smrx;
    ry = gRandom->Rndm() * 2.0*smry - smry;
    rz = gRandom->Rndm() * 2.0*smrz - smrz;
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
        dx = AliMathBase::TruncatedGaus(0.0,chdx,cutChdx); 
        dy = AliMathBase::TruncatedGaus(0.0,chdy,cutChdy); 
        dz = AliMathBase::TruncatedGaus(0.0,chdz,cutChdz); 
        rx = gRandom->Rndm() * 2.0*chrx - chrx;
        ry = gRandom->Rndm() * 2.0*chry - chry;
        rz = gRandom->Rndm() * 2.0*chrz - chrz;
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

