void MakeTRDResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for TRD
  //
  const char* macroname = "MakeTRDResMisAlignment.C";
  TClonesArray *array = new TClonesArray("AliAlignObjParams",1000);
  TClonesArray &alobj = *array;
   
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

  // sigmas for the chambers
  Double_t chdx=0.002; // 20 microns
  Double_t chdy=0.003; // 30 microns
  Double_t chdz=0.007; // 70 microns
  Double_t chrx=0.3/1000/TMath::Pi()*180; // 0.3 mrad
  Double_t chry=0.3/1000/TMath::Pi()*180; // 0.3 mrad
  Double_t chrz=0.1/1000/TMath::Pi()*180; // 0.1 mrad

  Double_t dx=0.,dy=0.,dz=0.,rx=0.,ry=0.,rz=0.;

  Int_t j=0;
  TRandom *ran = new TRandom(4357);
  UShort_t volid;
  const char* symname; 

  // create the supermodules' alignment objects
  for (int i; i<18; i++) {
    TString sm_symname(Form("TRD/sm%02d",i));
    new((*array)[j++])
      AliAlignObjParams(sm_symname.Data(),0,dx,dy,dz,rx,ry,rz,kTRUE);
  }
 
  // create the chambers' alignment objects
  for (Int_t iLayer = AliGeomManager::kTRD1; iLayer <= AliGeomManager::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {
      ran->Rannor(dx,rx);
      ran->Rannor(dy,ry);
      ran->Rannor(dz,rz);
      dx*=chdx;
      dy*=chdy;
      dz*=chdz;
      rx*=chrx;
      ry*=chry;
      rz*=chrz;
      volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      symname = AliGeomManager::SymName(volid);
      new(alobj[j++]) AliAlignObjParams(symname,volid,dx,dy,dz,rx,ry,rz,kFALSE);
    }
  }

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "TRDresidualMisalignment.root";
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
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Residual misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}


