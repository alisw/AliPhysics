void MakeTRDFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TRD
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",1000);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

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

  Double_t dx,dy,dz,rx,ry,rz;

  Int_t j=0;
  TRandom *ran = new TRandom(4357);
  UShort_t volid;
  const char *symname;

  TGeoManager::Import("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO.root"); // geometry where the BSEGMO volumes have been misaligned

  // create the supermodules' alignment objects
  for (int i; i<18; i++) {
    TString sm_symname(Form("TRD/sm%02d",i));
    ran->Rannor(dx,rx);
    ran->Rannor(dy,ry);
    ran->Rannor(dz,rz);
    dx*=smdx;
    dy*=smdy;
    dz*=smdz;
    rx*=smrx;
    ry*=smry;
    rz*=smrz;
    new(alobj[j++]) AliAlignObjAngles(sm_symname.Data(),0,dx,dy,dz,rx,ry,rz,kFALSE);
  }

  for(Int_t k=0; k<18; k++){
    AliAlignObjAngles* smobj = (AliAlignObjAngles*)array->UncheckedAt(k);
    if(!smobj->ApplyToGeometry()){
      cout<<"application of object "<<k<<" failed!"<<endl;
      return;
    }
  }
  gGeoManager->Export("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO_trdSM.root");
  gGeoManager=0x0;
  TGeoManager::Import("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO_trdSM.root");
  // create the chambers' alignment objects
  ran = new TRandom(4357);
  for (Int_t iLayer = AliAlignObj::kTRD1; iLayer <= AliAlignObj::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      ran->Rannor(dx,rx);
      ran->Rannor(dy,ry);
      ran->Rannor(dz,rz);
      dx*=chdx;
      dy*=chdy;
      dz*=chdz;
      rx*=chrx;
      ry*=chry;
      rz*=chrz;
      volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      symname = AliAlignObj::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(symname,volid,dx,dy,dz,rx,ry,rz,kFALSE);
    }
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TRDfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TRDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Full misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


