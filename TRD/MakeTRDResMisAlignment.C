void MakeTRDResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for TRD
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",1000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  // sigmas for the supermodules
  Double_t smdx=3; // 300 microns
  Double_t smdy=3; // 300 microns
  Double_t smdz=3; // 300 microns
  Double_t smrx=0.4/1000/TMath::Pi()*180; // 0.4 mrad
  Double_t smry=2.0/1000/TMath::Pi()*180; // 2 mrad
  Double_t smrz=0.4/1000/TMath::Pi()*180; // 0.4 mrad

  // sigmas for the chambers
  Double_t chdx=0.02; // 200 microns
  Double_t chdy=0.03; // 300 microns
  Double_t chdz=0.07; // 700 microns
  Double_t chrx=0.3/1000/TMath::Pi()*180; // 0.3 mrad
  Double_t chry=0.3/1000/TMath::Pi()*180; // 0.3 mrad
  Double_t chrz=0.1/1000/TMath::Pi()*180; // 0.1 mrad

  Double_t dx,dy,dz,rx,ry,rz;

  Int_t j=0;
  TRandom *ran = new TRandom(4357);
  UShort_t volid;
  const char *path;

  //generate the paths for the 18 supermodules
  Double_t r=325.;
  Double_t phi;
  Int_t isec;
  TString ts0;
  TString* sm_path = new TString[18];
  TRegexp regexp(".*BTRD[0-9][^/]*");

  for (isec=0; isec<18; isec++) {
    phi=(isec+0.5)/18*2*TMath::Pi();
    gGeoManager->FindNode(r*cos(phi),r*sin(phi),0.);
    ts0 = gGeoManager->GetPath();
    sm_path[isec] = ts0(regexp);
    sm_path[isec] +='\0';
  }

  // create the supermodules' alignment objects
  for (int i; i<18; i++) {
    ran.Rannor(dx,rx);
    ran.Rannor(dy,ry);
    ran.Rannor(dz,rz);
    dx*=smdx;
    dy*=smdy;
    dz*=smdz;
    rx*=smrx;
    ry*=smry;
    rz*=smrz;
    new(alobj[j++]) AliAlignObjAngles(sm_path[i].Data(),0,dx,dy,dz,rx,ry,rz,kFALSE);
  }

  // create the chambers' alignment objects
  for (Int_t iLayer = AliAlignObj::kTRD1; iLayer <= AliAlignObj::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      ran.Rannor(dx,rx);
      ran.Rannor(dy,ry);
      ran.Rannor(dz,rz);
      dx*=chdx;
      dy*=chdy;
      dz*=chdz;
      rx*=chrx;
      ry*=chry;
      rz*=chrz;
      volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      path = AliAlignObj::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(path,volid,dx,dy,dz,rx,ry,rz,kFALSE);
    }
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TRDresidualMisalignment.root","RECREATE");
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
    md->SetComment("Residual misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


