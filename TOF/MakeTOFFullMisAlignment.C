void MakeTOFFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TOF
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",2000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("/home/rgrosso/Prove/AliRoot/geometry.root");

  AliAlignObjAngles a;
  Double_t sfdpsi=0.,sfdtheta=0.,sfdphi=0.;
  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy vol id 

  // FRAME part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  const char* basepath ="ALIC_1/B077_1/BSEGMO";
  TString segmpath;

  // for dead weight !!!! carefull: in mm
  ifstream dw_file("/home/rgrosso/MacroAllineamento010207/disp_deadweight.txt", ios::in);
  if(!dw_file) {cerr<<"cannot open dead weight file for input\n"; return;}
  TArrayD dwdx(36);
  TArrayD dwdy(36);
  TArrayD dwdz(36);
  TArrayD dwdrot(36);
  // for additional load
  ifstream orig_file("/home/rgrosso/MacroAllineamento010207/disp_full.txt", ios::in);
//   ifstream orig_file("disp_full_nohmpid.txt", ios::in);
  if(!orig_file) {cerr<<"cannot open file for input\n"; return;}
  TArrayD adx(36);
  TArrayD ady(36);
  TArrayD adz(36);
  TArrayD adrot(36);
  // avarage displacements
  Double_t mean_dx[18];
  Double_t mean_dy[18];
  Double_t idx[18];
  Double_t idy[18];
  Double_t newx[18];
  Double_t newy[18];

  string buffer;
  Int_t point;
  Double_t dx,dy,dz,drot;

  // fill in from file the array for dead weight deformations
  dw_file>>point>>dx>>dy>>dz>>drot;
  dwdx.AddAt(dx,point-1);dwdy.AddAt(dy,point-1);dwdz.AddAt(dz,point-1);dwdrot.AddAt(drot,point-1);
  while(getline(dw_file,buffer))
    {
      dw_file>>point>>dx>>dy>>dz>>drot;
      dwdx.AddAt(dx,point-1);dwdy.AddAt(dy,point-1);dwdz.AddAt(dz,point-1);dwdrot.AddAt(drot,point-1);
    }

  // fill in from file the array for remaining load deformations
  orig_file>>point>>dx>>dy>>dz>>drot;
  adx.AddAt(dx,point-1);ady.AddAt(dy,point-1);adz.AddAt(dz,point-1);adrot.AddAt(drot,point-1);
  while(getline(orig_file,buffer))
    {
      orig_file>>point>>dx>>dy>>dz>>drot;
      adx.AddAt(dx,point-1);ady.AddAt(dy,point-1);adz.AddAt(dz,point-1);adrot.AddAt(drot,point-1);
    }

  //calculate the displacement of the center of the sm neglecting rotations,
  //thus as average of the four surrounding point for dw + other load
  // Prepare also to plot with respect to ideal circle.
  Double_t cx=0.5;  //just for plotting
  Double_t cy=0.5;
  Double_t rr=0.3;
  Double_t kk = TMath::Pi()*20./180.;
  Double_t scale = 100./4000.;
  TPolyMarker* iddu = new TPolyMarker(18,newx,newy);
  TPolyMarker* disp = new TPolyMarker(18,newx,newy);

  Int_t sm,outc,outclw; //outer number counterclock and clockwise
  for(sm=0; sm<18; sm++){
    outc=5-sm;
    if(outc<1) outc+=18;
    outclw=outc-1;
    if(outclw<1) outclw+=18;
    mean_dx[sm]=0.125*(adx[outc-1]+adx[outclw-1]+adx[outc+18-1]+adx[outclw+18-1]+dwdx[outc-1]+dwdx[outclw-1]+dwdx[outc+18-1]+dwdx[outclw+18-1]);
    mean_dy[sm]=0.125*(ady[outc-1]+ady[outclw-1]+ady[outc+18-1]+ady[outclw+18-1]+dwdy[outc-1]+dwdy[outclw-1]+dwdy[outc+18-1]+dwdy[outclw+18-1]);
    idx[sm]=cx+rr*TMath::Sin(sm*kk);
    idy[sm]=cy+rr*TMath::Cos(sm*kk);
    newx[sm]=idx[sm]+scale*mean_dx[sm];
    newy[sm]=idy[sm]+scale*mean_dy[sm];
    iddu->SetPoint(sm,idx[sm],idy[sm]);
    disp->SetPoint(sm,newx[sm],newy[sm]);
    
    segmpath=basepath;
    segmpath+=sm;
    segmpath+="_1";
    cout<<segmpath.Data()<<"  "<<dvoluid<<"  "<<mean_dx[sm]*0.1<<"  "<<mean_dy[sm]*0.1<<"  "<<dz<<"  "<<sfdpsi<<"  "<<sfdtheta<<"  "<<sfdphi<<endl;
    new(alobj[sm]) AliAlignObjAngles(segmpath.Data(), dvoluid, mean_dx[sm]*0.1,
		     mean_dy[sm]*0.1, dz, sfdpsi, sfdtheta, sfdphi, kTRUE);
  }
  
  for(Int_t k=0; k<18; k++){
    AliAlignObjAngles* smobj = (AliAlignObjAngles*)array->UncheckedAt(k);
    if(!smobj->ApplyToGeometry()){
      cout<<"application of object "<<k<<" failed!"<<endl;
      return;
    }
  }
  gGeoManager->Export("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO.root");
  gGeoManager=0x0;
  TGeoManager::Import("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO.root");

  // TOF part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Int_t nSMTOF = 18;
  Int_t i;
  Int_t j=18;
  Double_t tofdx, tofdy, tofdz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(2345);
  Double_t sigmatr = 0.4; // max shift in cm w.r.t. local ideal RS
  Double_t sigmarot = 0.06; // max rot in deg w.r.t. local ideal RS (~ 1 mrad)
  
  for(i=0; i<18; i++) {
    TString symname(Form("TOF/sm%02d",i));
    tofdx = rnd->Gaus(0.,sigmatr);
    tofdy = 0;
    tofdz = rnd->Gaus(0.,sigmatr);
    dpsi = 0;
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = 0.;
    new(alobj[j++]) AliAlignObjAngles(symname.Data(), dvoluid, tofdx, tofdy, tofdz, dpsi, dtheta, dphi, kFALSE);
  }
  for(Int_t k=18; k<36; k++){
    AliAlignObjAngles* smobj = (AliAlignObjAngles*)array->UncheckedAt(k);
    if(!smobj->ApplyToGeometry()){
      cout<<"application of object "<<k<<" failed!"<<endl;
      return;
    }
  }
  gGeoManager->Export("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO_tofSM.root");


  gGeoManager=0x0;
  TGeoManager::Import("/home/rgrosso/MacroAllineamento010207/geom_misalBSEGMO_tofSM.root");
  // tof strips as in residual
  AliAlignObj::ELayerID idTOF = AliAlignObj::kTOF;

  Double_t sdx=0.; 
  Double_t sdy=0.; 
  Double_t sdz=0.;
  Double_t sdpsi, sdtheta, sdphi;
  TRandom *rnds   = new TRandom(4357);
  sigmatr = 0.1; // max shift in cm w.r.t. local ideal RS

  for(i=0; i<AliAlignObj::LayerSize(idTOF); i++) {
    sdx = 0;
    sdy = rnds->Gaus(0.,sigmatr);
    sdz = rnds->Gaus(0.,sigmatr);
    sdpsi = 0.;
    sdtheta = 0.;
    sdphi = 0.;
    new(alobj[j++]) AliAlignObjAngles(AliAlignObj::SymName(idTOF,i), AliAlignObj::LayerToVolUID(idTOF,i), sdx, sdy, sdz, sdpsi, sdtheta, sdphi, kFALSE);
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TOFfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TOFAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Silvia Arcelli");
    md->SetComment("Full misalignment for TOF and FRAME");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TOF/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


