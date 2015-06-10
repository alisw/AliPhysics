// macro to create 
Double_t sgXMod=30e-4,sgYMod=30e-4,sgZMod=30e-4,sgThtMod=0.1,sgPsiMod=0.1,sgPhiMod=0.1;
Double_t sgXSta=30e-4,sgYSta=30e-4,sgZSta=30e-4,sgThtSta=0.1,sgPsiSta=0.1,sgPhiSta=0.1;
Double_t sgXLay=30e-4,sgYLay=30e-4,sgZLay=30e-4,sgThtLay=0.1,sgPsiLay=0.1,sgPhiLay=0.1;
Double_t sgXITS=100e-4,sgYITS=100e-4,sgZITS=2000e-4,sgThtITS=0.1,sgPsiITS=0.1,sgPhiITS=0.1;
//
TClonesArray* deltas=0;
//
void MakeITSUResMisAlignment()       			 
{
  gSystem->Load("libITSUpgradeBase");
  //
  // at the moment we don't want to write to official CDB
  gSystem->Setenv("TOCDB","kTRUE");
  gSystem->Setenv("STORAGE","local://");
  //
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://");
  cdb->SetRun(0);
  AliCDBStorage* storage;
  // 
  if (!gGeoManager) {
    if (!gSystem->AccessPathName("geometry.root")) {
      printf("Loading geometry.root from current directory\n");
      AliGeomManager::LoadGeometry("geometry.root"); //load geom from default CDB storage      
    }
  }
  //
  if (TString(gSystem->Getenv("TOCDB")) == TString("kTRUE")) {
    TString Storage = gSystem->Getenv("STORAGE");
    if (!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if (!storage) {
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    //
    if (!gGeoManager) {
      AliCDBPath path("GRP","Geometry","Data");
      AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
      if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
      entry->SetOwner(0);
      TGeoManager* geom = (TGeoManager*) entry->GetObject();
      AliGeomManager::SetGeometry(geom);
    }
  }
  //
  const UShort_t dummyVID = 0xffff;
  AliITSUGeomTGeo* gm0 = new AliITSUGeomTGeo(kTRUE);
  //
  deltas = new TClonesArray("AliAlignObjParams");
  double dx,dy,dz,dtht,dpsi,dphi;
  //
  TString sname;
  int idel = 0;
  //
  dx   = sgXITS*gRandom->Gaus();
  dy   = sgYITS*gRandom->Gaus();	
  dz   = sgZITS*gRandom->Gaus();
  dtht = sgThtITS*gRandom->Gaus(); 
  dpsi = sgPsiITS*gRandom->Gaus(); 
  dphi = sgPhiITS*gRandom->Gaus(); 	
  sname = gm0->ComposeSymNameITS();
  new( (*deltas)[idel++] ) AliAlignObjParams(sname.Data(),dummyVID,
					     dx,dy,dz,dtht,dpsi,dphi,kTRUE);
  //
  for (int ilr=0;ilr<gm0->GetNLayers();ilr++) {
    //
    dx   = sgXLay*gRandom->Gaus();
    dy   = sgYLay*gRandom->Gaus();	
    dz   = sgZLay*gRandom->Gaus();
    dtht = sgThtLay*gRandom->Gaus(); 
    dpsi = sgPsiLay*gRandom->Gaus(); 
    dphi = sgPhiLay*gRandom->Gaus(); 	
    sname = gm0->ComposeSymNameLayer(ilr);
    new( (*deltas)[idel++] ) AliAlignObjParams(sname.Data(),dummyVID,
					       dx,dy,dz,dtht,dpsi,dphi,kTRUE);
    //
    for (int ild=0;ild<gm0->GetNStaves(ilr);ild++) {
      //
      dx   = sgXSta*gRandom->Gaus();
      dy   = sgYSta*gRandom->Gaus();	
      dz   = sgZSta*gRandom->Gaus();
      dtht = sgThtSta*gRandom->Gaus(); 
      dpsi = sgPsiSta*gRandom->Gaus(); 
      dphi = sgPhiSta*gRandom->Gaus(); 	
      sname = gm0->ComposeSymNameStave(ilr,ild);
      new( (*deltas)[idel++] ) AliAlignObjParams(sname.Data(),dummyVID,
						 dx,dy,dz,dtht,dpsi,dphi,kTRUE);
      //
      for (int isn=0;isn<gm0->GetNChipsPerStave(ilr);isn++) {
	dx   = sgXMod*gRandom->Gaus();
	dy   = sgYMod*gRandom->Gaus();	
	dz   = sgZMod*gRandom->Gaus();
	dtht = sgThtMod*gRandom->Gaus(); 
	dpsi = sgPsiMod*gRandom->Gaus(); 
	dphi = sgPhiMod*gRandom->Gaus(); 	
	int mid = gm0->GetChipIndex(ilr,ild,isn);
	sname = gm0->ComposeSymNameChip(ilr,ild,-1,-1,isn);
	new( (*deltas)[idel++] ) AliAlignObjParams(sname.Data(),gm0->ChipVolUID(mid),
						   dx,dy,dz,dtht,dpsi,dphi,kTRUE);
      }
    }
  }
  //
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ITSUresidualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(deltas,"ITSUAlignObjs","kSingleKey");
    f.Close();
  } 
  else{
    // save in CDB storage
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("R.S.");
    md->SetComment("Residual misalignment for ITSU");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(deltas,id,md);
  }

  // apply to geometry
  printf("Applying created misalignment to geometry in memory\n");
  AliGeomManager::ApplyAlignObjsToGeom(*deltas);
  gGeoManager->LockGeometry();
  //
}
