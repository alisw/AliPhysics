void MakePMDZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for PMD
  //
  // Macro to randomly displace the 4 sectors of the PMD
  // in each plane. Each sector (to be misaligned) 
  // of PMD houses the following :
  // (a) 6 modules of preshower plane
  // (b) 6 modules of veto plane
  // (c) The FEE boards on back plane of each module
  // (d) 6 modules of convertor plates
  // The clustering is done module-wise
  // The actual amount displacement will be provided
  // by the survey data and has to be converted into
  // displacement in x,y,z,theta, phi and psi 
  
  
  // Now specify the path of the module to be misaligned
  // as followed in the PMD geant
  
  /*
     _____________
    |    |        |
    | 1  |   3    |
    |    |________|
    |____|___|    |
    |        | 2  |
    |   4    |    |
    |________|____|
    
    // Misalignment Matrix is expected to be
    // same for sectors 1 and 4 
    // and for the sectors 2 and 3
    // As these will be mounted on the same
    // Steel plate 
 */
  
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  //Create a TClonesArray of Align Object to store displacement Angles
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;
  
  AliAlignObjAngles o;
  
  Int_t iIndex=0; //  let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  Int_t i, j=0;

  for(i=1; i<=4; i++){
    TString snSector(Form("PMD/Sector%d",i));
    new(alobj[j++]) AliAlignObjAngles(snSector.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  if(!gSystem->Getenv("$TOCDB")){
    // Create a File to store the alignment data
    TFile f("PMDzeroMisalignment.root","RECREATE");
    if(!f) {cerr<<"cannot open file for output\n";}
    
    f.cd();
    f.WriteObject(array,"PMDAlignObjs ","kSingleKey");
    f.Close();
  }else{
  // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("");
    md->SetComment("Zero misalignment for PMD");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("PMD/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }
  array->Delete();

}
