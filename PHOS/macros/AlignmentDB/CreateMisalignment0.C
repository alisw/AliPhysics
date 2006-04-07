void CreateMisalignment0(const Int_t nModules=5)
{
  // *************************    1st step    ***************
  // Create TClonesArray of alignment objects for PHOS
  // with ideal geometry, i.e. zero displacement and disorientations
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",nModules);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  // null shifts and rotations

  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 

  TString basePath = "/ALIC_1/PHOS_"; 
  
  for (Int_t iModule = 1; iModule<=nModules; iModule++) {
    printf("Alignment object for %s is created\n",(basePath+iModule).Data());
    new(alobj[iModule-1]) AliAlignObjAngles((basePath+iModule).Data(),
					  dvoluid, dx, dy, dz, dpsi, dtheta, dphi);
  }

  // *************************    2nd step    ***************
  // Make CDB storage and put TClonesArray in
  // 
  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  AliCDBId id("PHOS/Align/Data",0,0);
  CDB->Put(array,id, md);
}





