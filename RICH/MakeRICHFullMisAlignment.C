void MakeRICHFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for RICH
  //
  Float_t sigmaTrans=0.1; // 1mm
  Float_t sigmaRot=0.001*180/TMath::Pi(); // 1 mrad
  Float_t dX, dY, dX;          Float_t dPsi, dTheta, dPhi;   //displacements

  TClonesArray *pCA = new TClonesArray("AliAlignObjMatrix",10);
  
  TRandom *pRnd   = new TRandom(4357);

  AliAlignObjMatrix o;
 
  Int_t idRICH =  AliAlignObj::kRICH;
  for (Int_t iCh = 0; iCh < 7; iCh++) {
    dX     = (pRnd->Uniform()-0.5)*sigmaTrans;    dY     = (pRnd->Uniform()-0.5)*sigmaTrans;    dZ     = (pRnd->Uniform()-0.5)*sigmaTrans;
    dPsi   = (pRnd->Uniform()-0.5)*sigmaRot;    dTheta = (pRnd->Uniform()-0.5)*sigmaRot;    dPhi   = (pRnd->Uniform()-0.5)*sigmaRot;
    new((*pCA)[iCh]) AliAlignObjMatrix(AliAlignObj::SymName(idRICH,iCh),
       AliAlignObj::LayerToVolUID(idRICH,iCh),dX,dY,dZ,dPsi,dTheta,dPhi,kTRUE);
  }

//   pCA->Print();
  
  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("RICHfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(pCA,"RICHAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData *pMeta= new AliCDBMetaData();  
    pMeta->SetResponsible("RICH Expert");
    pMeta->SetComment("Full alignment objects for RICH produced with sigmaTrans=1mm and sigmaRot=1mrad");
    pMeta->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("RICH/Align/Data",0,9999999);
    storage->Put(pCA,id,pMeta);
  }
  
  pCA->Delete();
}
