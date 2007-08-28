void MakeHMPIDResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for HMPID
  //
  Float_t sigmaTrans=0.1; // 1mm
  Float_t sigmaRot=0.001*180/TMath::Pi(); // 1 mrad
  Float_t dX, dY, dX;          Float_t dPsi, dTheta, dPhi;   //displacements

  TClonesArray *pCA = new TClonesArray("AliAlignObjMatrix",10);
  
  TRandom *pRnd   = new TRandom(4357);

  Int_t idHMPID =  AliGeomManager::kHMPID;
  for (Int_t iCh = 0; iCh < 7; iCh++) {
    dX     = (pRnd->Uniform()-0.5)*sigmaTrans;    dY     = (pRnd->Uniform()-0.5)*sigmaTrans;    dZ     = (pRnd->Uniform()-0.5)*sigmaTrans;
    dPsi   = (pRnd->Uniform()-0.5)*sigmaRot;    dTheta = (pRnd->Uniform()-0.5)*sigmaRot;    dPhi   = (pRnd->Uniform()-0.5)*sigmaRot;
    new((*pCA)[iCh]) AliAlignObjMatrix(AliGeomManager::SymName(idHMPID,iCh),AliGeomManager::LayerToVolUID(idHMPID,iCh),dX,dY,dZ,dPsi,dTheta,dPhi,kTRUE);
  }

//   pCA->Print();
  
  const char* macroname = "MakeHMPIDResMisAlignment.C";
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "HMPIDresidualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(pCA,"HMPIDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData *pMeta= new AliCDBMetaData();  
    pMeta->SetResponsible("HMPID Expert");
    pMeta->SetComment("Residual alignment objects for HMPID produced with sigmaTrans=1mm and sigmaRot=1mrad");
    pMeta->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("HMPID/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(pCA,id,pMeta);
  }
  
  pCA->Delete();
}
