void MakeHMPIDZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for HMPID
  //
  TClonesArray *pCA = new TClonesArray("AliAlignObjMatrix",10);
  
  Double_t dX=0.,dY=0.,dZ=0.,dPsi=0.,dTheta=0.,dPhi=0.;
 
  Int_t idHMPID =  AliGeomManager::kHMPID;
  for (Int_t iCh = 0; iCh < 7; iCh++) {
    new((*pCA)[iCh]) AliAlignObjMatrix(AliGeomManager::SymName(idHMPID,iCh),AliGeomManager::LayerToVolUID(idHMPID,iCh),dX,dY,dZ,dPsi,dTheta,dPhi,kTRUE);
  }

//   pCA->Print();
  const char* macroname = "MakeHMPIDZeroMisAlignment.C";  
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "HMPIDzeroMisalignment.root";
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
    pMeta->SetComment("Zero alignment objects for HMPID");
    pMeta->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("HMPID/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(pCA,id,pMeta);
  }
  
  pCA->Delete();
}
