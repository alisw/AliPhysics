void MakeHMPIDZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for HMPID
  //
  TClonesArray *pCA = new TClonesArray("AliAlignObjMatrix",10);
  
  AliAlignObjMatrix o;

  Double_t dX=0.,dY=0.,dZ=0.,dPsi=0.,dTheta=0.,dPhi=0.;
 
  Int_t idHMPID =  AliAlignObj::kHMPID;
  for (Int_t iCh = 0; iCh < 7; iCh++) {
    new((*pCA)[iCh]) AliAlignObjMatrix(AliAlignObj::SymName(idHMPID,iCh),AliAlignObj::LayerToVolUID(idHMPID,iCh),dX,dY,dZ,dPsi,dTheta,dPhi,kTRUE);
  }

//   pCA->Print();
  
  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("HMPIDzeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(pCA,"HMPIDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData *pMeta= new AliCDBMetaData();  
    pMeta->SetResponsible("HMPID Expert");
    pMeta->SetComment("Zero alignment objects for HMPID");
    pMeta->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("HMPID/Align/Data",0,9999999);
    storage->Put(pCA,id,pMeta);
  }
  
  pCA->Delete();
}
