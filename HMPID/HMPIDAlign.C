void RichAlign(Float_t sigmaTrans=0.1, Float_t sigmaRot=0.001)
{
  Float_t dX, dY, dZ;          Float_t dPsi, dTheta, dPhi;   //displacements

  TClonesArray *pCA = new TClonesArray("AliAlignObjMatrix",10000);
  
  TRandom *pRnd   = new TRandom(4357);

  AliAlignObjMatrix o;
 
  Int_t idHMPID =  AliAlignObj::kHMPID;
  for (Int_t iCh = 0; iCh < 7; iCh++) {
    dX     = pRnd->Gaus(0,sigmaTrans);    dY     = pRnd->Gaus(0,sigmaTrans);    dZ     = pRnd->Gaus(0,sigmaTrans);
    dPsi   = pRnd->Gaus(0,sigmaRot);      dTheta = pRnd->Gaus(0,sigmaRot);      dPhi   = pRnd->Gaus(0,sigmaRot);
    new((*pCA)[iCh]) AliAlignObjMatrix(AliAlignObj::GetVolPath(idHMPID,iCh), AliAlignObj::LayerToVolUID(idHMPID,iCh),dX,dY,dZ,dPsi,dTheta,dPhi);
  }

  pCA->Print();
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *pMeta= new AliCDBMetaData();  
  pMeta->SetResponsible("HMPID Expert");
  pMeta->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  AliCDBId id("HMPID/Align/Data",0,0); //you have to specify the run validity, although in the case of saving ideal objects makes not much sense
  AliCDBManager::Instance()->Put(pCA,id,pMeta);
  pCA->Delete();
}
