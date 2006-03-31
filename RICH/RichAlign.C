void RichAlign(Float_t sigmaTrans=0, Float_t sigmaRot=0)
{
  Float_t dX, dY, dX;          Float_t dPsi, dTheta, dPhi;   //displacements

  TClonesArray *pCA = new TClonesArray("AliAlignObjMatrix",10000);
  
  TRandom *pRnd   = new TRandom(4357);

  AliAlignObjMatrix o;
 
  Int_t idRICH =  AliAlignObj::kRICH;
  for (Int_t iCh = 0; iCh < 7; iCh++) {
    dX     = (pRnd->Uniform()-0.5)*sigmaTrans;    dY     = (pRnd->Uniform()-0.5)*sigmaTrans;    dZ     = (pRnd->Uniform()-0.5)*sigmaTrans;
    dPsi   = (pRnd->Uniform()-0.5)*sigmaRot;    dTheta = (pRnd->Uniform()-0.5)*sigmaRot;    dPhi   = (pRnd->Uniform()-0.5)*sigmaRot;
    new((*pCA)[iCh]) AliAlignObjMatrix(AliAlignObj::GetVolPath(idRICH,iCh), AliAlignObj::LayerToVolUID(idRICH,iCh),dX,dY,dZ,dPsi,dTheta,dPhi);
  }

  pCA->Print();
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  
  AliCDBMetaData *pMeta= new AliCDBMetaData();  
  pMeta->SetResponsible("RICH Expert");
  pMeta->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  AliCDBId id("RICH/Align/Data",0,0); //you have to specify the run validity, although in the case of saving ideal objects makes not much sense
  AliCDBManager::Instance()->Put(pCA,id,pMeta);
  pCA->Delete();
}
