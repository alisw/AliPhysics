red()
{
  gSystem->Load("libMinuit.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libCDB.so");

  gSystem->Load("libGed.so");
  gSystem->Load("libRGL.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");

  gSystem->Load("libReve.so");

  gSystem->Load("libRICHbase.so");
  gSystem->Load("libRICHsim.so");
  gSystem->Load("libRICHrec.so");
  
  TGeoManager::Import("geometry.root");
  
  AliRICHParam *pParam=AliRICHParam::Instance();
  
  AliRunLoader *pAL=AliRunLoader::Open(); pAL->LoadHits  ("RICH"); TTree *pHitT=pAL->GetTreeH("RICH", false);
                                          pAL->LoadDigits("RICH"); TTree *pDigT=pAL->GetTreeD("RICH", false); 
                                          
  Reve::PointSet *pHitPnt=new Reve::PointSet("Hits"));
  TPolyMarker3D  *pDigPnt=new TPolyMarker3D;  pDigPnt->SetName("Digits"); pDigPnt->SetMarkerColor(kGreen);iPntCnt=0;
                                          
  TPointSelector ps(pHitT,pHitPnt,"fX:fY:fZ",""); ps.Select();
  
  TClonesArray *pDigLst=new TClonesArray("AliRICHDigit"); //this is tmp dig list per chamber
 
  for(Int_t iCh=0;iCh<7;iCh++){
    pDigT->SetBranchAddress(Form("RICH%i",iCh+1),&pDigLst);
    pDigT->GetEntry(0);
    for(Int_t iDig=0;iDig<pDigLst->GetEntries();iDig++){
      AliRICHDigit *pDig=(AliRICHDigit*)pDigLst->At(iDig);    
      TVector2 lors=pParam->Pad2Loc(pDig->PadX(),pDig->PadY());
      TVector3 mars=pParam->Lors2Mars(iCh,lors.X(),lors.Y());
      pDigPnt->SetPoint(iPntCnt++,mars.X(),mars.Y(),mars.Z());
    }//digits loop for chamber        
  }//chambers loop
  
  if(!gReve) new Reve::RGTopFrame(0,0,0,2);
  gReve->AddGlobalRenderElement(new Reve::RenderElementObjPtr(pDigPnt));
  gReve->AddRenderElement(pHitPnt);
  gReve->Redraw3D();
}
