void Geom()
{
  if(gClassTable->GetID("TGeoManager")<0) gSystem->Load("libGeom");
  
  Int_t copy;
  Double_t cm=1,m=100,mm=0.1;
  Double_t dx,dy,dz,rmin,rmin1,rmin2,rmax,rmax1,rmax2;
  Double_t a,z,den;
  
  new TGeoManager("ALICE-RICH-GM", "ALICE-RICH geometry");
  
  TGeoMedium *pVacuum  =new TGeoMedium("Vacuum",  1,new TGeoMaterial("Vacuum",  a=0,z=0,den=0));
  TGeoMedium *pNeoceram=new TGeoMedium("Neoceram",1,new TGeoMaterial("Neoceram",a=0,z=0,den=0));
  TGeoMedium *pSiO2    =new TGeoMedium("SiO2",    1,new TGeoMaterial("SiO2",    a=0,z=0,den=0));
  TGeoMedium *pC6F14   =new TGeoMedium("C6F14",   1,new TGeoMaterial("C6F14",   a=0,z=0,den=0));

  TGeoVolume *pMother=gGeoManager->MakeBox("MOTHER",pVacuum,dx=1000,dy=1000,dz=3000);
  
  
  TGeoVolume *pRadiator =gGeoManager->MakeBox("Radiator" ,pC6F14,   dx=413*mm/2,dy=1330*mm/2,dz=24*mm/2);  
  TGeoVolume *pRadBack  =gGeoManager->MakeBox("RadBack"  ,pNeoceram,dx=413*mm/2,dy=1330*mm/2,dz= 4*mm/2);pRadBack->SetLineColor(kRed);
  TGeoVolume *pRadLong  =gGeoManager->MakeBox("RadLong"  ,pNeoceram,dx=  5*mm/2,dy=1330*mm/2,dz=15*mm/2);pRadLong->SetLineColor(kGreen);
  TGeoVolume *pRadShort =gGeoManager->MakeBox("RadShort" ,pNeoceram,dx=403*mm/2,dy=   5*mm/2,dz=15*mm/2);pRadShort->SetLineColor(kGreen);
  TGeoVolume *pRadWin   =gGeoManager->MakeBox("RadWin"   ,pSiO2,    dx=413*mm/2,dy=1330*mm/2,dz= 5*mm/2);pRadWin ->SetLineColor(kBlue);
  TGeoVolume *pRadSpacer=gGeoManager->MakeTube("RadSpacer",pSiO2,   rmin=0,rmax=10*mm/2,dz= 15*mm/2);pRadSpacer ->SetLineColor(kYellow);
  pRadiator->AddNode(pRadBack,copy=1,new TGeoTranslation(0,0,-10.0*mm));
  pRadiator->AddNode(pRadWin ,copy=1,new TGeoTranslation(0,0,  9.5*mm));
  pRadiator->AddNode(pRadLong,copy=1,new TGeoTranslation(-403*mm/2-2.5*mm,0,0));
  pRadiator->AddNode(pRadLong,copy=2,new TGeoTranslation( 403*mm/2+2.5*mm,0,0));
  pRadiator->AddNode(pRadShort,copy=1,new TGeoTranslation(0,-1330*mm/2+5*mm/2     ,0));
  pRadiator->AddNode(pRadShort,copy=2,new TGeoTranslation(0,-1330*mm/2+5*mm+5*mm/2,0));
  pRadiator->AddNode(pRadShort,copy=3,new TGeoTranslation(0, 1330*mm/2-5*mm-5*mm/2,0));
  pRadiator->AddNode(pRadShort,copy=4,new TGeoTranslation(0, 1330*mm/2-5*mm/2     ,0));
  
  for(int i=-1;i<=1;i++)
    for(int j=0;j<10;j++)
      pRadiator->AddNode(pRadSpacer,copy=j+1,new TGeoTranslation(i*105*mm,-1330*mm/2+116*mm+j*122*mm,0));
  
  
  pMother->AddNode(pRadiator,copy=1,new TGeoTranslation(-60*cm,0,-12*mm));//out sourface of quartz window at z=0
  pMother->AddNode(pRadiator,copy=2,new TGeoTranslation(  0   ,0,-12*mm));
  pMother->AddNode(pRadiator,copy=3,new TGeoTranslation( 60*cm,0,-12*mm));
  gGeoManager->SetTopVolume(pMother);  
  gGeoManager->CloseGeometry();   
  pMother->Draw();
  new TBrowser;
}//Geom()



//______________________________________________________________________________
void Dist()
{
  if(gClassTable->GetID("TGeoManager")<0) gSystem->Load("libGeom");
  new TGeoManager("GM", "kir test"); // you must have a manager
  TGeoShape *shape = new TGeoSphere(10, 20); // Rmin, Rmax
  Double_t point[3]={30,0,0};
  Double_t dir[3]={1,0,0};
  Bool_t isInside=shape->Contains(point);
  Double_t dist;
  if(isInside)
    dist = shape->DistToOut(point,dir,3);
  else
    dist = shape->DistToIn(point,dir,3);
  cout<< "dist="<<dist<<endl;
  delete gGeoManager;
}

