void Geom()
{
  if(gClassTable->GetID("TGeoManager")<0) gSystem->Load("libGeom");
  
  Int_t copy;
  Double_t cm=1,m=100,mm=0.1;
  Double_t dx,dy,dz,rmin,rmin1,rmin2,rmax,rmax1,rmax2;
  Double_t a,z,den;
  
  new TGeoManager("ALICE-RICH-GM", "ALICE-RICH geometry");
  
  TGeoMedium *pAir      =new TGeoMedium("Air"   ,  1,gGeoManager->GetMaterial("Air"));
  TGeoMedium *pNeoceram =new TGeoMedium("Neoceram",1,new TGeoMaterial("Neoceram",a=0,z=0,den=0));
  TGeoMedium *pSiO2     =new TGeoMedium("SiO2",    1,new TGeoMaterial("SiO2",    a=0,z=0,den=0));
  TGeoMedium *pC6F14    =new TGeoMedium("C6F14",   1,new TGeoMaterial("C6F14",   a=0,z=0,den=0));
  TGeoMedium *pCH4      =new TGeoMedium("CH4",     1,new TGeoMaterial("CH4",     a=0,z=0,den=0));

  TGeoVolume *pMother=gGeoManager->MakeBox("ALICE",pAir,dx=10*m/2,dy=10*m/2,dz=30*m/2);
  
  TGeoVolume *pRICH     =gGeoManager->MakeBox( "RICH" ,pCH4,   dx=2*m/2,dy=2*m/2,dz=10*cm/2);      
  TGeoVolume *pRad      =gGeoManager->MakeBox( "RICH radiator" ,pC6F14,   dx=413*mm/2,dy=1330*mm/2,dz=24*mm/2);  
  TGeoVolume *pRadBack  =gGeoManager->MakeBox( "RadBack"  ,pNeoceram,dx=413*mm/2,dy=1330*mm/2,dz= 4*mm/2);pRadBack   ->SetLineColor(kRed);
  TGeoVolume *pRadLong  =gGeoManager->MakeBox( "RadLong"  ,pNeoceram,dx=  5*mm/2,dy=1330*mm/2,dz=15*mm/2);pRadLong   ->SetLineColor(kGreen);
  TGeoVolume *pRadShort =gGeoManager->MakeBox( "RadShort" ,pNeoceram,dx=403*mm/2,dy=   5*mm/2,dz=15*mm/2);pRadShort  ->SetLineColor(kGreen);
  TGeoVolume *pRadWin   =gGeoManager->MakeBox( "RadWin"   ,pSiO2,    dx=413*mm/2,dy=1330*mm/2,dz= 5*mm/2);pRadWin    ->SetLineColor(kBlue);
  TGeoVolume *pRadSpacer=gGeoManager->MakeTube("RadSpacer",pSiO2,   rmin=0,rmax=10*mm/2,dz= 15*mm/2);    pRadSpacer ->SetLineColor(kYellow);
//Nodes  
  pRad->AddNode(pRadBack,copy=1,new TGeoTranslation(0,0,-10.0*mm));
  pRad->AddNode(pRadWin ,copy=1,new TGeoTranslation(0,0,  9.5*mm));
  pRad->AddNode(pRadLong,copy=1,new TGeoTranslation(-403*mm/2-2.5*mm,0,0));
  pRad->AddNode(pRadLong,copy=2,new TGeoTranslation( 403*mm/2+2.5*mm,0,0));
  pRad->AddNode(pRadShort,copy=1,new TGeoTranslation(0,-1330*mm/2+5*mm/2     ,0));
  pRad->AddNode(pRadShort,copy=2,new TGeoTranslation(0,-1330*mm/2+5*mm+5*mm/2,0));
  pRad->AddNode(pRadShort,copy=3,new TGeoTranslation(0, 1330*mm/2-5*mm-5*mm/2,0));
  pRad->AddNode(pRadShort,copy=4,new TGeoTranslation(0, 1330*mm/2-5*mm/2     ,0));
  
  for(int i=-1;i<=1;i++)
    for(int j=0;j<10;j++)
      pRad->AddNode(pRadSpacer,copy=i+j+1,new TGeoTranslation(i*105*mm,-1330*mm/2+116*mm+j*122*mm,0));
  
  AliRICHParam *pPar=new AliRICHParam;
  
  for(copy=1;copy<=7;copy++)  
    pMother->AddNode(pRICH,copy,new TGeoTranslation(pPar->C(copy)->X(),pPar->C(copy)->Y(),pPar->C(copy)->Z()));
  
  pRICH->AddNode(pRad,copy=1,new TGeoTranslation(-(413*mm+25*mm),0,-12*mm));//out sourface of quartz window at z=0
  pRICH->AddNode(pRad,copy=2,new TGeoTranslation(  0            ,0,-12*mm));
  pRICH->AddNode(pRad,copy=3,new TGeoTranslation(  413*mm+25*mm ,0,-12*mm));
  
//Closing and drawing  
  gGeoManager->SetTopVolume(pMother); 
  gGeoManager->CloseGeometry();     
  gGeoManager->SetVisOption(0);   
  gGeoManager->GetMasterVolume()->Draw();  
//axises  
  Double_t X[6]={0,0,0,300,0,0};  Double_t Y[6]={0,0,0,0,300,0};  Double_t Z[6]={0,0,0,0,0,300};  
  TPolyLine3D *pXaxis=new TPolyLine3D(2,X);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,Y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,Z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();  
  
  gPad->GetView()->RotateView(270,30);
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

