
TGeoManager *g=0;

void Hgeom(Bool_t isOnlyChambers=kFALSE)
{
  
  g=new TGeoManager("TestHMPID","Private HMPID geometry");
  Materials();  
  gGeoManager->MakeBox("ALIC",gGeoManager->GetMedium("HMPID_CH4"),3000/2,3000/2,3000/2); //arbitrary values  
  gGeoManager->SetTopVolume(gGeoManager->GetVolume("ALIC"));
    
  Hmpid(isOnlyChambers);
  
  gGeoManager->CloseGeometry();
  
//  gGeoManager->SetVisOption(0); gGeoManager->SetVisLevel(5); 
  
//  gGeoManager->GetMasterVolume()->Draw();
//  Axis();
//  gPad->GetView()->SetView(3,94,-70,0);
  new TBrowser;
}
//__________________________________________________________________________________________________
void Materials()
{
//Media for HMPID
  Double_t a=0,z=0,den=0,radlen=0,intlen=0;//tmp vars for material parameters
                                //  A       Z     rho
  TGeoMaterial *ar,*al,*cu,*w,*ro;  TGeoMixture *c6f14,*sio2,*ch4,*csi;
  al   =new TGeoMaterial("HMPID_Al"   ,  26.982 , 13 ,  2.700);  
  ro   =new TGeoMaterial("HMPID_Neoc" ,  12.01  ,  6 ,  0.1  ); 
  cu   =new TGeoMaterial("HMPID_Cu"   ,  63.546 , 29 ,  8.960); 
  w    =new TGeoMaterial("HMPID_W"    , 183.840 , 74 , 19.300); 
  ar   =new TGeoMaterial("HMPID_Ar"   ,  39.948 , 18 ,  1.396);  
  c6f14=new TGeoMixture ("HMPID_C6F14",2,1.68     );c6f14->DefineElement(0, 12.010, 6,6); c6f14->DefineElement(1, 18.998, 9,14);
  sio2 =new TGeoMixture ("HMPID_SiO2" ,2,2.20     );sio2 ->DefineElement(0, 28.085,14,1);  sio2->DefineElement(1, 15.999, 8, 2);
  ch4  =new TGeoMixture ("HMPID_CH4"  ,2,0.4224e-3);ch4  ->DefineElement(0, 12.010, 6,1);   ch4->DefineElement(1,  1.007, 1, 4);                
  csi  =new TGeoMixture ("HMPID_CsI"  ,2,1.8      );csi  ->DefineElement(0,132.900,55,1);   csi->DefineElement(1,126.900,53, 1);
  
  new TGeoMedium("HMPID_Al"   ,1,al);
  new TGeoMedium("HMPID_Neoc" ,2,ro);
  new TGeoMedium("HMPID_Cu"   ,3,cu);
  new TGeoMedium("HMPID_W"    ,4,w );
  new TGeoMedium("HMPID_C6F14",6,c6f14);
  new TGeoMedium("HMPID_SiO2" ,7,sio2); 
  new TGeoMedium("HMPID_CH4"  ,8,ch4);
  new TGeoMedium("HMPID_CsI"  ,9,csi);     
  new TGeoMedium("HMPID_Ar"   ,1,ar);
                                                          
}//Materials()
//__________________________________________________________________________________________________
void Hmpid(Bool_t isOnlyChambers)
{
  Double_t cm=1,mm=0.1*cm,um=0.0001*cm;//length units  
  TGeoMedium *al   =gGeoManager->GetMedium("HMPID_Al");    
  TGeoMedium *ch4  =gGeoManager->GetMedium("HMPID_CH4");    
  TGeoMedium *roha =gGeoManager->GetMedium("HMPID_Rohacell");    
  TGeoMedium *neoc =gGeoManager->GetMedium("HMPID_Neoceram");    
  TGeoMedium *c6f14=gGeoManager->GetMedium("HMPID_C6F14");    
  TGeoMedium *sio2 =gGeoManager->GetMedium("HMPID_SiO2");    
  TGeoMedium *cu   =gGeoManager->GetMedium("HMPID_Cu");    
  TGeoMedium *w    =gGeoManager->GetMedium("HMPID_W");    
  TGeoMedium *csi  =gGeoManager->GetMedium("HMPID_CsI");    
  TGeoMedium *ar   =gGeoManager->GetMedium("HMPID_Ar");    
  
  TGeoVolume *hmp=gGeoManager->MakeBox ("Hmp",ch4,1681*mm/2, 1466*mm/2,(2*80*mm+2*41*mm)/2);//2033P1  z from p84 TDR  
  const Double_t kAngHor=19.5; //  horizontal angle between chambers  19.5 grad
  const Double_t kAngVer=20;   //  vertical angle between chambers    20   grad     
  const Double_t kAngCom=30;   //  common HMPID rotation around z     30   grad     
  const Double_t trans[3]={490*cm,0*cm,0*cm}; //center of the chamber is on window- proximity gap surface
  for(Int_t ch=0;ch<=6;ch++){//place 7 chambers
    TGeoHMatrix *mat=new TGeoHMatrix;
    mat->RotateY(90);           //rotate around y since initial position is in XY plane -> now in YZ plane
    mat->SetTranslation(trans); //now plane in YZ is shifted along x 
    switch(ch){
      case 0:                mat->RotateY(kAngHor);  mat->RotateZ(-kAngVer);  break; //right and down 
      case 1:                                        mat->RotateZ(-kAngVer);  break; //down              
      case 2:                mat->RotateY(kAngHor);                           break; //right 
      case 3:                                                                 break; //no rotation
      case 4:                mat->RotateY(-kAngHor);                          break; //left   
      case 5:                                        mat->RotateZ(kAngVer);   break; //up
      case 6:                mat->RotateY(-kAngHor); mat->RotateZ(kAngVer);   break; //left and up 
    }
    mat->RotateZ(kAngCom);     //apply common rotation  in XY plane    
    gGeoManager->GetVolume("ALIC")->AddNode(hmp,ch,mat);
  }
  if(isOnlyChambers) return; //do not construct the detailed geometry  
  
  TGeoRotation *rot=new TGeoRotation("HwireRot"); rot->RotateY(90); //rotate wires around Y to be along X (initially along Z)
  TGeoVolume *sbo=gGeoManager->MakeBox ("Hsbo",ch4  , 1419*mm/2 , 1378.00*mm/2 ,   50.5*mm/2);//2072P1
  TGeoVolume *cov=gGeoManager->MakeBox ("Hcov",al   , 1419*mm/2 , 1378.00*mm/2 ,    0.5*mm/2);  
  TGeoVolume *hon=gGeoManager->MakeBox ("Hhon",roha , 1359*mm/2 , 1318.00*mm/2 ,   49.5*mm/2);  
  TGeoVolume *rad=gGeoManager->MakeBox ("Hrad",c6f14, 1330*mm/2 ,  413.00*mm/2 ,   24.0*mm/2); //2011P1
  TGeoVolume *neo=gGeoManager->MakeBox ("Hneo",neoc , 1330*mm/2 ,  413.00*mm/2 ,    4.0*mm/2); 
  TGeoVolume *win=gGeoManager->MakeBox ("Hwin",sio2 , 1330*mm/2 ,  413.00*mm/2 ,    5.0*mm/2); 
  TGeoVolume *si1=gGeoManager->MakeBox ("Hsi1",sio2 , 1330*mm/2 ,    5.00*mm/2 ,   15.0*mm/2);    
  TGeoVolume *si2=gGeoManager->MakeBox ("Hsi2",neoc ,   10*mm/2 ,  403.00*mm/2 ,   15.0*mm/2);    
  TGeoVolume *spa=gGeoManager->MakeTube("Hspa",sio2 ,    0*mm   ,    5.00*mm   ,   15.0*mm/2);         
  TGeoVolume *fr4=gGeoManager->MakeBox ("Hfr4",ch4  , 1407*mm/2 , 1366.00*mm/2 ,   15.0*mm/2);//2043P1 
  TGeoVolume *f4a=gGeoManager->MakeBox ("Hf4a",al   , 1407*mm/2 , 1366.00*mm/2 ,   10.0*mm/2); 
  TGeoVolume *f4i=gGeoManager->MakeBox ("Hf4i",ch4  , 1323*mm/2 , 1296.00*mm/2 ,   10.0*mm/2); 
  TGeoVolume *col=gGeoManager->MakeTube("Hcol",cu   ,    0*mm   ,  100.00*um   , 1323.0*mm/2);
  TGeoVolume *sec=gGeoManager->MakeBox ("Hsec",ch4  ,  648*mm/2 ,  411.00*mm/2 ,   45.5*mm/2);//sec=gap+ppf
  TGeoVolume *ppf=gGeoManager->MakeBox ("Hppf",al   ,  648*mm/2 ,  411.00*mm/2 ,   40.0*mm/2);//2001P2
  TGeoVolume *lar=gGeoManager->MakeBox ("Hlar",ar   ,  181*mm/2 ,   89.25*mm/2 ,   38.3*mm/2);//2001P2
  TGeoVolume *smo=gGeoManager->MakeBox ("Hsmo",ar   ,  114*mm/2 ,   89.25*mm/2 ,   38.3*mm/2);//2001P2
  TGeoVolume *gap=gGeoManager->MakeBox ("Hgap",ch4  ,  640*mm/2 ,  403.20*mm/2 ,    5.5*mm/2);//gap=pad+ano+cat
  TGeoVolume *cat=gGeoManager->MakeTube("Hcat",cu   ,    0*mm   ,   50.00*um   ,    8.0*mm/2); 
  TGeoVolume *ano=gGeoManager->MakeTube("Hano",w    ,    0*mm   ,   20.00*um   ,    8.0*mm/2); 
  TGeoVolume *pad=gGeoManager->MakeBox ("Hpad",csi  ,    8*mm/2 ,    8.40*mm/2 ,    1.0*mm/2);      
//
// ^ Y   z=         z=-12mm      z=98.25mm               ALIC->7xHmp (virtual)-->1xHsbo (virtual) --->2xHcov (real) 2072P1
// |  ____________________________________                                    |                   |-->1xHhon (real) 2072P1
// | |   ______     ____          ______  |                                   |
//   |  |      |   |    |   *    |      | |                                   |->3xHrad (virtual) --->1xHneo (real) 2011P1
//   |  |50.5mm|   |24mm|   *    |45.5mm| |                                   |                   |-->1xHwin (real) 2011P1
//   |  |      |   |    |   *    |      | |                                   |                   |-->2xHsi1 (real) 2011P1
//   |  |      |   |____|   *    |______| |                                   |                   |-->2xHsi2 (real) 2011P1
//   |  |      |    ____    *     ______  |                                   |                   |->30xHspa (real) 2011P1
//   |  |      |   |    |   *    |      | |                                   |
//   |  |      |   |    |   *    |      | |                                   |->1xHfr4 (vitual) --->1xHf4a (real)---->1xHf4i(virtual) 2043P1 
//   |  |  sb  |   | rad|   *    |      | |                                   |                  |-->322xHcol (real) 2043P1
//   |  |      |   |____|   *    |______| |                                   |
//   |  |      |    ____    *     ______  |                                   |->6xHsec (virtual) --> 1xHppf(real) ---->8xHlar (virtual) 2001P1
//   |  |      |   |    |   *    |      | |                                                                        |--->8xHsmo (virtual) 2001P1     
//   |  |      |   |    |   *    |      | |                                   |               
//   |  |      |   |    |   *    |      | |                                   |-> 1xHgap (virtual) --->48xHrow (virtual) -->80xHcel (virtual) -->4xHcat (real) from p84 TDR 
//   |  |______|   |____|   *    |______| |                                                                                                  |-->2xHano (real) from p84 TDR                                  
//   |____________________________________|                                                                                                  |-->1xHpad (real) from p84 TDR 
//                                                       --->Z 
  hmp->AddNode(sbo ,1,new TGeoTranslation(   0*mm,   0*mm, -73.75*mm));                     //p.84 TDR
     sbo->AddNode(hon ,1,new TGeoTranslation(  0*mm,0*mm,      0*mm)); //2072P1
     sbo->AddNode(cov ,1,new TGeoTranslation(  0*mm,0*mm,    +25*mm)); 
     sbo->AddNode(cov ,2,new TGeoTranslation(  0*mm,0*mm,    -25*mm)); 
  hmp->AddNode(rad,2,new TGeoTranslation(   0*mm,+434*mm, -12.00*mm)); 
  hmp->AddNode(rad,1,new TGeoTranslation(   0*mm,   0*mm, -12.00*mm)); 
  hmp->AddNode(rad,0,new TGeoTranslation(   0*mm,-434*mm, -12.00*mm)); 
    rad->AddNode(neo,1,new TGeoTranslation(   0*mm,   0*mm, -10.0*mm));
    rad->AddNode(win,1,new TGeoTranslation(   0*mm,   0*mm,   9.5*mm));
    rad->AddNode(si1,1,new TGeoTranslation(   0*mm,-204*mm,  -0.5*mm)); rad->AddNode(si1,2,new TGeoTranslation(   0*mm,+204*mm,  -0.5*mm));
    rad->AddNode(si2,1,new TGeoTranslation(-660*mm,   0*mm,  -0.5*mm)); rad->AddNode(si2,2,new TGeoTranslation(+660*mm,   0*mm,  -0.5*mm));
    for(Int_t i=0;i<3;i++) for(Int_t j=0;j<10;j++) rad->AddNode(spa,copy=10*i+j,new TGeoTranslation(-1330*mm/2+116*mm+j*122*mm,(i-1)*105*mm,-0.5*mm));
  hmp->AddNode(fr4,1,new TGeoTranslation(   0*mm,   0*mm,   9.00*mm));                     //p.84 TDR
  for(int i=1;i<=322;i++)  fr4->AddNode(col,i,new TGeoCombiTrans( 0*mm, -1296/2*mm+i*4*mm,-5*mm,rot)); //F4 2043P1
                           fr4->AddNode(f4a,1,new TGeoTranslation(   0*mm,0*mm, 2.5*mm));    
                                        f4a->AddNode(f4i,1,new TGeoTranslation(   0*mm,0*mm,   0*mm));
  hmp->AddNode(sec,4,new TGeoTranslation(-335*mm,+433*mm,  98.25*mm)); hmp->AddNode(sec,5,new TGeoTranslation(+335*mm,+433*mm,  98.25*mm));
  hmp->AddNode(sec,2,new TGeoTranslation(-335*mm,   0*mm,  98.25*mm)); hmp->AddNode(sec,3,new TGeoTranslation(+335*mm,   0*mm,  98.25*mm));
  hmp->AddNode(sec,0,new TGeoTranslation(-335*mm,-433*mm,  98.25*mm)); hmp->AddNode(sec,1,new TGeoTranslation(+335*mm,-433*mm,  98.25*mm));
    sec->AddNode(gap,1,new TGeoTranslation(0,0,-20.00*mm));
  TGeoVolume *row=          gap->Divide("Hrow",2,48,0,0);//along Y->48 rows
  TGeoVolume *cel=          row->Divide("Hcel",1,80,0,0);//along X->80 cells
      cel->AddNode(cat,1,new TGeoCombiTrans (0,  3.15*mm , -2.70*mm , rot)); //4 cathode wires
      cel->AddNode(ano,1,new TGeoCombiTrans (0,  2.00*mm , -0.29*mm , rot)); //2 anod wires
      cel->AddNode(cat,2,new TGeoCombiTrans (0,  1.05*mm , -2.70*mm , rot)); 
      cel->AddNode(cat,3,new TGeoCombiTrans (0, -1.05*mm , -2.70*mm , rot)); 
      cel->AddNode(ano,2,new TGeoCombiTrans (0, -2.00*mm , -0.29*mm , rot)); 
      cel->AddNode(cat,4,new TGeoCombiTrans (0, -3.15*mm , -2.70*mm , rot));   
      cel->AddNode(pad,1,new TGeoTranslation(0,  0.00*mm ,  2.25*mm));       //1 pad  
    sec->AddNode(ppf,1,new TGeoTranslation(0,0,  2.75*mm));
// ^ Y  single cell                                                5.5mm CH4 = 1*mm CsI + 4.45*mm CsI x cath +0.05*mm safety margin         
// |      ______________________________           
// |     |                              |          ^                            ||     
//       |                              | 1.05mm                                ||     
// 2.2*mm| xxxxxxxxxxxxxxxxxxxxxxxxxxxx |--              50um  x                || cat shift  x=0mm , y= 3.15mm , z=-2.70mm       
//       |                              |                                       ||     
//       |                              |                                       ||     
// __    |  ..........................  | 2.1mm                    20un .       ||  ano shift x=0mm , y= 2.00mm , z=-0.29mm   
//       |                              |                                       ||     
//       |                              |                                       ||     
//       | xxxxxxxxxxxxxxxxxxxxxxxxxxxx |--                    x                ||  cat shift x=0mm , y= 1.05mm , z=-2.70mm   
//       |                              |                                       ||     
//       |                              |         8.4mm                         ||   
// 4*mm  |                              | 2.1mm                                 ||  pad shift x=0mm , y= 0.00mm , z=2.25*mm   
//       |                              |                                       ||  
//       |                              |                                       ||  
//       | xxxxxxxxxxxxxxxxxxxxxxxxxxxx |--                    x                ||  cat shift x=0mm , y=-1.05mm , z=-2.70mm   
//       |                              |                                       ||  
//       |                              |                                       ||    
// __    |  ..........................  | 2.1mm                         . 2.04mm||  ano shift x=0mm , y=-2.00mm , z=-0.29mm   
//       |                              |                                       ||  
//       |                              |                                       ||  
//       | xxxxxxxxxxxxxxxxxxxxxxxxxxxx |--                    x    4.45mm      ||  cat shift x=0mm , y=-3.15mm , z=-2.70mm   
// 2.2*mm|                              |                                       ||  
//       |                              | 1.05mm                                ||         
//       |______________________________|          v                            ||    
//       <             8 mm             >                          
//                                   ----->X                                 ----->Z
  ppf->AddNode(lar,1,new TGeoTranslation(-224.5*mm,-151.875*mm,  0.85*mm));
  ppf->AddNode(lar,2,new TGeoTranslation(-224.5*mm,- 50.625*mm,  0.85*mm));
  ppf->AddNode(lar,3,new TGeoTranslation(-224.5*mm,+ 50.625*mm,  0.85*mm));
  ppf->AddNode(lar,4,new TGeoTranslation(-224.5*mm,+151.875*mm,  0.85*mm));
  ppf->AddNode(lar,5,new TGeoTranslation(+224.5*mm,-151.875*mm,  0.85*mm));
  ppf->AddNode(lar,6,new TGeoTranslation(+224.5*mm,- 50.625*mm,  0.85*mm));
  ppf->AddNode(lar,7,new TGeoTranslation(+224.5*mm,+ 50.625*mm,  0.85*mm));
  ppf->AddNode(lar,8,new TGeoTranslation(+224.5*mm,+151.875*mm,  0.85*mm));
  ppf->AddNode(smo,1,new TGeoTranslation(- 65.0*mm,-151.875*mm,  0.85*mm));
  ppf->AddNode(smo,2,new TGeoTranslation(- 65.0*mm,- 50.625*mm,  0.85*mm));
  ppf->AddNode(smo,3,new TGeoTranslation(- 65.0*mm,+ 50.625*mm,  0.85*mm));
  ppf->AddNode(smo,4,new TGeoTranslation(- 65.0*mm,+151.875*mm,  0.85*mm));
  ppf->AddNode(smo,5,new TGeoTranslation(+ 65.0*mm,-151.875*mm,  0.85*mm));
  ppf->AddNode(smo,6,new TGeoTranslation(+ 65.0*mm,- 50.625*mm,  0.85*mm));
  ppf->AddNode(smo,7,new TGeoTranslation(+ 65.0*mm,+ 50.625*mm,  0.85*mm));
  ppf->AddNode(smo,8,new TGeoTranslation(+ 65.0*mm,+151.875*mm,  0.85*mm)); 
}//Hmpid()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Axis()
{
// Draw axises  on top of geometry
  Double_t X[6]={0,0,0,300,0,0};  Double_t Y[6]={0,0,0,0,300,0};  Double_t Z[6]={0,0,0,0,0,300};
  TPolyLine3D *pXaxis=new TPolyLine3D(2,X);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,Y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,Z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
