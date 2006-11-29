
Int_t copy; //volume copy number 
Double_t dx,dy,dz,r1,r2;//tmp vars for volume dimentions
Double_t cm=1,m=100*cm,mm=0.1*cm,mkm=0.001*cm;//length units  
TGeoManager *g=0;

void RichGeom(Bool_t isOnlyChambers=kFALSE)
{
  
  g=new TGeoManager("HMPID","Private HMPID geometry");
  Materials();  
  gGeoManager->MakeBox("ALIC",gGeoManager->GetMedium("Air"),dx=30*m/2,dy=30*m/2,dz=30*m/2); //arbitrary values  
  gGeoManager->SetTopVolume(gGeoManager->GetVolume("ALIC"));
    
  Rich(isOnlyChambers);
  
//  RusGel();
  
//  Vhmpid();
  
  gGeoManager->CloseGeometry();
  
  gGeoManager->SetVisOption(0); gGeoManager->SetVisLevel(5); 
  
  gGeoManager->GetMasterVolume()->Draw();
  Axis();
//  gPad->GetView()->SetView(3,94,-70,0);
  new TBrowser;
}
//__________________________________________________________________________________________________
void Materials()
{
//Media for HMPID
  Double_t a=0,z=0,den=0,radlen=0,intlen=0;//tmp vars for material parameters
  new TGeoMaterial("Air"           ,a=26.98 ,z=13,den=0.4224                             ); new TGeoMedium("Air"          ,1,gGeoManager->GetMaterial("Air")); 
  new TGeoMaterial("HMPID_CH4"      ,a=26.98 ,z=13,den=0.4224                             ); new TGeoMedium("HMPID_CH4"     ,2,gGeoManager->GetMaterial("HMPID_CH4"));
  new TGeoMaterial("HMPID_CsI"      ,a=26.98 ,z=13,den=2.7 ,radlen=24.01*cm,intlen=70.6*cm); new TGeoMedium("HMPID_CsI"     ,3,gGeoManager->GetMaterial("HMPID_CsI"));
  new TGeoMaterial("HMPID_Al"       ,a=26.98 ,z=13,den=2.7 ,radlen=24.01*cm,intlen=70.6*cm); new TGeoMedium("HMPID_Al"      ,4,gGeoManager->GetMaterial("HMPID_Al")); 
  new TGeoMaterial("HMPID_W"        ,a=183.84,z=27,den=19.3,radlen= 9.59*cm,intlen=0.35*cm); new TGeoMedium("HMPID_W"       ,5,gGeoManager->GetMaterial("HMPID_W"));
  new TGeoMaterial("HMPID_Cu"       ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); new TGeoMedium("HMPID_Cu"      ,6,gGeoManager->GetMaterial("HMPID_Cu"));
  new TGeoMaterial("HMPID_Rohacell" ,a=12.01 ,z=6 ,den=0.1 ,radlen=18.8,intlen=0);           new TGeoMedium("HMPID_Rohacell",7,gGeoManager->GetMaterial("HMPID_Rohacell"));
  new TGeoMaterial("HMPID_SiO2"     ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_SiO2"    ,8,gGeoManager->GetMaterial("HMPID_SiO2")); 
  new TGeoMaterial("HMPID_C6F14"    ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_C6F14"   ,9,gGeoManager->GetMaterial("HMPID_C6F14"));
//Media for Sr90 source  
  new TGeoMaterial("HMPID_Perpex"   ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); new TGeoMedium("HMPID_Perpex"  ,10,gGeoManager->GetMaterial("HMPID_Perpex"));
  new TGeoMaterial("HMPID_Steel"    ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); new TGeoMedium("HMPID_Steel"   ,11,gGeoManager->GetMaterial("HMPID_Steel"));
  new TGeoMaterial("HMPID_Mylar"    ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); new TGeoMedium("HMPID_Mylar"   ,12,gGeoManager->GetMaterial("HMPID_Mylar"));
  new TGeoMaterial("HMPID_Sr90"     ,a=87.62 ,z=38,den=7.87,radlen=13.84*cm,intlen=82.8*cm); new TGeoMedium("HMPID_Sr90"    ,13,gGeoManager->GetMaterial("HMPID_Sr90"));
//Media for VHMPID Gas option  
  new TGeoMaterial("HMPID_CF4"      ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_CF4"     ,14,gGeoManager->GetMaterial("HMPID_CF4"));
  new TGeoMaterial("HMPID_C4F10"    ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_C4F10"   ,15,gGeoManager->GetMaterial("HMPID_C4F10"));
  new TGeoMaterial("HMPID_Ag"       ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Ag"      ,16,gGeoManager->GetMaterial("HMPID_Ag"));
//Media for VHMPID aerogel option  
  new TGeoMaterial("HMPID_Gel24"    ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Gel24"   ,17,gGeoManager->GetMaterial("HMPID_Gel24"));
  new TGeoMaterial("HMPID_Gel26"    ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Gel26"   ,18,gGeoManager->GetMaterial("HMPID_Gel26"));
  new TGeoMaterial("HMPID_Gel28"    ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Gel28"   ,19,gGeoManager->GetMaterial("HMPID_Gel28"));
  new TGeoMaterial("HMPID_Gel30"    ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Gel30"   ,20,gGeoManager->GetMaterial("HMPID_Gel30"));
  new TGeoMaterial("HMPID_Si"       ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Si"      ,21,gGeoManager->GetMaterial("HMPID_Si"));
  new TGeoMaterial("HMPID_Apd"      ,a=0     ,z=0 ,den=0);                                   new TGeoMedium("HMPID_Apd"     ,22,gGeoManager->GetMaterial("HMPID_Apd"));
}//Materials()
//__________________________________________________________________________________________________
void Rich(Bool_t isOnlyChambers)
{
//Rich  chamber
  TGeoVolume *pRich=gGeoManager->MakeBox("HMPID",gGeoManager->GetMedium("Air"),dx=(6*mm+1681*mm+6*mm)/2,  //main HMPID volume
                                                                              dy=(6*mm+1466*mm+6*mm)/2,
                                                                              dz=(80*mm+40*mm)*2/2);     //x,y taken from 2033P1  z from p84 TDR  
  const Double_t kAngHor=19.5; //  horizontal angle between chambers  19.5 grad
  const Double_t kAngVer=20;   //  vertical angle between chambers    20   grad     
  const Double_t kAngCom=30;   //  common HMPID rotation with respect to x axis  30   grad     
  const Double_t trans[3]={490,0,0}; //center of the chamber is on window-gap surface
  for(Int_t iCh=1;iCh<=7;iCh++){//place 7 chambers
    TGeoHMatrix *pMatrix=new TGeoHMatrix;
    pMatrix->RotateY(90);           //rotate around y since initial position is in XY plane -> now in YZ plane
    pMatrix->SetTranslation(trans); //now plane in YZ is shifted along x 
    switch(iCh){
      case 1:                pMatrix->RotateY(kAngHor);  pMatrix->RotateZ(-kAngVer);  break; //right and down 
      case 2:                                            pMatrix->RotateZ(-kAngVer);  break; //down              
      case 3:                pMatrix->RotateY(kAngHor);                               break; //right 
      case 4:                                                                         break; //no rotation
      case 5:                pMatrix->RotateY(-kAngHor);                              break; //left   
      case 6:                                            pMatrix->RotateZ(kAngVer);   break; //up
      case 7:                pMatrix->RotateY(-kAngHor); pMatrix->RotateZ(kAngVer);   break; //left and up 
    }
    pMatrix->RotateZ(kAngCom);     //apply common rotation  in XY plane    
    gGeoManager->GetVolume("ALIC")->AddNode(pRich,iCh,pMatrix);
  }
  if(isOnlyChambers) return; //do not construct the detailed geometry  
//Pad Panel frame  
  TGeoVolume *pPpf     =gGeoManager->MakeBox("Rppf"  ,gGeoManager->GetMedium("HMPID_Al")  ,dx=648*mm/2,dy=  411*mm/2 ,dz=40  *mm/2);//PPF 2001P2 inner size of the slab by 1mm more
  TGeoVolume *pPpfLarge=gGeoManager->MakeBox("Rppf1" ,gGeoManager->GetMedium("Air")      ,dx=181*mm/2,dy=89.25*mm/2 ,dz=38.3*mm/2);     //large whole
  TGeoVolume *pPpfSmall=gGeoManager->MakeBox("Rppf2" ,gGeoManager->GetMedium("Air")      ,dx=114*mm/2,dy=89.25*mm/2 ,dz=38.3*mm/2);//small whole
  TGeoVolume *pPc      =gGeoManager->MakeBox("Rpc"   ,gGeoManager->GetMedium("HMPID_CsI") ,dx=644*mm/2,dy=  407*mm/2 ,dz= 1.7*mm/2);//by 0.2 mm more then actual size (PCB 2006P1)
  
  pRich->AddNode(pPpf,copy=1,new TGeoTranslation(-335*mm,-433*mm,8*cm+20*mm));//F1 2040P1 z p.84 TDR
  pRich->AddNode(pPpf,copy=2,new TGeoTranslation(+335*mm,-433*mm,8*cm+20*mm));
  pRich->AddNode(pPpf,copy=3,new TGeoTranslation(-335*mm,   0*mm,8*cm+20*mm));
  pRich->AddNode(pPpf,copy=4,new TGeoTranslation(+335*mm,   0*mm,8*cm+20*mm));
  pRich->AddNode(pPpf,copy=5,new TGeoTranslation(-335*mm,+433*mm,8*cm+20*mm));
  pRich->AddNode(pPpf,copy=6,new TGeoTranslation(+335*mm,+433*mm,8*cm+20*mm));
  pPpf->AddNode( pPc ,copy=1,new TGeoTranslation(   0*mm,   0*mm,-19.15*mm));//PPF 2001P2 
  pPpf->AddNode(pPpfLarge,copy=1,new TGeoTranslation(-224.5*mm,-151.875*mm,  0.85*mm));
  pPpf->AddNode(pPpfLarge,copy=2,new TGeoTranslation(-224.5*mm,- 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfLarge,copy=3,new TGeoTranslation(-224.5*mm,+ 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfLarge,copy=4,new TGeoTranslation(-224.5*mm,+151.875*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=1,new TGeoTranslation(- 65.0*mm,-151.875*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=2,new TGeoTranslation(- 65.0*mm,- 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=3,new TGeoTranslation(- 65.0*mm,+ 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=4,new TGeoTranslation(- 65.0*mm,+151.875*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=5,new TGeoTranslation(+ 65.0*mm,-151.875*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=6,new TGeoTranslation(+ 65.0*mm,- 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=7,new TGeoTranslation(+ 65.0*mm,+ 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfSmall,copy=8,new TGeoTranslation(+ 65.0*mm,+151.875*mm,  0.85*mm)); 
  pPpf->AddNode(pPpfLarge,copy=5,new TGeoTranslation(+224.5*mm,-151.875*mm,  0.85*mm));
  pPpf->AddNode(pPpfLarge,copy=6,new TGeoTranslation(+224.5*mm,- 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfLarge,copy=7,new TGeoTranslation(+224.5*mm,+ 50.625*mm,  0.85*mm));
  pPpf->AddNode(pPpfLarge,copy=8,new TGeoTranslation(+224.5*mm,+151.875*mm,  0.85*mm));
//gap - anod wires
  TGeoVolume *pGap =gGeoManager->MakeBox ("Rgap" ,gGeoManager->GetMedium("HMPID_CH4") ,dx=648*mm/2,dy=  411*mm/2 ,dz=4.45*mm/2);//xy as PPF 2001P2 z WP 2099P1
  TGeoVolume *pAnod=gGeoManager->MakeTube("Rano" ,gGeoManager->GetMedium("HMPID_W")   ,r1=  0*mm  ,r2=  20*mkm/2 ,dz=648*mm/2); //WP 2099P1 z = gap x PPF 2001P2
  TGeoRotation *pAnodRot=new TGeoRotation("RanW",90,90,0);
  
  pRich->AddNode(pGap,copy=1,new TGeoTranslation(-335*mm,-433*mm,8*cm-2.225*mm)); //F1 2040P1 z WP 2099P1
  pRich->AddNode(pGap,copy=2,new TGeoTranslation(+335*mm,-433*mm,8*cm-2.225*mm)); 
  pRich->AddNode(pGap,copy=3,new TGeoTranslation(-335*mm,   0*mm,8*cm-2.225*mm)); 
  pRich->AddNode(pGap,copy=4,new TGeoTranslation(+335*mm,   0*mm,8*cm-2.225*mm)); 
  pRich->AddNode(pGap,copy=5,new TGeoTranslation(-335*mm,+433*mm,8*cm-2.225*mm)); 
  pRich->AddNode(pGap,copy=6,new TGeoTranslation(+335*mm,+433*mm,8*cm-2.225*mm)); 
  for(int i=1;i<=96;i++)
    pGap->AddNode(pAnod,copy=i,new TGeoCombiTrans( 0*mm, -411/2*mm+i*4*mm, 0.185*mm,pAnodRot)); //WP 2099P1  
//frame 3- cathode wires      
  TGeoVolume *pCath=gGeoManager->MakeTube("RcaW"  ,gGeoManager->GetMedium("Cu")  ,r1=0  ,r2=100*mkm/2,dz=1323*mm/2);//r WP 2099P1 z F3 2041P1       
  TGeoRotation *pCathRot=new TGeoRotation("CathRot",90,90,0);
  for(int i=1;i<=618;i++)
    pRich->AddNode(pCath,copy=i,new TGeoCombiTrans( 0*mm, -649.5*mm+i*2.1*mm, 75*mm,pCathRot)); //WP 2099P1    
//Frame 4- collection wires
  TGeoVolume *pF4  =gGeoManager->MakeBox( "Rfr4"       ,gGeoManager->GetMedium("HMPID_CH4")   ,dx=1407*mm/2 ,dy=1366*mm/2  ,dz=  15*mm/2);//F4 2043P1 
  TGeoVolume *pF4al=gGeoManager->MakeBox( "Rfr4al"     ,gGeoManager->GetMedium("HMPID_Al")    ,dx=1407*mm/2 ,dy=1366*mm/2  ,dz=  10*mm/2); 
  TGeoVolume *pF4in=gGeoManager->MakeBox( "Rfr4in"     ,gGeoManager->GetMedium("HMPID_CH4")   ,dx=1323*mm/2 ,dy=1296*mm/2  ,dz=  10*mm/2); 
  TGeoVolume *pColl=gGeoManager->MakeTube("RcoW"       ,gGeoManager->GetMedium("HMPID_Cu")    ,r1=   0*mm   ,r2=100*mkm/2  ,dz=1323*mm/2);
  TGeoRotation *pCollRot=new TGeoRotation("RcoRot",90,90,0);
  
  pRich->AddNode(pF4      ,copy=1,new TGeoTranslation(   0*mm,0*mm,   9*mm)); //F4 to Rich p.84 TDR
    pF4  ->AddNode(pF4al    ,copy=1,new TGeoTranslation(   0*mm,0*mm, 2.5*mm)); //F4 al to F4 2043P1 
    pF4al->AddNode(pF4in    ,copy=1,new TGeoTranslation(   0*mm,0*mm,   0*mm)); //F4 whole F4 al 2043P1   
    for(int i=1;i<=322;i++)
      pF4->AddNode(pColl,copy=i,new TGeoCombiTrans( 0*mm, -1296/2*mm+i*4*mm, -5*mm,pCollRot)); //F4 2043P1
//radiators
  TGeoVolume *pRad      =gGeoManager->MakeBox( "Rad"      ,gGeoManager->GetMedium("HMPID_C6F14")    ,dx=1330*mm/2 ,dy= 413*mm/2  ,dz=  24*mm/2); // Rad 2011P1
  TGeoVolume *pRadFront =gGeoManager->MakeBox( "RadFront" ,gGeoManager->GetMedium("HMPID_Neoceram") ,dx=1330*mm/2 ,dy= 413*mm/2  ,dz=   4*mm/2); 
  TGeoVolume *pRadWin   =gGeoManager->MakeBox( "RadWin"   ,gGeoManager->GetMedium("HMPID_SiO2")     ,dx=1330*mm/2 ,dy= 413*mm/2  ,dz=   5*mm/2); 
  TGeoVolume *pRadLong  =gGeoManager->MakeBox( "RadLong"  ,gGeoManager->GetMedium("HMPID_Neoceram") ,dx=1330*mm/2 ,dy=   5*mm/2  ,dz=  15*mm/2);    
  TGeoVolume *pRadShort =gGeoManager->MakeBox( "RadShort" ,gGeoManager->GetMedium("HMPID_Neoceram") ,dx=  10*mm/2 ,dy= 403*mm/2  ,dz=  15*mm/2);    
  TGeoVolume *pRadSpacer=gGeoManager->MakeTube("RadSpacer",gGeoManager->GetMedium("HMPID_SiO2")     ,r1= 0        ,r2=10*mm/2  ,dz=  15*mm/2);         
    
  pRich->AddNode(pRad      ,copy=1,new TGeoTranslation(   0*mm,-434*mm,   -12*mm)); 
  pRich->AddNode(pRad      ,copy=2,new TGeoTranslation(   0*mm,   0*mm,   -12*mm)); 
  pRich->AddNode(pRad      ,copy=3,new TGeoTranslation(   0*mm,+434*mm,   -12*mm)); 
    
  pRad ->AddNode(pRadFront ,copy=1,new TGeoTranslation(   0*mm,   0*mm, -10.0*mm));
  pRad ->AddNode(pRadWin   ,copy=1,new TGeoTranslation(   0*mm,   0*mm,   9.5*mm));
  pRad ->AddNode(pRadLong  ,copy=1,new TGeoTranslation(   0*mm,-204*mm,  -0.5*mm));
  pRad ->AddNode(pRadLong  ,copy=2,new TGeoTranslation(   0*mm,+204*mm,  -0.5*mm));
  pRad ->AddNode(pRadShort ,copy=1,new TGeoTranslation(-660*mm,   0*mm,  -0.5*mm));
  pRad ->AddNode(pRadShort ,copy=2,new TGeoTranslation(+660*mm,   0*mm,  -0.5*mm));
  for(int i=0;i<3;i++) for(int j=0;j<10;j++)  pRad->AddNode(pRadSpacer,copy=10*i+j,new TGeoTranslation(-1330*mm/2+116*mm+j*122*mm,(i-1)*105*mm,-0.5*mm));
//sandbox  
  TGeoVolume *pSandBox  =gGeoManager->MakeBox( "RSandBox"  ,gGeoManager->GetMedium("Air")  ,dx=1419*mm/2 ,dy=1378*mm/2   ,dz=50.5*mm/2);  //2072P1   
  TGeoVolume *pSandCover=gGeoManager->MakeBox( "RSandCover",gGeoManager->GetMedium("HMPID_Al")       ,dx=1419*mm/2 ,dy=1378*mm/2   ,dz= 0.5*mm/2);  
  TGeoVolume *pSandComb =gGeoManager->MakeBox( "RSandComb" ,gGeoManager->GetMedium("HMPID_Rohacell") ,dx=1359*mm/2 ,dy=1318*mm/2   ,dz=49.5*mm/2);  
  
  pRich->AddNode(pSandBox,copy=1,new TGeoTranslation(  0*mm,0*mm, -73.75*mm)); //p.84 TDR
    pSandBox->AddNode(pSandComb  ,copy=1,new TGeoTranslation(  0*mm,0*mm,      0*mm)); //2072P1
    pSandBox->AddNode(pSandCover ,copy=1,new TGeoTranslation(  0*mm,0*mm,    +25*mm)); 
    pSandBox->AddNode(pSandCover ,copy=2,new TGeoTranslation(  0*mm,0*mm,    -25*mm)); 
}//Rich()
//__________________________________________________________________________________________________
void Sr90(TGeoVolume *pTop)
{
    pSrc               =gGeoManager->MakeTube("Src"              ,gGeoManager->GetMedium("Air")        , 0 , 70*mm/2 ,  30*mm/2);       //top container
      pAlGlass         =gGeoManager->MakeTube("SrcAlGlass"       ,gGeoManager->GetMedium("HMPID_Al")    , 0 , 38*mm/2 ,21.8*mm/2);       //Al glass wall        
        pPerpexPlug    =gGeoManager->MakeTube("SrcPerpex"        ,gGeoManager->GetMedium("HMPID_Perpex"), 0 , 34*mm/2 ,  20*mm/2);       //Perpex plug         
          pScrewCentral=gGeoManager->MakeTube("SrcScrewCentral"  ,gGeoManager->GetMedium("HMPID_Steel") , 0 ,  5*mm/2 ,  15*mm/2);       //Steel screw in the center        
          pScrewSr90   =gGeoManager->MakeTube("SrcScrewSr90"     ,gGeoManager->GetMedium("HMPID_Steel") , 0 ,  2*mm/2 ,  10*mm/2);       //Steel screw to support Sr90 
            pSr90      =gGeoManager->MakeTube("SrcSr90"          ,gGeoManager->GetMedium("HMPID_Sr90")  , 0 ,  1*mm/2 ,   1*mm/2);       //Sr90 source
          pHolePerpex  =gGeoManager->MakeTube("SrcHolePerpex"    ,gGeoManager->GetMedium("Air")        , 0 ,  4*mm/2 ,  10*mm/2);       //Air hole in perpex plug      
        pHoleAl        =gGeoManager->MakeTube("SrcHoleAl"        ,gGeoManager->GetMedium("Air")        , 0 ,  5*mm/2 , 1.8*mm/2);       //Air hole in Al glass bottom
    pMylarFoil         =gGeoManager->MakeTube("SrcMylarFoil"     ,gGeoManager->GetMedium("HMPID_Mylar") , 0 , 30*mm/2 , 50*mkm/2);       //Mylar foil                
                
    pTop->AddNode(pSrc,1,new TGeoTranslation(30*cm,0,1*cm));
      pSrc ->AddNode(pMylarFoil,1,new TGeoTranslation(0,0,21.8*mm/2+50*mkm/2));
      pSrc ->AddNode(pAlGlass,1,new TGeoTranslation(0,0,0));                           //Al glass to fake Src volume
        pAlGlass->AddNode(       pHoleAl      ,1,new TGeoTranslation(6*mm,0, -10*mm));
        pAlGlass->AddNode(       pPerpexPlug  ,1,new TGeoTranslation(0*mm,0, 0.9*mm));
          pPerpexPlug->AddNode(  pHolePerpex  ,1,new TGeoTranslation(6*mm,0,  -5*mm));      
          pPerpexPlug->AddNode(  pScrewCentral,1,new TGeoTranslation(0   ,0, 2.5*mm));  
          pPerpexPlug->AddNode(  pScrewSr90   ,1,new TGeoTranslation(6*mm,0,   5*mm));  
            pScrewSr90->AddNode( pSr90        ,1,new TGeoTranslation(0   ,0,-4.5*mm));  
}//Sr90()    
//__________________________________________________________________________________________________
void RusGel()
{
//Defines VHMPID aerogel option geometry. 
//                 top view normal position                                   side view normal position
//                                                                                                 ^ y
//                  --------                                                                       | 
//                  |      |                                                                       |
//                  |      |          z<-----* y                                          z<-------* x  ----> MUON side
//                  |      |                 |                                  ----                    
//                  |      |                 |                                  |  | 
//                  ________                 v x                                |  |
//                                                                              |--|
//                                                                              |  |
//                                                                              ----       
//Chamber consists from Al box filled with air where 4 aerogel blocks and APD wall are positioned.
//  ------------------------------
//  |-|           |-| |-| |-| |-||
//  |-|           | | | | | | | ||  
//  |-|           | | | | | | | ||
//  |-| APD wall  | | | | | | | ||    z<---* y  top view   
//  |-|           | | | | | | | ||         |
//  |-|           | | | | | | | ||         | 
//  |-|           |_| |_| |_| |_||         v x
//  ------------------------------          
// 
//             ALIC
//               |
//             Vbox (Al)
//               |
//             Vair (Air) 
//        _______|________
//        |              |
//       4*Vgel         Vwall
//                       |
//                      Vcolumn (division along X)
//                       |
//                      Vcell   (division along Y)
//                       |
//                      Vapd

  Double_t cm=1 , m=100 , mm=0.1                                                     ;//dimentions, default is cm
  
  Int_t    iNapdsX        =10                                    ;//number of APDs along x
  Int_t    iNapdsY        =16                                    ;//number of APDs along y
  Double_t dCellX         =1.5*mm             *0.5                ;//cell X half size
  Double_t dCellY         =1.5*mm             *0.5                ;//cell Y half size
  Double_t dCellZ         =0.5*mm             *0.5                ;//APD wall thickness   
  Double_t dWallX         =                         iNapdsX*dCellX;//APD wall X half size
  Double_t dWallY         =                         iNapdsY*dCellY;//APD wall Y half size
  Double_t dWallZ         =                                 dCellZ;//APD wall half thickness  
  Double_t dApdR          =0.5*mm                                 ;//APD radius 
  Double_t dApdZ          =                                 dCellZ;//APD Z half size
  Double_t dGelX          =                                 dWallX;//gel block X half size
  Double_t dGelY          =                                 dWallY;//gel block Y half size
  Double_t dGelZ          =10*mm              *0.5                ;//gel block Z half size
  Double_t dProxGap       =50*cm                                  ;//half distance between APD wall and last aerogel block
  Double_t dAirX          =                                 dWallX;//internal air X hald size
  Double_t dAirY          =                                 dWallY;//internal air Y hald size
  Double_t dAirZ          =                dWallZ+dProxGap+7*dGelZ;//internal air Z hald size
  Double_t dBoxWall       =2*mm               *0.5                ;//Al box walls thickness
  Double_t dBoxX          =                         dAirX+dBoxWall;//Al box x half size 
  Double_t dBoxY          =                         dAirY+dBoxWall;//Al box y half size 
  Double_t dBoxZ          =                         dAirZ+dBoxWall;//Al box z half size
  
  Int_t copy;    Double_t rmin,rmax,dx,dy,dz;
//make external Al box
  TGeoVolume *pBox=gGeoManager->MakeBox("Gbox",gGeoManager->GetMedium("HMPID_Al"),dx=dBoxX,dy=dBoxY,dz=dBoxZ);
  TGeoRotation *pRot=new TGeoRotation("GboxRot"); pRot->RotateX(90);
  gGeoManager->GetVolume("ALIC")->AddNode(pBox,copy=1,new TGeoCombiTrans(0*m,-5.2*m,2.5*m,pRot));//normal position
//position Air to Al box   
  TGeoVolume *pAir=gGeoManager->MakeBox( "Gair",gGeoManager->GetMedium("Air"),dx=dAirX,dy=dAirY,dz=dAirZ);                   
  pBox->AddNode(pAir,copy=1); 
//position 4 gel blocks to Air
  TGeoVolume *pGel24=gGeoManager->MakeBox( "Ggel24",gGeoManager->GetMedium("HMPID_Gel24"),dx=dGelX,dy=dGelY,dz=dGelZ);     
    pAir->AddNode(pGel24,copy=1,new TGeoTranslation(0,0,-dAirZ+1*dGelZ)); 
  TGeoVolume *pGel26=gGeoManager->MakeBox( "Ggel26",gGeoManager->GetMedium("HMPID_Gel26"),dx=dGelX,dy=dGelY,dz=dGelZ);     
    pAir->AddNode(pGel26,copy=1,new TGeoTranslation(0,0,-dAirZ+5*dGelZ)); 
  TGeoVolume *pGel28=gGeoManager->MakeBox( "Ggel28",gGeoManager->GetMedium("HMPID_Gel28"),dx=dGelX,dy=dGelY,dz=dGelZ);     
    pAir->AddNode(pGel28,copy=1,new TGeoTranslation(0,0,-dAirZ+9*dGelZ)); 
  TGeoVolume *pGel30=gGeoManager->MakeBox( "Ggel30",gGeoManager->GetMedium("HMPID_Gel30"),dx=dGelX,dy=dGelY,dz=dGelZ);     
    pAir->AddNode(pGel30,copy=1,new TGeoTranslation(0,0,-dAirZ+13*dGelZ)); 
//position APD wall to air
  TGeoVolume   *pWall     =gGeoManager->MakeBox ("Gwall",gGeoManager->GetMedium("HMPID_Si"),dx=dWallX , dy=dWallY , dz=dWallZ );  
  pAir->AddNode(pWall,copy=1,new TGeoTranslation(0,0,dAirZ-dWallZ)); 
//divide wall into cells
  Int_t axis,ndiv; Double_t start,step;
  TGeoVolume *pWallCol =pWall      ->Divide("Gcol",axis=1,ndiv=iNapdsX,start=0,step=0);//divide VhGap along X by NpadsX columns
  TGeoVolume *pWallCell=pWallCol   ->Divide("Gcel",axis=2,ndiv=iNapdsY,start=0,step=0);//divide VhGapCol along Y by NpadsY cells
//position APD to wall cell
  TGeoVolume *pApd=gGeoManager->MakeTube("Gapd",gGeoManager->GetMedium("HMPID_Apd"),rmin=0,rmax=dApdR,dz=dApdZ); pWallCell->AddNode(pApd,copy=1); 
}//RusGel()
//__________________________________________________________________________________________________
void Vhmpid()
{
//Defines VHMPID geometry for TIC option.
//                 top view normal position                                   side view normal position
//                                                                                                 ^ y
//                  --------                                                                       | 
//                  |------|                                                                       |
//                  |      |          z<-----* y                                          z<-------* x  ----> MUON side
//                  |      |                 |                                  ----                    
//                  |      |                 |                                  |  | 
//                  ________                 v x                                |  |
//                                                                              |--|
//                                                                              |  |
//                                                                              ----       
//Chamber consists from Al box filled with radiator CF4 and C4F10 quartz window in between , Al mirror and MWPC.
//  ---------------------------------------------------  top view chamber in test position
//  |  -------------- MWPC     | quartz window        |
//  |   . . . . . .            |                      |  
//  |   \                      |                      |
//  |    \             CF4     |    C4F10             |  z<-----* y   
//  |     \ mirror             |                      |         |
//  |      \                   |                      |         | 
//  |       \                  |                      |         v x  
//  ---------------------------------------------------    
//                                                         z ^
//                                                           |
//                                                           | 
//    <-Y0-> X <--Y1--> X <--Y1--> X <--Y1--> X <-Y0->       |            cath wires: r 50mkm; shift Y0=1.05m;, pitch Y1=2.1; center to PC  4.45mm; material Cu
//                                                         x *------->y            
//        
//    <--Y0--> x <-----------Y1-----------> x <--Y0-->                    anod wires: r 20mkm; shift Y0=2.2mm; pitch Y1=4.0mm; center to PC  2.04mm; material W
//                                                                       
//                                                                                                
//   |________________________________________________|                   pad size y 8.4mm
// 
//             ALIC
//               |
//             Vbox
//        _______|______
//       |       |      |
//     Vc4f    Vwin   Vcf4
//                 _____|________
//                |              |
//              Vmir             |
//                               |
//                             Vgap
//                               |
//                             Vcol (column of gap cells) X
//                               |
//                             Vcel  cell in the column   Y                  
//                         ______|______
//                        |      |      |
//                      Vpad   Vano   Vcat
  Double_t cm=1 , m=100 , mm=0.1 , um=1e-4                                                      ;//dimentions, default is cm
  
  Int_t    iNpadsX        =                                    AliHMPIDParam::NpadsX()           ;//number of pads along x parametrised
  Int_t    iNpadsY        =                                    AliHMPIDParam::NpadsY()           ;//number of pads along y parametrised
  Double_t wCathR         =50  *um                                                              ;//cathode wire radius defined by USER
  Double_t wCathShift     =1.05*mm                                                              ;//cathode wire shift from pad edge defined by USER
  Double_t wCathPitch     =2.1 *mm                                                              ;//cathode wire pitch defined by USER
  Double_t wCathPc        =4.45*mm                                                              ;//distance from pc to cathode wire defined by USER
  Double_t wAnodR         =20  *um                                                              ;//anod wire radius defined by USER
  Double_t wAnodShift     =2.2 *mm                                                              ;//anod wire shift from pad edge defined by USER
  Double_t wAnodPc        =2.04*mm                                                              ;//distance from anod wire center to pc defined by USER
  Double_t dPadX          =                    0.5*                   AliHMPIDParam::PadSizeX()  ;//pad X half size parametrised
  Double_t dPadY          =                    0.5*                   AliHMPIDParam::PadSizeY()  ;//pad Y half size parametrised
  Double_t dPadZ          =1.0 *mm            *0.5                                              ;//CsI film thickness 
  Double_t dGapX          =                    0.5*                          iNpadsX*2*dPadX    ;//gap x half size n. pads x * pad size
  Double_t dGapY          =                    0.5*                          iNpadsY*2*dPadY    ;//gap y half size
  Double_t dGapZ          =                    0.5*                    (2*dPadZ+wCathPc+wCathR) ;//gap half thickness
  Double_t dMirX          =                    0.5*     2*dGapX/TMath::Cos(45*TMath::DegToRad());//Ag mirror x half size defined by gap size and angle 45 degrees
  Double_t dMirY          =                    0.5*                           2*dGapY           ;//Ag mirror y half size defined by gap size
  Double_t dMirZ          =1.0 *mm            *0.5                                              ;//Ag mirror z half size defined by USER
  Double_t wBoxWall       =2.0 *mm                                                              ;//Al box walls thickness defined by USER
  Double_t dBoxX          =                    0.5*                         (2*dGapX+2*cm)      ;//Al box x half size defined by gap size 2 cm for tolerance
  Double_t dBoxY          =                    0.5*                         (2*dGapY+2*cm)      ;//Al box y half size defined by gap size 2 cm for tolerance 
  Double_t dBoxZ          =1.8*m              *0.5                                              ;//Al box z half size defined by USER
  Double_t dWinX          =                                                  (dBoxX-wBoxWall)   ;//SiO2 window x half size defined by box size
  Double_t dWinY          =                                                  (dBoxY-wBoxWall)   ;//SiO2 window y half size defined by box size
  Double_t dWinZ          =1.0*cm             *0.5                                              ;//SiO2 window z half size defined by USER
  Double_t dCF4X          =                                                  (dBoxX-wBoxWall)   ;//CF4 radiator x half size defined by box size
  Double_t dCF4Y          =                                                  (dBoxY-wBoxWall)   ;//CF4 radiator y half size defined by box size
  Double_t dCF4Z          =                                                  0.4*dBoxZ          ;//CF4 radiator z half size defined by box size or by USER
  Double_t dC4F10X        =                                                  (dBoxX-wBoxWall)   ;//C4F10 radiator x half size defined by box size
  Double_t dC4F10Y        =                                                  (dBoxY-wBoxWall)   ;//C4F10 radiator y half size defined by box size
  Double_t dC4F10Z        =                                      (dBoxZ-dWinZ-dCF4Z-wBoxWall)   ;//C4F10 radiator z half size defined by box, CF4 and window sizes
  
  Int_t copy;  
  Double_t rmin,rmax,dx,dy,dz;
  
//make VHMPID type 2 volume  (2 radiators)
  TGeoVolume *pBox=gGeoManager->MakeBox("Vbox",gGeoManager->GetMedium("HMPID_Al"),dx=dBoxX,dy=dBoxY,dz=dBoxZ);
  
  TGeoRotation *pRot=new TGeoRotation("VboxRot"); pRot->RotateX(90);//normal position
  gGeoManager->GetVolume("ALIC")->AddNode(pBox,copy=1,new TGeoCombiTrans(0*m,-5.2*m,2.5*m,pRot));
//position C4F10 radiator to Al box   
  TGeoVolume *pC4F10=gGeoManager->MakeBox("Vc4f",gGeoManager->GetMedium("HMPID_C4F10"),dx=dC4F10X,dy=dC4F10Y,dz=dC4F10Z);      
  pBox->AddNode(pC4F10,copy=1,new TGeoTranslation(0*cm,0*cm,-dBoxZ+wBoxWall+dC4F10Z)); 
//position quartz window  to Al box  
  TGeoVolume *pWindow=gGeoManager->MakeBox( "Vwin",gGeoManager->GetMedium("HMPID_SiO2"),dx=dWinX,dy=dWinY,dz=dWinZ); 
  pBox->AddNode(pWindow,copy=1,new TGeoTranslation(0*cm,0*cm,-dBoxZ+wBoxWall+2*dCF4Z+dWinZ)); 
//position CF4 radiator to Al box   
  TGeoVolume *pCF4=gGeoManager->MakeBox( "Vcf4",gGeoManager->GetMedium("HMPID_CF4"),dx=dCF4X,dy=dCF4Y,dz=dCF4Z);                   
  pBox->AddNode(pCF4,copy=1,new TGeoTranslation(0*cm,0*cm,dBoxZ-wBoxWall-dCF4Z)); 
//position mirror to CF4 radiator   
  TGeoVolume *pMirror=gGeoManager->MakeBox( "Vmir",gGeoManager->GetMedium("HMPID_Ag"),dx=dMirX,dy=dMirY,dz=dMirZ);     
  TGeoRotation *pMirrorRot=new TGeoRotation("VmirRot"); pMirrorRot->RotateY(45);   
  pCF4->AddNode(pMirror,copy=1,new TGeoCombiTrans(0*cm,0*cm,dCF4Z-1*cm-dGapX,pMirrorRot)); 
//position gap to  CF4 radiator   
  TGeoVolume   *pGap     =gGeoManager->MakeBox ("Vgap"     ,gGeoManager->GetMedium("HMPID_CF4"),dx=dGapX , dy=dGapY , dz=dGapZ );  
  TGeoRotation *pMwpcRot=new TGeoRotation("VmpcRot"); pMwpcRot->RotateY(90);  
  pCF4->AddNode(pGap,copy=1,new TGeoCombiTrans(-dBoxX+1*cm,0*cm,dCF4Z-1*cm-dGapX,pMwpcRot)); 
//divide gap into 80x48 cells
  Int_t axis,ndiv; Double_t start,step;
  TGeoVolume *pGapCol =pGap      ->Divide("Vcol",axis=1,ndiv=iNpadsX,start=0,step=0);//divide VhGap along X by NpadsX columns
  TGeoVolume *pGapCell=pGapCol   ->Divide("Vcel",axis=2,ndiv=iNpadsY,start=0,step=0);//divide VhGapCol along Y by NpadsY cells
//position pad to gap cell
  TGeoVolume *pPad=gGeoManager->MakeBox ("Vpad",gGeoManager->GetMedium("HMPID_CsI"),dx=dPadX,dy=dPadY,dz=dPadZ);      
  pGapCell->AddNode(pPad,copy=1,new TGeoTranslation(0,0,-dGapZ+dPadZ)); 
//define wire rotation common for both anod and cathode wires
  TGeoRotation *pWireRot=new TGeoRotation("VwireRot"); pWireRot->RotateY(90); //rotate wires around Y to be along X (initially along Z)
//position 2 anod  wires to gap cell
  TGeoVolume *pAnodWire =gGeoManager->MakeTube("Vano",gGeoManager->GetMedium("HMPID_W")  ,rmin=0   , rmax=wAnodR   , dz=dPadX );  
  pGapCell->AddNode(pAnodWire,copy=1,new TGeoCombiTrans (0, -dPadY+wAnodShift             , -dGapZ+wAnodPc+2*dPadZ , pWireRot)); 
  pGapCell->AddNode(pAnodWire,copy=2,new TGeoCombiTrans (0,  dPadY-wAnodShift             , -dGapZ+wAnodPc+2*dPadZ , pWireRot)); 
//position 4 cathode wires to gap cell  
  TGeoVolume *pCathWire =gGeoManager->MakeTube("Vcat",gGeoManager->GetMedium("HMPID_Cu") ,rmin=0   , rmax=wCathR , dz=dPadX );
  pGapCell->AddNode(pCathWire,copy=1,new TGeoCombiTrans (0, -dPadY+wCathShift            , -dGapZ+wCathPc+2*dPadZ , pWireRot)); 
  pGapCell->AddNode(pCathWire,copy=2,new TGeoCombiTrans (0, -dPadY+wCathShift+wCathPitch , -dGapZ+wCathPc+2*dPadZ , pWireRot)); 
  pGapCell->AddNode(pCathWire,copy=3,new TGeoCombiTrans (0,  dPadY-wCathShift-wCathPitch , -dGapZ+wCathPc+2*dPadZ , pWireRot)); 
  pGapCell->AddNode(pCathWire,copy=4,new TGeoCombiTrans (0,  dPadY-wCathShift            , -dGapZ+wCathPc+2*dPadZ , pWireRot)); 
}//Vhmpid()
//__________________________________________________________________________________________________
void Axis()
{
// Draw axises  on top of geometry    
  Double_t X[6]={0,0,0,300,0,0};  Double_t Y[6]={0,0,0,0,300,0};  Double_t Z[6]={0,0,0,0,0,300};  
  TPolyLine3D *pXaxis=new TPolyLine3D(2,X);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,Y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,Z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();  
}
//__________________________________________________________________________________________________
