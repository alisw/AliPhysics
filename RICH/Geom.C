TGeoManager  *g=0;
AliRICHParam *p=0;
Int_t copy; //volume copy number 
Double_t dx,dy,dz,r1,r2;//tmp vars for volume dimentions
Double_t cm=1,m=100*cm,mm=0.1*cm,mkm=0.001*cm;//length units  

void Geom()
{
  AliRICHParam::SetRadioSrc(kTRUE);
  p=new AliRICHParam;
  
  //gSystem->Load("libGeom.so");  
  g=new TGeoManager("GEL","Aerogel test configuration");
  Materials();  
  TGeoVolume *pMother=g->MakeBox("Mother",g->GetMedium("Air"),dx=30*m/2,dy=30*m/2,dz=30*m/2); //arbitrary values  
  g->SetTopVolume(g->GetVolume("Mother"));
    
  Rich(pMother);
  
  Colors();  
  g->CloseGeometry();
  
  g->SetVisOption(0); g->SetVisLevel(5); 
  
  g->GetMasterVolume()->Draw();
  Axis();
//  gPad->GetView()->SetView(3,94,-70,0);
  new TBrowser;
}
//__________________________________________________________________________________________________
void Aerogel()
{
  
  TGeoVolume *pGel     =g->MakeBox("RGEL" ,g->GetMedium("Air"),dx=10*cm/2,dy= 10*cm/2 ,dz=10*cm/2);//10x10x10 cm aerogel cubic
  for(int i=1;i<=7;i++)//put 7 cubics
    g->GetVolume("Mother")->AddNode(pGel,copy=i,new TGeoCombiTrans(p->C(i)->Center().X(), p->C(i)->Center().Y(), p->C(i)->Center().Z(),
                                                new TGeoRotation(Form("GelRot%i",i),p->C(i)->ThetaXd(),p->C(i)->PhiXd(),
                                                                                    p->C(i)->ThetaYd(),p->C(i)->PhiYd(),
                                                                                    p->C(i)->ThetaZd(),p->C(i)->PhiZd())));

}
//__________________________________________________________________________________________________
void Materials()
{
// Defines all the materials and tracking media
  Double_t a=0,z=0,den=0,radlen=0,intlen=0;//tmp vars for material parameters
  new TGeoMaterial("Air"      ,a=26.98 ,z=13,den=0.4224                             );        new TGeoMedium("Air"      ,1,g->GetMaterial("Air")); 
  new TGeoMaterial("CH4"      ,a=26.98 ,z=13,den=0.4224                             ); pCH4  =new TGeoMedium("CH4"      ,2,g->GetMaterial("CH4"));
  new TGeoMaterial("CsI"      ,a=26.98 ,z=13,den=2.7 ,radlen=24.01*cm,intlen=70.6*cm); pCsI  =new TGeoMedium("CsI"      ,3,g->GetMaterial("CsI"));
  new TGeoMaterial("Al"       ,a=26.98 ,z=13,den=2.7 ,radlen=24.01*cm,intlen=70.6*cm); pAl   =new TGeoMedium("Al"       ,4,g->GetMaterial("Al")); 
  new TGeoMaterial("W"        ,a=183.84,z=27,den=19.3,radlen= 9.59*cm,intlen=0.35*cm); pW    =new TGeoMedium("W"        ,5,g->GetMaterial("W"));
  new TGeoMaterial("Cu"       ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); pCu   =new TGeoMedium("Cu"       ,6,g->GetMaterial("Cu"));
  new TGeoMaterial("Rohacell" ,a=12.01 ,z=6 ,den=0.1 ,radlen=18.8,intlen=0);           pRoha =new TGeoMedium("Rohacell" ,7,g->GetMaterial("Rohacell"));
  new TGeoMaterial("SiO2"     ,a=0     ,z=0 ,den=0);                                   pSiO2 =new TGeoMedium("SiO2"     ,8,g->GetMaterial("SiO2"));
  new TGeoMaterial("C6F14"    ,a=0     ,z=0 ,den=0);                                   pC6F14=new TGeoMedium("C6F14"    ,9,g->GetMaterial("C6F14"));
//Medium for Sr90 source  
  new TGeoMaterial("Perpex"  ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); pPerpex=new TGeoMedium("Perpex"  ,10,g->GetMaterial("Perpex"));
  new TGeoMaterial("Steel"   ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); pSteel =new TGeoMedium("Steel"   ,11,g->GetMaterial("Steel"));
  new TGeoMaterial("Mylar"   ,a=55.845,z=26,den=7.87,radlen=13.84*cm,intlen=82.8*cm); pMylar =new TGeoMedium("Mylar"   ,12,g->GetMaterial("Mylar"));
  new TGeoMaterial("Sr90"    ,a=87.62 ,z=38,den=7.87,radlen=13.84*cm,intlen=82.8*cm); pSr    =new TGeoMedium("Sr90"    ,13,g->GetMaterial("Sr90"));
}//Materials()
//__________________________________________________________________________________________________
void Rich(TGeoVolume *pTop)
{
//Rich  
  TGeoVolume *pRich=g->MakeBox( "Rich"     ,g->GetMedium("CH4"),dx=(6*mm+1681*mm+6*mm)/2,  //main RICH volume
                                                                dy=(6*mm+1466*mm+6*mm)/2,
                                                                dz=(80*mm+40*mm)*2/2);     //x,y taken from 2033P1  z from p84 TDR  
  for(int i=1;i<=p->Nchambers();i++)//put 7 chambers
    pTop->AddNode(pRich,copy=i,new TGeoCombiTrans(p->C(i)->Center().X(), p->C(i)->Center().Y(), p->C(i)->Center().Z(),
                                                 new TGeoRotation(Form("RichRot%i",i),p->C(i)->ThetaXd(),p->C(i)->PhiXd(),
                                                                                      p->C(i)->ThetaYd(),p->C(i)->PhiYd(),
                                                                                      p->C(i)->ThetaZd(),p->C(i)->PhiZd()))); //Rich to Mother    
  
                 PadPanelFrame(pRich); //photocathode
                 Gap(pRich);           //anod wires
                 Frame3(pRich);        //cathode wires
                 Frame4(pRich);        //collection wires  
  if(!p->IsRadioSrc())    Radiators(pRich);
  else                    Sr90(pRich);
                 Sandbox(pRich);
}//Rich()
//__________________________________________________________________________________________________
void PadPanelFrame(TGeoVolume *pTop)
{
//Pad Panel frame  
  TGeoVolume *pPpf     =g->MakeBox("PPF"      ,g->GetMedium("Al")  ,dx=648*mm/2,dy=  411*mm/2 ,dz=40  *mm/2);//PPF 2001P2 inner size of the slab by 1mm more
  TGeoVolume *pPpfLarge=g->MakeBox("PPFlarge" ,g->GetMedium("Air") ,dx=181*mm/2,dy=89.25*mm/2 ,dz=38.3*mm/2);//large whole
  TGeoVolume *pPpfSmall=g->MakeBox("PPFsmall" ,g->GetMedium("Air") ,dx=114*mm/2,dy=89.25*mm/2 ,dz=38.3*mm/2);//small whole
  TGeoVolume *pPc      =g->MakeBox("PC"       ,g->GetMedium("CsI") ,dx=644*mm/2,dy=  407*mm/2 ,dz= 1.7*mm/2);//by 0.2 mm more then actual size (PCB 2006P1)
  
  pTop->AddNode(pPpf,copy=1,new TGeoTranslation(-335*mm,-433*mm,8*cm+20*mm));//F1 2040P1 z p.84 TDR
  pTop->AddNode(pPpf,copy=2,new TGeoTranslation(+335*mm,-433*mm,8*cm+20*mm));
  pTop->AddNode(pPpf,copy=3,new TGeoTranslation(-335*mm,   0*mm,8*cm+20*mm));
  pTop->AddNode(pPpf,copy=4,new TGeoTranslation(+335*mm,   0*mm,8*cm+20*mm));
  pTop->AddNode(pPpf,copy=5,new TGeoTranslation(-335*mm,+433*mm,8*cm+20*mm));
  pTop->AddNode(pPpf,copy=6,new TGeoTranslation(+335*mm,+433*mm,8*cm+20*mm));
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
}//PadPanelFrame()
//__________________________________________________________________________________________________
void Gap(TGeoVolume *pTop)
{
//Gap - anod wires
  TGeoVolume *pGap =g->MakeBox ("Gap"  ,g->GetMedium("CH4") ,dx=648*mm/2,dy=  411*mm/2 ,dz=4.45*mm/2);//xy as PPF 2001P2 z WP 2099P1
  TGeoVolume *pAnod=g->MakeTube("Anod" ,g->GetMedium("W")   ,r1=  0*mm  ,r2=  20*mkm/2 ,dz=648*mm/2);  //WP 2099P1 z = gap x PPF 2001P2
  TGeoRotation *pAnodRot=new TGeoRotation("AnodRot",90,90,0);
  
  pTop->AddNode(pGap,copy=1,new TGeoTranslation(-335*mm,-433*mm,8*cm-2.225*mm)); //F1 2040P1 z WP 2099P1
  pTop->AddNode(pGap,copy=2,new TGeoTranslation(+335*mm,-433*mm,8*cm-2.225*mm)); 
  pTop->AddNode(pGap,copy=3,new TGeoTranslation(-335*mm,   0*mm,8*cm-2.225*mm)); 
  pTop->AddNode(pGap,copy=4,new TGeoTranslation(+335*mm,   0*mm,8*cm-2.225*mm)); 
  pTop->AddNode(pGap,copy=5,new TGeoTranslation(-335*mm,+433*mm,8*cm-2.225*mm)); 
  pTop->AddNode(pGap,copy=6,new TGeoTranslation(+335*mm,+433*mm,8*cm-2.225*mm)); 
  for(int i=1;i<=96;i++)
    pGap->AddNode(pAnod,copy=i,new TGeoCombiTrans( 0*mm, -411/2*mm+i*4*mm, 0.185*mm,pAnodRot)); //WP 2099P1  
}//Gap()
//__________________________________________________________________________________________________
void Frame3(TGeoVolume *pTop)
{//Frame 3- cathode wires      
  TGeoVolume *pCath=g->MakeTube("Cath"  ,g->GetMedium("Cu")  ,r1=0  ,r2=100*mkm/2,dz=1323*mm/2);//r WP 2099P1 z F3 2041P1       
  TGeoRotation *pCathRot=new TGeoRotation("CathRot",90,90,0);
  for(int i=1;i<=618;i++)
    pTop->AddNode(pCath,copy=i,new TGeoCombiTrans( 0*mm, -649.5*mm+i*2.1*mm, 75*mm,pCathRot)); //WP 2099P1    
}//Frame3()
//__________________________________________________________________________________________________
void Frame4(TGeoVolume *pTop)
{
//Frame 4- collection wires
  TGeoVolume *pF4  =g->MakeBox( "F4"       ,g->GetMedium("CH4")   ,dx=1407*mm/2 ,dy=1366*mm/2  ,dz=  15*mm/2);//F4 2043P1 
  TGeoVolume *pF4al=g->MakeBox( "F4al"     ,g->GetMedium("Al")    ,dx=1407*mm/2 ,dy=1366*mm/2  ,dz=  10*mm/2); 
  TGeoVolume *pF4in=g->MakeBox( "F4in"     ,g->GetMedium("CH4")   ,dx=1323*mm/2 ,dy=1296*mm/2  ,dz=  10*mm/2); 
  TGeoVolume *pColl=g->MakeTube("Coll"     ,g->GetMedium("Cu")    ,r1=   0*mm   ,r2=100*mkm/2  ,dz=1323*mm/2);
  TGeoRotation *pCollRot=new TGeoRotation("CollRot",90,90,0);
  
  pTop->AddNode(pF4      ,copy=1,new TGeoTranslation(   0*mm,0*mm,   9*mm)); //F4 to Rich p.84 TDR
    pF4  ->AddNode(pF4al    ,copy=1,new TGeoTranslation(   0*mm,0*mm, 2.5*mm)); //F4 al to F4 2043P1 
    pF4al->AddNode(pF4in    ,copy=1,new TGeoTranslation(   0*mm,0*mm,   0*mm)); //F4 whole F4 al 2043P1   
    for(int i=1;i<=322;i++)
      pF4->AddNode(pColl,copy=i,new TGeoCombiTrans( 0*mm, -1296/2*mm+i*4*mm, -5*mm,pCollRot)); //F4 2043P1
}//void Frame4()
//__________________________________________________________________________________________________
void Radiators(TGeoVolume *pTop)
{
  TGeoVolume *pRad      =g->MakeBox( "Rad"      ,g->GetMedium("C6F14")    ,dx=1330*mm/2 ,dy= 413*mm/2  ,dz=  24*mm/2); // Rad 2011P1
  TGeoVolume *pRadFront =g->MakeBox( "RadFront" ,g->GetMedium("Neoceram") ,dx=1330*mm/2 ,dy= 413*mm/2  ,dz=   4*mm/2); 
  TGeoVolume *pRadWin   =g->MakeBox( "RadWin"   ,g->GetMedium("SiO2")     ,dx=1330*mm/2 ,dy= 413*mm/2  ,dz=   5*mm/2); 
  TGeoVolume *pRadLong  =g->MakeBox( "RadLong"  ,g->GetMedium("Neoceram") ,dx=1330*mm/2 ,dy=   5*mm/2  ,dz=  15*mm/2);    
  TGeoVolume *pRadShort =g->MakeBox( "RadShort" ,g->GetMedium("Neoceram") ,dx=  10*mm/2 ,dy= 403*mm/2  ,dz=  15*mm/2);    
  TGeoVolume *pRadSpacer=g->MakeTube("RadSpacer",g->GetMedium("SiO2")     ,r1= 0        ,r2=10*mm/2  ,dz=  15*mm/2);         
    
  pTop->AddNode(pRad      ,copy=1,new TGeoTranslation(   0*mm,-434*mm,   -12*mm)); 
  pTop->AddNode(pRad      ,copy=2,new TGeoTranslation(   0*mm,   0*mm,   -12*mm)); 
  pTop->AddNode(pRad      ,copy=3,new TGeoTranslation(   0*mm,+434*mm,   -12*mm)); 
    
  pRad ->AddNode(pRadFront ,copy=1,new TGeoTranslation(   0*mm,   0*mm, -10.0*mm));
  pRad ->AddNode(pRadWin   ,copy=1,new TGeoTranslation(   0*mm,   0*mm,   9.5*mm));
  pRad ->AddNode(pRadLong  ,copy=1,new TGeoTranslation(   0*mm,-204*mm,  -0.5*mm));
  pRad ->AddNode(pRadLong  ,copy=2,new TGeoTranslation(   0*mm,+204*mm,  -0.5*mm));
  pRad ->AddNode(pRadShort ,copy=1,new TGeoTranslation(-660*mm,   0*mm,  -0.5*mm));
  pRad ->AddNode(pRadShort ,copy=2,new TGeoTranslation(+660*mm,   0*mm,  -0.5*mm));
  for(int i=0;i<3;i++)
    for(int j=0;j<10;j++)
      pRad->AddNode(pRadSpacer,copy=10*i+j,new TGeoTranslation(-1330*mm/2+116*mm+j*122*mm,(i-1)*105*mm,-0.5*mm));
}//Radiators()  
//__________________________________________________________________________________________________
void Sr90(TGeoVolume *pTop)
{
    pSrc               =g->MakeTube("Src"              ,g->GetMedium("CH4")   , 0 , 70*mm/2 ,  30*mm/2);       //top container
      pAlGlass         =g->MakeTube("SrcAlGlass"       ,g->GetMedium("Al")    , 0 , 38*mm/2 ,21.8*mm/2);       //Al glass wall        
        pPerpexPlug    =g->MakeTube("SrcPerpex"        ,g->GetMedium("Perpex"), 0 , 34*mm/2 ,  20*mm/2);       //Perpex plug         
          pScrewCentral=g->MakeTube("SrcScrewCentral"  ,g->GetMedium("Steel") , 0 ,  5*mm/2 ,  15*mm/2);       //Steel screw in the center        
          pScrewSr90   =g->MakeTube("SrcScrewSr90"     ,g->GetMedium("Steel") , 0 ,  2*mm/2 ,  10*mm/2);       //Steel screw to support Sr90 
            pSr90      =g->MakeTube("SrcSr90"          ,g->GetMedium("Sr90")  , 0 ,  1*mm/2 ,   1*mm/2);       //Sr90 source
          pHolePerpex  =g->MakeTube("SrcHolePerpex"    ,g->GetMedium("Air")   , 0 ,  4*mm/2 ,  10*mm/2);       //Air hole in perpex plug      
        pHoleAl        =g->MakeTube("SrcHoleAl"        ,g->GetMedium("Air")   , 0 ,  5*mm/2 , 1.8*mm/2);       //Air hole in Al glass bottom
    pMylarFoil         =g->MakeTube("SrcMylarFoil"     ,g->GetMedium("Mylar") , 0 , 30*mm/2 , 50*mkm/2);       //Mylar foil                
                
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
void Sandbox(TGeoVolume *pTop)
{//Sandbox  
  TGeoVolume *pSandBox  =g->MakeBox( "SandBox"  ,g->GetMedium("Air")  ,dx=1419*mm/2 ,dy=1378*mm/2   ,dz=50.5*mm/2);  //2072P1   
  TGeoVolume *pSandCover=g->MakeBox( "SandCover",g->GetMedium("Al")   ,dx=1419*mm/2 ,dy=1378*mm/2   ,dz= 0.5*mm/2);  
  TGeoVolume *pSandComb =g->MakeBox( "SandComb" ,g->GetMedium("Roha") ,dx=1359*mm/2 ,dy=1318*mm/2   ,dz=49.5*mm/2);  
  
  pTop->AddNode(pSandBox,copy=1,new TGeoTranslation(  0*mm,0*mm, -73.75*mm)); //p.84 TDR
    pSandBox->AddNode(pSandComb  ,copy=1,new TGeoTranslation(  0*mm,0*mm,      0*mm)); //2072P1
    pSandBox->AddNode(pSandCover ,copy=1,new TGeoTranslation(  0*mm,0*mm,    +25*mm)); 
    pSandBox->AddNode(pSandCover ,copy=2,new TGeoTranslation(  0*mm,0*mm,    -25*mm)); 
}//Sandbox()
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
void Colors()
{
//Set volume colors  
  TGeoVolume *pVol=0;
  pVol=g->GetVolume("Pc")       ; if(pVol) pVol->SetLineColor(kGreen);
  pVol=g->GetVolume("PPFlarge") ; if(pVol) pVol->SetLineColor(kYellow);
  pVol=g->GetVolume("PPFsmall") ; if(pVol) pVol->SetLineColor(kYellow);
  pVol=g->GetVolume("RadFront") ; if(pVol) pVol->SetLineColor(kRed);
  pVol=g->GetVolume("RadLong")  ; if(pVol) pVol->SetLineColor(46);
  pVol=g->GetVolume("RadShort") ; if(pVol) pVol->SetLineColor(kMagenta);
  pVol=g->GetVolume("RadWin")   ; if(pVol) pVol->SetLineColor(kBlue);
  pVol=g->GetVolume("RadSpacer"); if(pVol) pVol->SetLineColor(kYellow);
  pVol=g->GetVolume("SrcSr90")       ; if(pVol) pVol->SetLineColor(kRed);
  pVol=g->GetVolume("SrcScrewSr90")  ; if(pVol) pVol->SetLineColor(kGreen);
  pVol=g->GetVolume("SrcHolePerpex") ; if(pVol) pVol->SetLineColor(26);
  pVol=g->GetVolume("SrcHoleAl")     ; if(pVol) pVol->SetLineColor(27);
}
