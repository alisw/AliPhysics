void ss()
{
  if(rl->LoadSDigits()) return;
  rl->TreeS()->GetEntry(0);
  r->SDigits()->Print();
  Info("sd","totally %i",r->SDigits()->GetEntries());
  rl->UnloadSDigits();
}

void sd()
{
  if(rl->LoadDigits()) return;
  rl->TreeD()->GetEntry(0);
  for(int i=1;i<=7;i++) r->Digits(i)->Print();
  rl->UnloadDigits();
}

void sc()
{
  if(rl->LoadRecPoints()) return;
  rl->TreeR()->GetEntry(0);
  for(int i=1;i<=7;i++) r->Clusters(i)->Print();
  rl->UnloadRecPoints();
}

Double_t r2d = TMath::RadToDeg();
Double_t d2r = TMath::DegToRad();

void DisplFast()
{
  AliRICHDisplFast *d = new AliRICHDisplFast();

  for (int nev=0; nev< a->GetEventsPerRun(); nev++) {    // Event Loop
    al->GetEvent(nev);
    cout <<endl<< "Processing event:" <<nev<<endl;
    d->Display();
  } // event loop  
}  


void Digits2Recos()
{
  AliRICHRecon *detect = new AliRICHRecon("RICH patrec algorithm","Reconstruction");
    

  for (int nev=0; nev< a->GetEventsPerRun(); nev++) {    // Event Loop
    al->GetEvent(nev);
    cout <<endl<< "Processing event:" <<nev<<endl;
    detect->StartProcessEvent();
  } // event loop  
  delete detect;
}  



void MainTranck()
{
  TStopwatch sw;TDatime time;
  OLD_S_SD(); SD_D(); D_C();
  cout<<"\nInfo in <hits->digits>: Start time: ";time.Print();
  cout<<"Info in <hits->digits>: Stop  time: ";time.Set();  time.Print();
  cout<<"Info in <hits->digits>: Time  used: ";sw.Print();
}

void D_C()
{
  AliRICHClusterFinder *z=new AliRICHClusterFinder(r);
  z->Exec();  
}
//__________________________________________________________________________________________________
void OLD_S_SD()
{
  Info("OLD_S_SD","Start.");
  
  rl->LoadHits(); 
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
    
    rl->MakeTree("S");  r->MakeBranch("S");
    r->ResetSDigits();  r->ResetSpecialsOld();

    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t i=0;i<r->Specials()->GetEntries();i++){//specials loop          
        Int_t padx= ((AliRICHSDigit*)r->Specials()->At(i))->PadX();
        Int_t pady= ((AliRICHSDigit*)r->Specials()->At(i))->PadY();
        Double_t q=  ((AliRICHSDigit*)r->Specials()->At(i))->QPad();
        Int_t hitN= ((AliRICHSDigit*)r->Specials()->At(i))->HitNumber()-1;//!!! important -1
        Int_t chamber=((AliRICHhit*)r->Hits()->At(hitN))->C();
        Int_t track=((AliRICHhit*)r->Hits()->At(hitN))->GetTrack();
        r->AddSDigit(chamber,padx+r->Param()->NpadsX()/2,pady+r->Param()->NpadsY()/2,q,track);
      }//specials loop
    }//prims loop
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop  
    rl->UnloadHits();     rl->UnloadSDigits();  
  Info("OLD_S_SD","Stop.");    
}//Specials2Sdigits()
//__________________________________________________________________________________________________
void H_SD()
{
  Info("H_SD","Start.");
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
  
    if(!rl->TreeH()) rl->LoadHits();//from
    if(!rl->TreeS()) rl->MakeTree("S");    r->MakeBranch("S");//to
      
    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<3;iHitN++){//hits loop  ???
        AliRICHhit *pHit=r->Hits()->At(iHitN);        
        TVector3 globX3(pHit->X(),pHit->Y(),pHit->Z());        
        TVector3 locX3=r->C(pHit->C())->Glob2Loc(globX3);
        
        Int_t sector;
        Int_t iTotQdc=r->Param()->Loc2TotQdc(locX3,pHit->Eloss(),pHit->Pid(),sector);
        
        Int_t iPadXmin,iPadXmax,iPadYmin,iPadYmax;
        r->Param()->Loc2Area(locX3,iPadXmin,iPadYmin,iPadXmax,iPadYmax);
        cout<<"left-down=("<<iPadXmin<<","<<iPadYmin<<") right-up=("<<iPadXmax<<','<<iPadYmax<<')'<<endl;
        for(Int_t iPadY=iPadYmin;iPadY<=iPadYmax;iPadY++)
          for(Int_t iPadX=iPadXmin;iPadX<=iPadXmax;iPadX++){
            Double_t padQdc=iTotQdc*r->Param()->Loc2PadFrac(locX3,iPadX,iPadY);
            if(padQdc>0.1)r->AddSDigit(pHit->C(),iPadX,iPadY,padQdc,pHit->GetTrack());
          }            
      }//hits loop
    }//prims loop
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop
  
  rl->UnloadHits();
  rl->UnloadSDigits();  
  Info("H_SD","Stop.");  
}//Hits2Sdigits()
//__________________________________________________________________________________________________
void SD_D()
{
  Info("SD_D","Start.");  

  rl->LoadSDigits();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
    
    rl->MakeTree("D");r->MakeBranch("D"); //create TreeD with RICH branches 
    r->ResetSDigits();r->ResetDigits();//reset lists of sdigits and digits
    rl->TreeS()->GetEntry(0);  
    r->SDigits()->Sort();
  
    Int_t kBad=-101;
    Int_t chamber,x,y,tr[3],id;
    Double_t q=kBad;
    chamber=x=y=tr[0]=tr[1]=tr[2]=id=kBad;
    Int_t iNdigitsPerPad=kBad;//how many sdigits for a given pad
        
    for(Int_t i=0;i<r->SDigits()->GetEntries();i++){//sdigits loop (sorted)
      AliRICHdigit *pSdig=(AliRICHdigit*)r->SDigits()->At(i);
      if(pSdig->Id()==id){//still the same pad
        iNdigitsPerPad++;
        q+=pSdig->Q();
        if(iNdigitsPerPad<=3)
          tr[iNdigitsPerPad-1]=pSdig->T(0);
        else
          Info("","More then 3 sdigits for the given pad");
      }else{//new pad, add the pevious one
        if(id!=kBad) r->AddDigit(chamber,x,y,q,tr[0],tr[1],tr[2]);//ch-xpad-ypad-qdc-tr1-2-3
        chamber=pSdig->C();x=pSdig->X();y=pSdig->Y();q=pSdig->Q();tr[0]=pSdig->T(0);id=pSdig->Id();
        iNdigitsPerPad=1;tr[1]=tr[2]=kBad;
      }
    }//sdigits loop (sorted)
  
    if(r->SDigits()->GetEntries())r->AddDigit(chamber,x,y,q,tr[0],tr[1],tr[2]);//add the last digit
        
    rl->TreeD()->Fill();  
    rl->WriteDigits("OVERWRITE");
  }//events loop
  rl->UnloadSDigits();     rl->UnloadDigits();  
  r->ResetSDigits();r->ResetDigits();//reset lists of sdigits and digits
  Info("SD_D","Stop.");  
}



void OLD_SD_D()
{
  Info("SD_DOLD","Start.");  

  rl->LoadSDigits();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
    
    rl->MakeTree("D");r->MakeBranch("D"); //create TreeD with RICH branches 
    r->ResetSDigits();r->ResetDigitsOld();//reset lists of sdigits and digits
    rl->TreeS()->GetEntry(0);  
    r->SDigits()->Sort();
  
    Int_t kBad=-101;
    
    Int_t tr[3],q[3],dig[5]; for(Int_t i=0;i<3;i++) tr[i]=q[i]=kBad;    for(Int_t i=0;i<5;i++) dig[i]=kBad;        
    Int_t chamber=kBad,id=kBad,iNdigitsPerPad=kBad;//how many sdigits for a given pad
        
    for(Int_t i=0;i<r->SDigits()->GetEntries();i++){//sdigits loop (sorted)
      AliRICHdigit *pSdig=(AliRICHdigit*)r->SDigits()->At(i);
      if(pSdig->Id()==id){//still the same pad
        iNdigitsPerPad++;
        dig[2]+=pSdig->Q();//sum up qdc
        if(iNdigitsPerPad<=3)
          tr[iNdigitsPerPad-1]=pSdig->T(0);
        else
          Info("","More then 3 sdigits for the given pad");
      }else{//new pad, add the pevious one
        if(id!=kBad) r->AddDigitOld(chamber,tr,q,dig);
        chamber=pSdig->C();dig[0]=pSdig->X();dig[1]=pSdig->Y();dig[2]=pSdig->Qdc();tr[0]=pSdig->T(0);id=pSdig->Id();
        iNdigitsPerPad=1;tr[1]=tr[2]=kBad;
      }
    }//sdigits loop (sorted)
    r->AddDigitOld(chamber,tr,q,dig);//add the last digit
        
    rl->TreeD()->Fill();  
    rl->WriteDigits("OVERWRITE");
  }//events loop
  rl->UnloadSDigits();     rl->UnloadDigits();  
  r->ResetSDigits();r->ResetDigitsOld();//reset lists of sdigits and digits
  Info("SD_DOLD","Stop.");  
}


void Show3()
{  
  cout<<endl;
  al->LoadHeader();  al->LoadKinematics();
  
  rl->LoadHits();  Bool_t isSdigits=!rl->LoadSDigits();  Bool_t isDigits=!rl->LoadDigits();//loaders
  
  cout<<endl;  cout<<endl;  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    Int_t iNparticles=a->GetEvent(iEventN);
    Int_t iNprims=rl->TreeH()->GetEntries();
    
    Int_t iTotalHits=0,iTotalCerenkovs=0,iTotalSpecials=0;
    for(Int_t iPrimN=0;iPrimN<iNprims;iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);      
      iTotalHits+=r->Hits()->GetEntries();
      iTotalCerenkovs+=r->Cerenkovs()->GetEntries();
      iTotalSpecials+=r->Specials()->GetEntries();
      TParticle *pPrim=al->Stack()->Particle(iPrimN);
      Info("Show","Evt %4i prim %4i has %4i hits %5i cerenkovs and %5i specials from %s (,%7.2f,%7.2f)",
                           iEventN,
                                    iPrimN,
                                             r->Hits()->GetEntries(),
                                                      r->Cerenkovs()->GetEntries(),
                                                                        r->Specials()->GetEntries(),
                                                                                         pPrim->GetName(),
                                                                                 pPrim->Theta()*r2d,pPrim->Phi()*r2d);
    }//prims loop
    Info("Show-HITS","Evt %i total:  %i particles %i primaries %i hits %i cerenkovs %i specials",
                        iEventN,   iNparticles, iNprims,     iTotalHits,iTotalCerenkovs,iTotalSpecials);
    if(isSdigits){
      rl->TreeS()->GetEntry(0);
      Info("Show-SDIGITS","Evt %i contains %5i sdigits",iEventN,r->SDigits()->GetEntries());
    }
    if(isDigits){
      rl->TreeD()->GetEntry(0);
      for(int i=1;i<=7;i++)
        Info("Show-DIGITS","Evt %i chamber %i contains %5i NEW digits and %5i OLD",
                                 iEventN,   i,           r->Digits(i)->GetEntries(),r->DigitsOld(i)->GetEntries());
    }
    cout<<endl;
  }//events loop
  rl->UnloadHits();  if(isSdigits) rl->UnloadSDigits(); if(isDigits) rl->UnloadDigits();
  al->UnloadHeader();
  al->UnloadKinematics();
  cout<<endl;
}//void Show()



void DebugOFF()
{
  Info("DebugOFF","");
  a->SetDebug(0);
  r->SetDebug(0);
  AliLoader::SetDebug(0);
}//void DebugOFF()

void DebugON()
{
  Info("DebugON","");  
  a->SetDebug(1);
  r->SetDebug(1);
  AliLoader::SetDebug(1);
}//void DebugON()


AliRun *a;
AliRunLoader *al;
AliLoader *rl,*tl,*il;

AliRICH *r;

Bool_t ReadAlice()
{
  Info("ReadAlice","Tring to read ALICE from SIMULATED FILE.");
  AliLoader::SetDebug(0);
  if(gAlice) delete gAlice;      
  if(!(al=AliRunLoader::Open("galice.root","AlicE","update"))){
    gSystem->Exec("rm -rf *.root *.dat");
    Error("ReadAlice","galice.root broken, removing all this garbage then init new one");
    new AliRun("gAlice","Alice experiment system");
    gAlice->SetDebug(-1);
    gAlice->Init("ConfigRich.C");
    r=(AliRICH*)gAlice->GetDetector("RICH");
    return kFALSE;
  }
  al->LoadgAlice();
  if(!gAlice) Fatal("ReadAlice","No gAlice in file");
  a=al->GetAliRun();
  a->SetDebug(0);    
//RICH      
  if(!(r=(AliRICH*)gAlice->GetDetector("RICH"))) Warning("RICH/menu.C::ReadAlice","No RICH in file");
  r->SetDebug(0);
  if(!(rl=al->GetLoader("RICHLoader")))          Warning("RICH/menu.C::ReadAlice","No RICH loader in file");        
        
  Info("ReadAlice","Run contains %i event(s)",gAlice->GetEventsPerRun());      
  return kTRUE;
}
//__________________________________________________________________________________________________
void GeoTest()
{

  TBRIK *pAliceBRIK=new TBRIK("aliceBRIK","ALICE mother volume","void",500,500,500);
  TBRIK *pArmBRIK=new TBRIK("armBRIK","RICH arm1","void",pRICH->GetSizeX(),pRICH->GetSizeY(),pRICH->GetSizeZ());
   
   TNode *pAliceNode=new TNode("aliceNode","Mother volume","aliceBRIK");
   pAliceNode->cd();

// ARM 1 LEFT      
   TRotation rot1;
   TVector3  v1(0,pRICH->GetOffset(),0);
   
         
   rot1.RotateX(pRICH->GetYZAngle()*kDegrad);        v1.RotateX(pRICH->GetYZAngle()*kDegrad);
   rot1.RotateZ(pRICH->GetXYAngle()*kDegrad);        v1.RotateZ(pRICH->GetXYAngle()*kDegrad);
   rot1.RotateZ(pRICH->GetRotAngle()*kDegrad);       v1.RotateZ(pRICH->GetRotAngle()*kDegrad);
         
   TRotMatrix *pArm1RotMatrix=new TRotMatrix("rotArm1","rotArm1",rot1.ThetaX()*kRaddeg, rot1.PhiX()*kRaddeg,
                                                                 rot1.ThetaY()*kRaddeg, rot1.PhiY()*kRaddeg,
                                                                 rot1.ThetaZ()*kRaddeg, rot1.PhiZ()*kRaddeg);

   TNode *pArm1Node=new TNode("arm1Node","Left arm","armBRIK",v1.X(),v1.Y(),v1.Z(),"rotArm1");
   arm1Node->SetLineColor(kRed);
   
// ARM 2 LEFT      
   TRotation rot2;
   TVector3  v2(0,pRICH->GetOffset(),0);
   
         
   rot2.RotateX( pRICH->YZAngle()*kDegrad);        v2.RotateX(pRICH->GetYZAngle()*kDegrad);
   rot2.RotateZ(-pRICH->XYAngle()*kDegrad);        v2.RotateZ(-pRICH->GetXYAngle()*kDegrad);
   rot2.RotateZ( pRICH->RotAngle()*kDegrad);        v2.RotateZ(pRICH->GetRotAngle()*kDegrad);
         
   TRotMatrix *pArm2RotMatrix=new TRotMatrix("rotArm2","rotArm2",rot2.ThetaX()*kRaddeg, rot2.PhiX()*kRaddeg,
                                                                 rot2.ThetaY()*kRaddeg, rot2.PhiY()*kRaddeg,
                                                                 rot2.ThetaZ()*kRaddeg, rot2.PhiZ()*kRaddeg);

   TNode *pArm2Node=new TNode("arm2Node","Left arm","armBRIK",v2.X(),v2.Y(),v2.Z(),"rotArm2");
   arm2Node->SetLineColor(kBlue);
   
   aliceNode->Draw();
}//void GeoTest()
//__________________________________________________________________________________________________
void PrintGeo(Float_t rotDeg=0)
{
  AliRICHParam *p=new AliRICHParam;  
  Double_t r=p->Offset();
  Double_t kP=p->AngleXY();
  Double_t kT=p->AngleYZ();
  Double_t kRot;
  
  if(rotDeg==0)
    kRot=p->AngleRot();
  else
    kRot=rotDeg*deg;
        
  cout<<endl;
  Double_t  phi=90*deg+kRot+kP,theta=90*deg+kT;
  Info("   menu for          1","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot+kP,theta=90*deg;
  Info("   menu for          2","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot,theta=90*deg-kT;    
  Info("   menu for          3","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  
  phi=90*deg+kRot,theta=90*deg;
  Info("   menu for          4","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));

  
  phi=90*deg+kRot,theta=90*deg+kT;
  Info("   menu for          5","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                                r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  
  phi=90*deg+kRot-kP,theta=90*deg;
  Info("   menu for          6","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot-kP,theta=90*deg+kT;
  Info("   menu for          7","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));

  delete p;
}//void PrintGeo()



//__________________________________________________________________________________________________
void TestGain()
{
  AliRICHParam *pParam=new AliRICHParam;
  AliRICHResponse *pRes=new AliRICHResponse;
  
  TLegend *pLegend=new TLegend(0.6,0.3,0.85,0.5);  
  TH1F *pH0=new TH1F("pH1","Gain",100,0,600); 
  TH1F *pH10=new TH1F("pH10","Gain",100,0,600);
  TH1F *pH20=new TH1F("pH20","Gain",100,0,600);
  TH1F *pH30=new TH1F("pH30","Gain",100,0,600);
  TH1F *pHold=new TH1F("pHold","Mip Charge",100,0,2000);
  for(int i=0;i<1000;i++){
    pH0 ->Fill(pParam->Gain(0));
    pH10->Fill(pParam->Gain(10));
    pH20->Fill(pParam->Gain(20));
    pH30->Fill(pParam->Gain(30));
    pHold->Fill(pRes->IntPH(30));
  }
  pH0->Draw();
  pH10->Draw("same");
  pH20->Draw("same");
  pH30->Draw("same");
  pHold->Draw("same");
  pLegend->AddEntry(pH0,"y=0");  
  pLegend->AddEntry(pH10,"y=10");pH10->SetLineColor(kRed);
  pLegend->AddEntry(pH20,"y=20");pH20->SetLineColor(kBlue);
  pLegend->AddEntry(pH30,"y=30");pH30->SetLineColor(kGreen);
  pLegend->AddEntry(pHold,"res");pHold->SetLineColor(kMagenta);
  pLegend->Draw();
}//void TestGain()
//__________________________________________________________________________________________________
void TestMipCharge()
{
  AliRICHParam *pParam=new AliRICHParam;
  AliRICHResponse *pRes=new AliRICHResponse;
  
  TLegend *pLegend=new TLegend(0.6,0.3,0.85,0.5);  
  TH1F *pH0= new TH1F("pH1", "Mip Charge",100,0,500); 
  TH1F *pH10=new TH1F("pH10","Mip Charge",100,0,500);
  TH1F *pH20=new TH1F("pH20","Mip Charge",100,0,500);
  TH1F *pH30=new TH1F("pH30","Mip Charge",100,0,500);
  TH1F *pHold=new TH1F("pHold","Mip Charge",100,0,500);
  for(int i=0;i<1000;i++){
    pH0 ->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,0));
    pH10->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,10));
    pH20->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,20));
    pH30->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,30));
    pHold->Fill(pRes->IntPH(0.5e-9,-30));
  }
  pH0->Draw();
  pH10->Draw("same");
  pH20->Draw("same");
  pH30->Draw("same");
  pHold->Draw("same");
  pLegend->AddEntry(pH0,"y=0");  
  pLegend->AddEntry(pH10,"y=10");pH10->SetLineColor(kRed);
  pLegend->AddEntry(pH20,"y=20");pH20->SetLineColor(kBlue);
  pLegend->AddEntry(pH30,"y=30");pH30->SetLineColor(kGreen);
  pLegend->AddEntry(pHold,"res");pHold->SetLineColor(kMagenta);
  pLegend->Draw();
}//void TestGain()
//__________________________________________________________________________________________________
void TestDigitsOLD()
{
  Info("TestDigitsOLD","Creating test digits.");
  rl->MakeTree("D");r->MakeBranch("D");


  Int_t t[10];
  Int_t c[10];
  Int_t d[5];
  t[0]=100;t[1]=200;t[2]=300;
  c[0]=10;c[1]=20;c[2]=30;


  d[0]=1;d[1]=1;d[2]=10;d[3]=3;d[4]=4;r->AddDigitOld(1,t,c,d);

  d[0]=2;d[1]=2;d[2]=10;d[3]=3;d[4]=4;r->AddDigitOld(2,t,c,d);
  d[0]=2;d[1]=3;d[2]=10;d[3]=3;d[4]=4;r->AddDigitOld(2,t,c,d);

  d[0]=2;d[1]=2;d[2]=100;d[3]=3;d[4]=4;r->AddDigitOld(3,t,c,d);
  d[0]=2;d[1]=3;d[2]= 50;d[3]=3;d[4]=4;r->AddDigitOld(3,t,c,d);
  d[0]=2;d[1]=4;d[2]=200;d[3]=3;d[4]=4;r->AddDigitOld(3,t,c,d);

  d[0]=2;d[1]=2;d[2]=100;d[3]=3;d[4]=4;r->AddDigitOld(6,t,c,d);
  d[0]=2;d[1]=3;d[2]= 50;d[3]=3;d[4]=4;r->AddDigitOld(6,t,c,d);
  d[0]=2;d[1]=4;d[2]=200;d[3]=3;d[4]=4;r->AddDigitOld(6,t,c,d);
  d[0]=2;d[1]=5;d[2]= 50;d[3]=3;d[4]=4;r->AddDigitOld(6,t,c,d);
  d[0]=2;d[1]=6;d[2]=300;d[3]=3;d[4]=4;r->AddDigitOld(6,t,c,d);
  d[0]=2;d[1]=7;d[2]= 50;d[3]=3;d[4]=4;r->AddDigitOld(6,t,c,d);

  rl->TreeD()->Fill();
  rl->WriteDigits("OVERWRITE");
  rl->UnloadDigits();
  r->ResetDigitsOld();
  Info("TestDigitsOLD","Stop.");
}//void TestDigits()
//__________________________________________________________________________________________________
void TestSD()
{
  Info("TestSD","Creating test sdigits.");
  rl->MakeTree("S");r->MakeBranch("S");
  
    for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
  
    if(!rl->TreeH()) rl->LoadHits();//from
    if(!rl->TreeS()) rl->MakeTree("S");    r->MakeBranch("S");//to
      
    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<3;iHitN++){//hits loop  ???
        AliRICHhit *pHit=r->Hits()->At(iHitN);        
        TVector3 globX3(pHit->X(),pHit->Y(),pHit->Z());        
        TVector3 locX3=r->C(pHit->C())->Glob2Loc(globX3);
        
        Int_t sector;
        Int_t iTotQdc=r->Param()->Loc2TotQdc(locX3,pHit->Eloss(),pHit->Pid(),sector);
        
        Int_t iPadXmin,iPadXmax,iPadYmin,iPadYmax;
        r->Param()->Loc2Area(locX3,iPadXmin,iPadYmin,iPadXmax,iPadYmax);
        cout<<"left-down=("<<iPadXmin<<","<<iPadYmin<<") right-up=("<<iPadXmax<<','<<iPadYmax<<')'<<endl;
        for(Int_t iPadY=iPadYmin;iPadY<=iPadYmax;iPadY++)
          for(Int_t iPadX=iPadXmin;iPadX<=iPadXmax;iPadX++){
            Float_t iPadQdc=iTotQdc*r->Param()->Loc2PadFrac(locX3,iPadX,iPadY);
            Int_t padx,pady; r->Param()->Loc2Pad(locX3.X(),locX3.Y(),padx,pady);
            cout<<"hit="<<iHitN<<" ("<<locX3.X()<<','<<locX3.Y()<<")("<<padx<<','<<pady<<") cur pad("<<iPadX<<","<<iPadY<<") qtot="<<iTotQdc<<" qfrac="<<r->Param()->Loc2PadFrac(locX3,iPadX,iPadY)<<endl;
          }
//            r->AddSdigit(pHit->C(),padx,pady,r->Param()->Local2PadQdc(localX3,padx,pady),pHit->GetTrack());
      }//hits loop
    }//prims loop
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop
  
  rl->UnloadHits();
  rl->UnloadSDigits();  

  rl->TreeS()->Fill();
  rl->WriteSDigits("OVERWRITE");
  rl->UnloadSDigits();
  cout<<endl;r->Sdigits()->Print();
  r->ResetSDigits();
  Info("TestSdigits","Stop.");
}//void TestSdigits()
//__________________________________________________________________________________________________
void TestClustersOLD()
{
  Info("TestClusters","Creating test clusters.");
  rl->MakeTree("R");r->MakeBranch("R");
  
  AliRICHRawCluster c;
  r->AddClusterOld(1,c);  
  rl->TreeR()->Fill();
  rl->WriteRecPoints("OVERWRITE");
  rl->UnloadRecPoints();
  r->ResetRawClusters();
  
  Info("TestClusters","Stop.");
}//void TestClustersOLD()
//__________________________________________________________________________________________________
void TestSeg()
{
  AliRICHParam *p=r->Param();
  Int_t padx,pady;
  Double_t x,y;
  Double_t dz=p->DeadZone();
  Double_t sx=p->SectorSizeX(); Double_t sy=p->SectorSizeY();  Double_t px=p->PcSizeX(); Double_t py=p->PcSizeY();
  cout<<endl;
  Info("  1-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-px/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 48-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2-dz , -py/2        ,padx,pady),padx,pady);
  Info(" 49-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 96-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 97-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2+dz , -py/2        ,padx,pady),padx,pady);
  Info("144-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad( px/2    , -py/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-px/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 48- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2-dz , -dz/2        ,padx,pady),padx,pady);
  Info(" 49- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 96- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 97- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2+dz , -dz/2        ,padx,pady),padx,pady);
  Info("144- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad( px/2    , -dz/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-px/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 48- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2-dz ,  dz/2        ,padx,pady),padx,pady);
  Info(" 49- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 96- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 97- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2+dz ,  dz/2        ,padx,pady),padx,pady);
  Info("144- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad( px/2    ,  dz/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-px/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 48-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2-dz ,  py/2        ,padx,pady),padx,pady);
  Info(" 49-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad(-sx/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 96-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 97-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad( sx/2+dz ,  py/2        ,padx,pady),padx,pady);
  Info("144-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad( px/2    ,  py/2        ,padx,pady),padx,pady);  
  cout<<endl;
  Info(" 73-160","sec=%i padx=%3i pady=%3i",p->Loc2Pad(    0    ,  py/2      ,padx,pady),padx,pady);    
  Info(" 73- 81","sec=%i padx=%3i pady=%3i",p->Loc2Pad(    0    ,  dz/2      ,padx,pady),padx,pady);    
  Info("0-0dead","sec=%i padx=%3i pady=%3i",p->Loc2Pad(    0    ,   0        ,padx,pady),padx,pady);    
  Info(" 73- 80","sec=%i padx=%3i pady=%3i",p->Loc2Pad(    0    , -dz/2      ,padx,pady),padx,pady);    
  Info(" 73-  1","sec=%i padx=%3i pady=%3i",p->Loc2Pad(    0    , -py/2      ,padx,pady),padx,pady);    
  cout<<endl;
  p->Pad2Loc(padx=  1,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 48,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 49,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 96,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 97,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx=144,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  cout<<endl;
  p->Pad2Loc(padx=  1,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 48,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 49,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 96,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 97,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx=144,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  cout<<endl;
  p->Pad2Loc(padx=  1,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 48,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 49,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 96,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 97,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx=144,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  cout<<endl;
  p->Pad2Loc(padx=  1,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 48,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 49,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 96,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx= 97,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Loc(padx=144,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
}//void TestSeg()
//__________________________________________________________________________________________________
void TestMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH test");
  pMenu->AddButton("Test segmentation",  "TestSeg()",         "Test AliRICHParam::L2P() method");
  pMenu->AddButton("Test transform",     "TestTransform()",   "Test ALiRICHChamber::L2G() and G2L methods");
  pMenu->AddButton("Test gain",          "TestGain()",        "Test AliRICHParam::Gain() method");
  pMenu->AddButton("Test MIP charge",    "TestMipCharge()",   "Test AliRICHParam::TotalCharge() method");
  pMenu->AddButton("Test sdigits",       "TestSD()",          "Create test set of sdigits");
  pMenu->AddButton("Test Digits OLD",    "TestDigitsOLD()",   "Create test set of OLD digits");
  pMenu->AddButton("Test Clusters OLD",  "TestClustersOLD()", "Create test set of OLD clusters");
  pMenu->Show();  
}//TestMenu()
//__________________________________________________________________________________________________
void GeoMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH draw");
  pMenu->AddButton("RICH Isometry", "gMC->Gdraw(\"ALIC\", 60,40,0, 10,10, 0.008,0.008)","Draws ALIC volume in isometry");
  pMenu->AddButton("RICH Front XY", "gMC->Gdraw(\"ALIC\", 0,0,0, 10,10, 0.01,0.01)","Draws ALIC volume in XY view");
  pMenu->AddButton("RICH Side YZ",  "gMC->Gdraw(\"ALIC\",90,180, 0, 10,10, 0.01,0.01)","Draws ALIC volume in YZ view");
  pMenu->AddButton("RICH Top XZ",   "gMC->Gdraw(\"ALIC\",90, 90, 0, 10,10, 0.01,0.01)","Draws ALIC volume in XZ view");
  pMenu->AddButton("Module Isometry","gMC->Gdraw(\"SRIC\", 30,60,0, 10,10, 0.1,0.1)","Draws SRIC volume in isometry");
  pMenu->AddButton("Module Front XY","gMC->Gdraw(\"SRIC\", 0,0,0, 10,10, 0.1,0.1)","Draws SRIC volume in XY view");
  pMenu->AddButton("Module Top XZ", "gMC->Gdraw(\"SRIC\",90, 90, 0, 10,10, 0.1,0.1)","Draws SRIC volume in XZ view");
  pMenu->AddButton("ALICE Tree", "((TGeant3*)gMC)->Gdtree(\"ALIC\")","Draws ALICE tree");      
  pMenu->AddButton("RICH Tree",  "((TGeant3*)gMC)->Gdtree(\"RICH\")","Draws RICH tree");      
  pMenu->AddButton("Geo test",  "GeoTest()",   "Draw RICH geo as a macro");
  pMenu->AddButton("Print ref", "PrintGeo()",  "Print RICH chambers default position");
  pMenu->AddButton("AliRICH::Print", "r->Print();", "Print RICH chambers default position");
  pMenu->AddButton("Test transform","TestTransform()","Test L2G and G2L methods");
  pMenu->Show();  
}//GeoMenu()
//__________________________________________________________________________________________________
void menu()
{ 
  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
  pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
  if(ReadAlice()){//it's from file, reconstruct
    pMenu->AddButton("hits->sdigits->digits","MainTranck()","Convert");
    
    pMenu->AddButton("hits->sdigits",    "H_SD()",       "Perform first phase converstion");
    pMenu->AddButton("sdigits->digits",  "SD_D()",       "Perform first phase converstion");
    pMenu->AddButton("digits->clusters", "D_C()",        "Perform first phase converstion");

    pMenu->AddButton("OLD Show","Show3()","Shows the structure of events in files");
    pMenu->AddButton("OLD specials->sdigits",          "OLD_S_SD()",       "Perform first phase converstion");
    pMenu->AddButton("OLD sdigits->digits",         "OLD_SD_D()","Perform second phase converstion");
    pMenu->AddButton("OLD digits->clusters",     "OLD_D_C()",  "Perform second phase converstion");
    
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
  }
  pMenu->AddButton("Geo submenu",     "GeoMenu()",            "Shows geomentry submenu");
  pMenu->AddButton("Geo GUI", "new G3GeometryGUI;","Create instance of G4GeometryGUI"); 
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",         "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Display Fast",    "DisplFast()",           "Display Fast");
  pMenu->AddButton("Quit",            ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//menu()
//__________________________________________________________________________________________________
