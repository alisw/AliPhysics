
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



void Hits2Digits()
{
  TStopwatch sw;TDatime time;
  r->Hits2SDigits();r->SDigits2Digits();
  cout<<"\nInfo in <hits->digits>: Start time: ";time.Print();
  cout<<"Info in <hits->digits>: Stop  time: ";time.Set();  time.Print();
  cout<<"Info in <hits->digits>: Time  used: ";sw.Print();
}

void DigitsOLD2RawClustersOLD()
{
  AliRICHClusterFinder *z=new AliRICHClusterFinder(r);
  z->Exec();  
}

void Specials2DigitsOLD()
{
  Info("OLDspec2d","Start.");    
  
//  delete gAlice;
  
  AliRunDigitizer *pManager = new AliRunDigitizer(1,1);
  pManager->SetDebug(10);
  pManager->SetInputStream(0,"galice.root");
  new AliRICHDigitizer(pManager);
  pManager->Exec("deb");
  delete pManager;
  Info("OLDspec2d","Stop.");
}

void Specials2Sdigits()
{
  Info("Specials2Sdigits","Start.");
  
  rl->LoadHits(); 
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
    
    rl->MakeTree("S");  r->MakeBranch("S");
    r->ResetSdigits();  r->ResetSpecialsOld();

    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t i=0;i<r->Specials()->GetEntries();i++){//specials loop          
        Int_t padx= ((AliRICHSDigit*)r->Specials()->At(i))->PadX();
        Int_t pady= ((AliRICHSDigit*)r->Specials()->At(i))->PadY();
        Int_t qdc=  ((AliRICHSDigit*)r->Specials()->At(i))->QPad();
        Int_t hitN= ((AliRICHSDigit*)r->Specials()->At(i))->HitNumber()-1;//!!! important -1
        Int_t chamber=((AliRICHhit*)r->Hits()->At(hitN))->C();
        Int_t track=((AliRICHhit*)r->Hits()->At(hitN))->GetTrack();
        r->AddSdigit(chamber,padx+r->Param()->NpadsX()/2,pady+r->Param()->NpadsY()/2,qdc,track);
      }//specials loop
    }//prims loop
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop  
    rl->UnloadHits();     rl->UnloadSDigits();  
  Info("Specials2Sdigits","Stop.");    
}//Specials2Sdigits()
//__________________________________________________________________________________________________
void Hits2Sdigits()
{
  Info("Hits2Sdigits","Start.");
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
  
    if(!rl->TreeH()) rl->LoadHits();//from
    if(!rl->TreeS()) rl->MakeTree("S");    r->MakeBranch("S");//to
      
    NOT YET DONE!
    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<r->Hits()->GetEntries();iHitN++){//hits loop  
        AliRICHhit *pHit=r->Hits()->At(iHitN);
        r->Param()->G2Px
        r->AddSdigit(pHit->C(),padx+r->Param()->NpadsX()/2,pady+r->Param()->NpadsY()/2,qdc,track);
      }
    }
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop
  
  rl->UnloadHits();
  rl->UnloadSDigits();  
  Info("Hits2Sdigits","Stop.");  
}//Hits2Sdigits()
//__________________________________________________________________________________________________
void Sdigits2Digits()
{
  Info("Sdigits2Digits","Start.");  

  rl->LoadSDigits();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
    
    rl->MakeTree("D");r->MakeBranch("D"); //create TreeD with RICH branches 
    r->ResetSdigits();r->ResetDigits();//reset lists of sdigits and digits
    rl->TreeS()->GetEntry(0);  
    r->Sdigits()->Sort();
  
    Int_t kBad=-101;
    Int_t chamber,x,y,qdc,tr[3],id;
    chamber=x=y=qdc=tr[0]=tr[1]=tr[2]=id=kBad;
    Int_t iNdigitsPerPad=kBad;//how many sdigits for a given pad
        
    for(Int_t i=0;i<r->Sdigits()->GetEntries();i++){//sdigits loop (sorted)
      AliRICHdigit *pSdig=(AliRICHdigit*)r->Sdigits()->At(i);
      if(pSdig->Id()==id){//still the same pad
        iNdigitsPerPad++;
        qdc+=pSdig->Qdc();
        if(iNdigitsPerPad<=3)
          tr[iNdigitsPerPad-1]=pSdig->T(0);
        else
          Info("","More then 3 sdigits for the given pad");
      }else{//new pad, add the pevious one
        if(id!=kBad) r->AddDigit(chamber,x,y,qdc,tr[0],tr[1],tr[2]);//ch-xpad-ypad-qdc-tr1-2-3
        chamber=pSdig->C();x=pSdig->X();y=pSdig->Y();qdc=pSdig->Qdc();tr[0]=pSdig->T(0);id=pSdig->Id();
        iNdigitsPerPad=1;tr[1]=tr[2]=kBad;
      }
    }//sdigits loop (sorted)
  
    r->AddDigit(chamber,x,y,qdc,tr[0],tr[1],tr[2]);//add the last digit
        
    rl->TreeD()->Fill();  
    rl->WriteDigits("OVERWRITE");
  }//events loop
  rl->UnloadSDigits();     rl->UnloadDigits();  
  r->ResetSdigits();r->ResetDigits();//reset lists of sdigits and digits
  Info("Sdigits2Digits","Stop.");  
}



void Sdigits2DigitsOLD()
{
  Info("Sdigits2DigitsOLD","Start.");  

  rl->LoadSDigits();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
    
    rl->MakeTree("D");r->MakeBranch("D"); //create TreeD with RICH branches 
    r->ResetSdigits();r->ResetDigitsOld();//reset lists of sdigits and digits
    rl->TreeS()->GetEntry(0);  
    r->Sdigits()->Sort();
  
    Int_t kBad=-101;
    
    Int_t tr[3],q[3],dig[5]; for(Int_t i=0;i<3;i++) tr[i]=q[i]=kBad;    for(Int_t i=0;i<5;i++) dig[i]=kBad;        
    Int_t chamber=kBad,id=kBad,iNdigitsPerPad=kBad;//how many sdigits for a given pad
        
    for(Int_t i=0;i<r->Sdigits()->GetEntries();i++){//sdigits loop (sorted)
      AliRICHdigit *pSdig=(AliRICHdigit*)r->Sdigits()->At(i);
      if(pSdig->Id()==id){//still the same pad
        iNdigitsPerPad++;
        dig[2]+=pSdig->Qdc();//sum up qdc
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
  r->ResetSdigits();r->ResetDigitsOld();//reset lists of sdigits and digits
  Info("Sdigits2DigitsOLD","Stop.");  
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
      Info("Show-SDIGITS","Evt %i contains %5i sdigits",iEventN,r->Sdigits()->GetEntries());
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

Bool_t CheckAlice()
{
  if(gAlice){//it's aliroot
    if(gSystem->Exec("ls galice.root")){
      Info("CheckAlice","It's AliRoot, and no galice.root: SIMULATION");
      gAlice->SetDebug(-1);
      gAlice->Init("ConfigRich.C");
      r=(AliRICH*)gAlice->GetDetector("RICH");
      return kFALSE;
    }else{//galice.root is present we want to read alice from file
      ReadAlice();
      return kTRUE;
    }
  }else{//it's root with ALICE libs loaded
    ReadAlice();
    return kTRUE;
  }       
}//void CheckAlice()         

void ReadAlice()
{
  Info("ReadAlice","Reading ALICE from SIMULATED FILE.");
  AliLoader::SetDebug(0);
  if(gAlice) delete gAlice;      
  if(!(al=AliRunLoader::Open("galice.root","AlicE","update"))){
    gSystem->Exec("rm -rf *.root *.dat");
    Fatal("ReadAlice","galice.root broken, removing all this garbage");
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
}
//__________________________________________________________________________________________________
void RingViewer()
{
  gStyle->SetPalette(1);
  TCanvas *view=new TCanvas("Display","ALICE RICH Display",0,0,600,600);
  
  TH2F *pH2=new TH2F("pH2F","RICH DISPLAY",r->Param()->Nx(),0,r->Param()->Nx(),r->Param()->Ny(),0,r->Param()->Ny());
  pH2->SetStats(0);
  pH2->SetMaximum(100);

  Int_t Nevents = gAlice->GetEventsPerRun();
}
//______________________________________________________________________________
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
  AliRICHResponseV0 *pRes=new AliRICHResponseV0;
  
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
  AliRICHResponseV0 *pRes=new AliRICHResponseV0;
  
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
void TestSdigits()
{
  Info("TestSdigits","Creating test sdigits.");
  rl->MakeTree("S");r->MakeBranch("S");
//totally 19 must be trasformd to 6 digits
  r->AddSdigit(1,40,40,10,40); r->AddSdigit(1,40,40,10,41); r->AddSdigit(1,40,40,10,42); r->AddSdigit(1,40,40,10,43);
  r->AddSdigit(1,45,45,10,45); r->AddSdigit(1,45,45,10,45); r->AddSdigit(1,45,45,10,45); r->AddSdigit(1,45,45,10,45);

  r->AddSdigit(1,20,20,10,20);  r->AddSdigit(1,20,20,10,21);
  r->AddSdigit(1,25,25,10,25);  r->AddSdigit(1,25,25,10,25);
  
  r->AddSdigit(1,30,30,10,30); r->AddSdigit(1,30,30,10,31); r->AddSdigit(1,30,30,10,32);
  r->AddSdigit(1,35,35,10,35); r->AddSdigit(1,35,35,10,35); r->AddSdigit(1,35,35,10,35);
  

  r->AddSdigit(1,10,10,10,10);
  
  
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
  Float_t x,y;
  Float_t dz=p->DeadZone();
  Float_t sx=p->SectorSizeX(); Float_t sy=p->SectorSizeY();  Float_t px=p->PcSizeX();    Float_t py=p->PcSizeY();
  cout<<endl;
  Info("  1-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad(-px/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 48-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2-dz , -py/2        ,padx,pady),padx,pady);
  Info(" 49-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 96-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 97-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2+dz , -py/2        ,padx,pady),padx,pady);
  Info("144-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad( px/2    , -py/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad(-px/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 48- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2-dz , -dz/2        ,padx,pady),padx,pady);
  Info(" 49- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 96- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 97- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2+dz , -dz/2        ,padx,pady),padx,pady);
  Info("144- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad( px/2    , -dz/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad(-px/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 48- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2-dz ,  dz/2        ,padx,pady),padx,pady);
  Info(" 49- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 96- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 97- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2+dz ,  dz/2        ,padx,pady),padx,pady);
  Info("144- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad( px/2    ,  dz/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1-160","sec=%i padx=%3i pady=%3i",p->Local2Pad(-px/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 48-160","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2-dz ,  py/2        ,padx,pady),padx,pady);
  Info(" 49-160","sec=%i padx=%3i pady=%3i",p->Local2Pad(-sx/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 96-160","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 97-160","sec=%i padx=%3i pady=%3i",p->Local2Pad( sx/2+dz ,  py/2        ,padx,pady),padx,pady);
  Info("144-160","sec=%i padx=%3i pady=%3i",p->Local2Pad( px/2    ,  py/2        ,padx,pady),padx,pady);  
  cout<<endl;
  Info(" 73-160","sec=%i padx=%3i pady=%3i",p->Local2Pad(    0    ,  py/2      ,padx,pady),padx,pady);    
  Info(" 73- 81","sec=%i padx=%3i pady=%3i",p->Local2Pad(    0    ,  dz/2      ,padx,pady),padx,pady);    
  Info("0-0dead","sec=%i padx=%3i pady=%3i",p->Local2Pad(    0    ,   0        ,padx,pady),padx,pady);    
  Info(" 73- 80","sec=%i padx=%3i pady=%3i",p->Local2Pad(    0    , -dz/2      ,padx,pady),padx,pady);    
  Info(" 73-  1","sec=%i padx=%3i pady=%3i",p->Local2Pad(    0    , -py/2      ,padx,pady),padx,pady);    
  cout<<endl;
  p->Pad2Local(padx=  1,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 48,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 49,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 96,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 97,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx=144,pady=1,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  cout<<endl;
  p->Pad2Local(padx=  1,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 48,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 49,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 96,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 97,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx=144,pady=80,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  cout<<endl;
  p->Pad2Local(padx=  1,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 48,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 49,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 96,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 97,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx=144,pady=81,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  cout<<endl;
  p->Pad2Local(padx=  1,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 48,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 49,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 96,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx= 97,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
  p->Pad2Local(padx=144,pady=160,x,y);  cout<<"padx="<<padx<<" pady="<<pady<<" x="<<x<<" y="<<y<<endl;
}//void TestSeg()
//__________________________________________________________________________________________________
void TestMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH test");
  pMenu->AddButton("Test segmentation",  "TestSeg()",         "Test AliRICHParam::L2P() method");
  pMenu->AddButton("Test transform",     "TestTransform()",   "Test ALiRICHChamber::L2G() and G2L methods");
  pMenu->AddButton("Test gain",          "TestGain()",        "Test AliRICHParam::Gain() method");
  pMenu->AddButton("Test MIP charge",    "TestMipCharge()",   "Test AliRICHParam::TotalCharge() method");
  pMenu->AddButton("Test Sdigits",       "TestSdigits()",     "Create test set of sdigits");
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
  pMenu->AddButton("Geo GUI", "new G3GeometryGUI;","Create instance of G4GeometryGUI"); 
  pMenu->Show();  
}//GeoMenu()
//__________________________________________________________________________________________________
void menu()
{ 
  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
  pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
  if(CheckAlice()){//it's from file, reconstruct
    pMenu->AddButton("Hits->Sdigits->Digits","Hits2Digits()","Convert");
    pMenu->AddButton("Digits->Recos",         "Digits2Recos()","Convert");
    pMenu->AddButton("Show","Show3()","Shows the structure of events in files");
    pMenu->AddButton("Hits->Sdigits",    "Hits2Sdigits()",       "Perform first phase converstion");
    pMenu->AddButton("Specials->Sdigits","Specials2Sdigits()",    "Perform first phase converstion");
    pMenu->AddButton("Sdigits->Digits",  "Sdigits2Digits()",       "Perform first phase converstion");
    pMenu->AddButton("Digits->Clusters", "Digits2Clusters()",        "Perform first phase converstion");

    pMenu->AddButton("Sdigits->DigitsOLD",        "Sdigits2DigitsOLD()","Perform second phase converstion");
    pMenu->AddButton("DigitsOLD->RawClustersOLD", "DigitsOLD2RawClustersOLD()",  "Perform second phase converstion");
    
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
  }
  pMenu->AddButton("Geo submenu",     "GeoMenu()",            "Shows geomentry submenu");
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",         "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Display Fast",    "DisplFast()",           "Display Fast");
  pMenu->AddButton("Quit",            ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//menu()
//__________________________________________________________________________________________________
