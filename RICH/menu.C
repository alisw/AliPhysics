//__________________________________________________________________________________________________
void H_SD()
{
  Info("H_SD","Start.");
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);
  
    if(!rl->TreeH()) rl->LoadHits();  al->LoadHeader(); al->LoadKinematics();//from
    if(!rl->TreeS()) rl->MakeTree("S");    r->MakeBranch("S");//to
          
    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t iHitN=0;iHitN<r->Hits()->GetEntries();iHitN++){//hits loop  ???
        AliRICHhit *pHit=r->Hits()->At(iHitN);        
        
        TVector2 x2 = r->Param()->ShiftToWirePos(r->C(pHit->C())->Glob2Loc(pHit->OutX3()));        
        
        Int_t iTotQdc=r->Param()->TotQdc(x2,pHit->Eloss());
        
        Int_t iPadXmin,iPadXmax,iPadYmin,iPadYmax;
        r->Param()->Loc2Area(x2,iPadXmin,iPadYmin,iPadXmax,iPadYmax);
        cout<<"left-down=("<<iPadXmin<<","<<iPadYmin<<") right-up=("<<iPadXmax<<','<<iPadYmax<<')'<<endl;
        for(Int_t iPadY=iPadYmin;iPadY<=iPadYmax;iPadY++)
          for(Int_t iPadX=iPadXmin;iPadX<=iPadXmax;iPadX++){
            Double_t padQdc=iTotQdc*r->Param()->FracQdc(x2,iPadX,iPadY);
            cout<<padQdc<<endl;
            if(padQdc>0.1) r->AddSDigit(pHit->C(),iPadX,iPadY,padQdc,al->Stack()->Particle(pHit->GetTrack())->GetPdgCode(),pHit->GetTrack());
          }            
      }//hits loop
    }//prims loop
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop
  
  rl->UnloadHits(); al->UnloadHeader(); al->UnloadKinematics();
  rl->UnloadSDigits();  
  Info("H_SD","Stop.");  
}//H_SD()
//__________________________________________________________________________________________________

Int_t countContrib[7][3];

void ControlPlots()
{  
  Int_t iChamber=1;
  
  TFile *pFile = new TFile("$(HOME)/plots.root","RECREATE");   
  TH1F *pCqH1=new TH1F("ClusQ",   "Cluster Charge all chambers;q [QDC]",r->P()->MaxQdc(),0,r->P()->MaxQdc());
  TH1F *pCsH1=new TH1F("ClusSize","Cluster size all chambers;size [number of pads in cluster]",100,0,100);
  TH2F *pCmH2=new TH2F("ClusMap", "Cluster map;x [cm];y [cm]",1000,-r->P()->PcSizeX()/2,r->P()->PcSizeX()/2,
                                                             1000,-r->P()->PcSizeY()/2,r->P()->PcSizeY()/2);
  
  TH1F *pCqMipH1=new TH1F("MipClusQ",   "MIP Cluster Charge all chambers;q [QDC]",r->P()->MaxQdc(),0,r->P()->MaxQdc());
  TH1F *pCsMipH1=new TH1F("MipClusSize","MIP Cluster size all chambers;size [number of pads in cluster]",100,0,100);
  TH2F *pCmMipH2=new TH2F("MipClusMap", "MIP Cluster map;x [cm];y [cm]",1000,-r->P()->PcSizeX()/2,r->P()->PcSizeX()/2,
                                                             1000,-r->P()->PcSizeY()/2,r->P()->PcSizeY()/2);
  
  TH1F *pCqCerH1=new TH1F("CerClusQ",   "Cerenkov Cluster Charge all chambers;q [QDC]",r->P()->MaxQdc(),0,r->P()->MaxQdc());
  TH1F *pCsCerH1=new TH1F("CerClusSize","Cernekov Cluster size all chambers;size [number of pads in cluster]",100,0,100);
  TH2F *pCmCerH2=new TH2F("CerClusMap", "Cerenkov Cluster map;x [cm];y [cm]",1000,-r->P()->PcSizeX()/2,r->P()->PcSizeX()/2,
                                                             1000,-r->P()->PcSizeY()/2,r->P()->PcSizeY()/2);
  TH1F *pCqFeeH1=new TH1F("FeeClusQ",   "Feedback Cluster Charge all chambers;q [QDC]",r->P()->MaxQdc(),0,r->P()->MaxQdc());
  TH1F *pCsFeeH1=new TH1F("FeeClusSize","Feedback Cluster size all chambers;size [number of pads in cluster]",100,0,100);
  TH2F *pCmFeeH2=new TH2F("FeeClusMap", "Feedback Cluster map;x [cm];y [cm]",1000,-r->P()->PcSizeX()/2,r->P()->PcSizeX()/2,
                                                             1000,-r->P()->PcSizeY()/2,r->P()->PcSizeY()/2);
  Bool_t isClusters=!rl->LoadRecPoints();
  r->SetTreeAddress();  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);    
    if(isClusters){
      rl->TreeR()->GetEntry(0);
      Int_t iTotalClusters=0;
      for(int i=1;i<=7;i++){//chambers loop
        iTotalClusters+=r->Clusters(i)->GetEntries();    
        for(Int_t iClusterN=0;iClusterN<r->Clusters(i)->GetEntries();iClusterN++){//clusters loop
          AliRICHcluster *pClus=(AliRICHcluster*)r->Clusters(i)->At(iClusterN);
          
          countContrib[i-1][0] += pClus->Nmips();
          countContrib[i-1][1] += pClus->Ncerenkovs();
          countContrib[i-1][2] += pClus->Nfeedbacks();
              
          pCqH1->Fill(pClus->Q());             
          pCsH1->Fill(pClus->Size());           
          pCmH2->Fill(pClus->X(),pClus->Y());  
          
          if(pClus->IsPureMip()){ //Pure Mips
            pCqMipH1->Fill(pClus->Q());
            pCsMipH1->Fill(pClus->Size()); 
            pCmMipH2->Fill(pClus->X(),pClus->Y());
          }
          if(pClus->IsPureCerenkov()){ //Pure Photons
            pCqCerH1->Fill(pClus->Q());
            pCsCerH1->Fill(pClus->Size()); 
            pCmCerH2->Fill(pClus->X(),pClus->Y());
          }
          if(pClus->IsPureFeedback()){ //Pure Feedbacks
            pCqFeeH1->Fill(pClus->Q());
            pCsFeeH1->Fill(pClus->Size()); 
            pCmFeeH2->Fill(pClus->X(),pClus->Y());
          }
        }//clusters loop
      }//chambers loop
    }//isClusters
    Info("ControlPlots","Event %i processed.",iEventN);
  }//events loop 
  if(isClusters) rl->UnloadRecPoints();
  
  pFile->Write();
  delete pFile;
  for(Int_t i=0;i<7;i++)
    cout <<" chamber " << i+1 << " n. mips " << countContrib[i][0] << " n. ckovs " << countContrib[i][1] << " n. fdbks " << countContrib[i][2] << endl;
}//void ControlPlots()
//__________________________________________________________________________________________________
void MainTrank()
{
  TStopwatch sw;TDatime time;
  H_SD(); SD_D();   AliRICHClusterFinder *z=new AliRICHClusterFinder(r); z->Exec();//delete z;  
  cout<<"\nInfo in <MainTrank>: Start time: ";time.Print();
    cout<<"Info in <MainTrank>: Stop  time: ";time.Set();  time.Print();
    cout<<"Info in <MainTrank>: Time  used: ";sw.Print();
}
//__________________________________________________________________________________________________
void sh()
{
  if(rl->LoadHits()) return;
  
  Int_t iTotalHits=0;
  for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
    rl->TreeH()->GetEntry(iPrimN);      
    r->Hits()->Print();
    iTotalHits+=r->Hits()->GetEntries();
  }
  Info("sh","totally %i hits",iTotalHits);
  rl->UnloadHits();
}
//__________________________________________________________________________________________________
void ssp()
{
  if(rl->LoadHits()) return;
  
  Int_t iTotalSpecials=0;
  for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
    rl->TreeH()->GetEntry(iPrimN);      
    r->Specials()->Print();
    iTotalSpecials+=r->Specials()->GetEntries();
  }
  Info("ssp","totally %i specials",iTotalSpecials);
  rl->UnloadHits();
}
//__________________________________________________________________________________________________
void ss()
{
  if(rl->LoadSDigits()) return;
  rl->TreeS()->GetEntry(0);
  r->SDigits()->Print();
  Info("ss","totally %i sdigits",r->SDigits()->GetEntries());
  rl->UnloadSDigits();
}
//__________________________________________________________________________________________________
void sd()
{
  if(rl->LoadDigits()) return;
  rl->TreeD()->GetEntry(0);
  Int_t iTotalDigits=0;
  for(int i=1;i<=7;i++){
    r->Digits(i)->Print();
    iTotalDigits+=r->Digits(i)->GetEntries();
  }
  Info("sd","totally %i digits",iTotalDigits);
  rl->UnloadDigits();
}

void sc()
{
  if(rl->LoadRecPoints()) return;
  r->SetTreeAddress();
  rl->TreeR()->GetEntry(0);
  for(int i=1;i<=7;i++) r->Clusters(i)->Print();
  rl->UnloadRecPoints();
}

Double_t r2d = TMath::RadToDeg();
Double_t d2r = TMath::DegToRad();

void DisplFast(){ AliRICHDisplFast *d = new AliRICHDisplFast();  d->Exec();}  


void C_R()
{
  AliRICHRecon *detect = new AliRICHRecon("RICH patrec algorithm","Reconstruction");
    

  for (int nev=0; nev< a->GetEventsPerRun(); nev++) {    // Event Loop
    al->GetEvent(nev);
    cout <<endl<< "Processing event:" <<nev<<endl;
    detect->StartProcessEvent();
  } // event loop  
  delete detect;
}  
//__________________________________________________________________________________________________
void D_C()
{
  TStopwatch sw;TDatime time;

   AliRICHClusterFinder *z=new AliRICHClusterFinder(r); z->Exec();

   cout << endl;
   cout << "Info in Digits->Clusters: Start time: ";time.Print();
   cout << "Info in Digits->Clusters: Stop  time: ";time.Set();  time.Print();
   cout << "Info in Digits->Clusters: Time  used: ";sw.Print();
}
//__________________________________________________________________________________________________
void OLD_S_SD()
{
  Info("OLD_S_SD","Start.");  
  rl->LoadHits();   
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);    Info("OLD_S_SD","Processing event %i",iEventN);  
    
    rl->MakeTree("S");  r->MakeBranch("S");
    r->ResetSDigits();  r->ResetSpecialsOld();

    for(Int_t iPrimN=0;iPrimN<rl->TreeH()->GetEntries();iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      for(Int_t i=0;i<r->Specials()->GetEntries();i++){//specials loop
        Int_t padx=1+ ((AliRICHSDigit*)r->Specials()->At(i))->PadX()+r->Param()->NpadsX()/2;
        Int_t pady=1+ ((AliRICHSDigit*)r->Specials()->At(i))->PadY()+r->Param()->NpadsY()/2;
        Double_t q=  ((AliRICHSDigit*)r->Specials()->At(i))->QPad();
        
        Int_t hitN= ((AliRICHSDigit*)r->Specials()->At(i))->HitNumber()-1;//!!! important -1
        Int_t chamber=((AliRICHhit*)r->Hits()->At(hitN))->C();
        Int_t tid=((AliRICHhit*)r->Hits()->At(hitN))->GetTrack();
        Int_t pid=((AliRICHhit*)r->Hits()->At(hitN))->Pid();
        if(padx<1 || padx>r->Param()->NpadsX() ||pady<1 || pady>r->Param()->NpadsY())
          Warning("OLD_S_SD","pad is out of valid range padx= %i pady=%i event %i",padx,pady,iEventN);
        else
          r->AddSDigit(chamber,padx,pady,q,pid,tid);
      }//specials loop
    }//prims loop
    rl->TreeS()->Fill();
    rl->WriteSDigits("OVERWRITE");
  }//events loop  
    rl->UnloadHits();     rl->UnloadSDigits();  
  Info("OLD_S_SD","Stop.");    
}//OLD_S_SD()
//__________________________________________________________________________________________________
void SD_D()
{
  Info("SD_D","Start.");  
  extern Int_t kBad; 
  r->Param()->GenSigmaThMap();
  rl->LoadSDigits();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    al->GetEvent(iEventN);    cout<<"Event "<<iEventN<<endl;  
    rl->MakeTree("D");r->MakeBranch("D"); //create TreeD with RICH branches 
    r->ResetSDigits();r->ResetDigits();//reset lists of sdigits and digits
    rl->TreeS()->GetEntry(0);  
    r->SDigits()->Sort();
      
    Int_t combiPid,chamber,x,y,tid[3],id; Double_t q;
    Int_t iNdigitsPerPad;//how many sdigits for a given pad
    const int kBad=-101;//??? to be removed in code    
    for(Int_t i=0;i<r->SDigits()->GetEntries();i++){//sdigits loop (sorted)
      AliRICHdigit *pSdig=(AliRICHdigit*)r->SDigits()->At(i);
      if(pSdig->Id()==id){//still the same pad
        iNdigitsPerPad++;
        q+=pSdig->Q();
        combiPid+=pSdig->CombiPid();
        if(iNdigitsPerPad<=3)
          tid[iNdigitsPerPad-1]=pSdig->Tid(0);
//        else
//          Info("","More then 3 sdigits for the given pad");
      }else{//new pad, add the pevious one
        if(id!=kBad&&r->Param()->IsOverTh(chamber,x,y,q)) {
           r->AddDigit(chamber,x,y,q,combiPid,tid);
         }
        combiPid=pSdig->CombiPid();chamber=pSdig->C();id=pSdig->Id();
        x=pSdig->X();y=pSdig->Y();
        q=pSdig->Q();
        tid[0]=pSdig->Tid(0);
        iNdigitsPerPad=1;tid[1]=tid[2]=kBad;
      }
    }//sdigits loop (sorted)
  
    if(r->SDigits()->GetEntries()&&r->Param()->IsOverTh(chamber,x,y,q))
      r->AddDigit(chamber,x,y,q,combiPid,tid);//add the last digit
        
    rl->TreeD()->Fill();  
    rl->WriteDigits("OVERWRITE");
  }//events loop
  rl->UnloadSDigits();     rl->UnloadDigits();  
  r->ResetSDigits();r->ResetDigits();//reset lists of sdigits and digits
  Info("SD_D","Stop.");  
}//SD_D()
//__________________________________________________________________________________________________
void Show()
{  
  cout<<endl;
  al->LoadHeader();  al->LoadKinematics();
  
  rl->LoadHits();  
    Bool_t isSdigits=!rl->LoadSDigits();  
      Bool_t isClusters=!rl->LoadRecPoints();
        Bool_t isDigits=!rl->LoadDigits();//loaders
  
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
        Info("Show-DIGITS","Evt %i chamber %i contains %5i digits",
                                 iEventN,   i,           r->Digits(i)->GetEntries());
    }
    if(isClusters){
      rl->TreeR()->GetEntry(0);
      for(int i=1;i<=7;i++)
        Info("Show-CLUSTERS","Evt %i chamber %i contains %5i clusters",
                                 iEventN,   i,           r->Clusters(i)->GetEntries());
    }
    cout<<endl;
  }//events loop
  rl->UnloadHits();  
    if(isSdigits) rl->UnloadSDigits(); 
      if(isDigits) rl->UnloadDigits(); 
        if(isClusters) rl->UnloadRecPoints();
  al->UnloadHeader();
  al->UnloadKinematics();
  cout<<endl;
}//void Show()
//__________________________________________________________________________________________________


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
    gAlice->Init("Config.C");
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
}//PrintGeo()
//__________________________________________________________________________________________________

Double_t Gain(Double_t *x,Double_t *par)
{
  return AliRICHParam::GainSag(x[0],par[0]);
}

void TestResponse()
{
  TCanvas *pC=new TCanvas("c","Amplification test",900,800);
  pC->Divide(1,2);
  pC->cd(1);
  TF1 *pF1=new TF1("f1",Gain,-70,70,1);  pF1->SetParameters(1,1);pF1->SetParNames("Sector");
  TF1 *pF2=new TF1("f2",Gain,-70,70,1);  pF2->SetParameters(2,1);pF2->SetParNames("Sector");
  pF1->Draw();pF2->Draw("same");
  
  pC->cd(2);
  
  const Int_t nPoints=8;
  THStack *pStack=new THStack("stack","photons");
  TLegend *pLeg=new TLegend(0.6,0.2,0.9,0.5,"legend");    
  TH1F *apH[nPoints];
  
  Double_t starty=AliRICHParam::DeadZone()/2;
  Double_t deltay=AliRICHParam::SectorSizeY()/nPoints;
  
  for(int i=0;i<nPoints;i++){
    apH[i]=new TH1F(Form("h%i",i),"Qdc for Photon;QDC;Counts",500,0,500); apH[i]->SetLineColor(i);
    pStack->Add(apH[i]);                 
    pLeg->AddEntry(apH[i],Form("@(0,%5.2f->%5.2f)",starty+i*deltay,starty+i*deltay-AliRICHParam::SectorSizeY()/2));
  }
        
  
  TVector2 x2(0,0);  
//  AliRICHParam::ResetWireSag();
  for(Int_t i=0;i<10000;i++){//events loop
    for(int j=0;j<nPoints;j++){
      x2.Set(0,starty-j*deltay);
      apH[j]->Fill(AliRICHParam::TotQdc(x2,0));
    }
  }
  pStack->Draw("nostack");
  pLeg->Draw();
}//TestResponse()
//__________________________________________________________________________________________________
void TestSD()
{
  Info("TestSD","Creating test sdigits.");
  TVector3 hit(426.55,246.28,17.21);        
  TVector2 x2=r->C(4)->Glob2Loc(hit);        
  Int_t iTotQdc=r->Param()->TotQdc(x2,624e-9);        
  Int_t iPadXmin,iPadXmax,iPadYmin,iPadYmax;
  Int_t padx,pady;
  Int_t sec=r->Param()->Loc2Pad(x2,padx,pady);
  r->Param()->Loc2Area(x2,iPadXmin,iPadYmin,iPadXmax,iPadYmax);
  Info("TestSD","Initial hit (%7.2f,%7.2f,%7.2f)->(%7.2f,%7.2f)->(%4i,%4i,%4i) gives %i charge",
                          hit.X(),hit.Y(),hit.Z(),x2.X(),x2.Y(),sec,padx,pady,iTotQdc);
  
  cout<<"left-down=("<<iPadXmin<<","<<iPadYmin<<") right-up=("<<iPadXmax<<','<<iPadYmax<<')'<<endl;
  for(Int_t iPadY=iPadYmin;iPadY<=iPadYmax;iPadY++)
    for(Int_t iPadX=iPadXmin;iPadX<=iPadXmax;iPadX++)
       cout<<r->Param()->FracQdc(x2,iPadX,iPadY)<<endl;
  Info("TestSD","Stop.");
}//void TestSdigits()
//__________________________________________________________________________________________________
void TestC()
{
  Info("TestC","Creating test clusters.");
  rl->MakeTree("R");r->MakeBranch("R");
  
  AliRICHcluster c;
  c.AddDigit(new AliRICHdigit(1,20,21,200,1,2,3));
  c.AddDigit(new AliRICHdigit(1,22,21,250,1,2,3));
  c.CoG();
  
  r->AddCluster(c);  
  
  rl->TreeR()->Fill();
  rl->WriteRecPoints("OVERWRITE");
  rl->UnloadRecPoints();
  r->ResetClusters();
  
  Info("TestC","Stop.");
}//TestC()
//__________________________________________________________________________________________________
void TestSeg()
{
  AliRICHParam *p=r->Param();
  Int_t padx,pady,sec;
  Double_t x,y;
  Double_t eps=0.0000001;
  Double_t x1=-0.5*p->PcSizeX()+eps; Double_t x2=-0.5*p->SectorSizeX()-p->DeadZone()-eps; Double_t x3=-0.5*p->SectorSizeX()+eps;
  Double_t x6= 0.5*p->PcSizeX()-eps; Double_t x5= 0.5*p->SectorSizeX()+p->DeadZone()+eps; Double_t x4= 0.5*p->SectorSizeX()-eps;
  Double_t y1=-0.5*p->PcSizeY()+eps; Double_t y2=-0.5*p->DeadZone()-eps;
  Double_t y4= 0.5*p->PcSizeY()-eps; Double_t y3= 0.5*p->DeadZone()+eps;
  TVector2 v2;
  
  AliRICHParam::Print();
  
  sec=p->Loc2Pad(TVector2(x= 0,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 73-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x= 0,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 73- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x= 0,y= 0),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" dead  ","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x= 0,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 73- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x= 0,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 73-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)\n",x,y,sec,padx,pady,v2.X(),v2.Y());

  sec=p->Loc2Pad(TVector2(x=x1,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info("  1-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x2,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 48-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x3,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 49-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x4,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 96-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x5,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 97-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x6,y=y4),padx,pady); v2=p->Pad2Loc(padx,pady); Info("144-160","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)\n",x,y,sec,padx,pady,v2.X(),v2.Y());
  
  sec=p->Loc2Pad(TVector2(x=x1,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info("  1- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x2,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 48- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x3,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 49- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x4,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 96- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x5,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 97- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x6,y=y3),padx,pady); v2=p->Pad2Loc(padx,pady); Info("144- 81","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)\n",x,y,sec,padx,pady,v2.X(),v2.Y());
  
  sec=p->Loc2Pad(TVector2(x=x1,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info("  1- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x2,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 48- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x3,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 49- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x4,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 96- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x5,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 97- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x6,y=y2),padx,pady); v2=p->Pad2Loc(padx,pady); Info("144- 80","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)\n",x,y,sec,padx,pady,v2.X(),v2.Y());
  
  sec=p->Loc2Pad(TVector2(x=x1,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info("  1-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x2,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 48-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x3,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 49-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x4,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 96-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x5,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info(" 97-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)",x,y,sec,padx,pady,v2.X(),v2.Y());
  sec=p->Loc2Pad(TVector2(x=x6,y=y1),padx,pady); v2=p->Pad2Loc(padx,pady); Info("144-  1","(%7.2f,%7.2f)->(%5i,%5i,%5i) center=(%7.2f,%7.2f)\n",x,y,sec,padx,pady,v2.X(),v2.Y());
   
}//void TestSeg()
//__________________________________________________________________________________________________
void TestMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH test");
  pMenu->AddButton("Test segmentation",  "TestSeg()",         "Test AliRICHParam segmentation methods");
  pMenu->AddButton("Test response",      "TestResponse()",    "Test AliRICHParam response methods");
  pMenu->AddButton("Test sdigits",       "TestSD()",          "Create test set of sdigits");
  pMenu->AddButton("Test clusters",      "TestC()",           "Create test set of clusters");
  pMenu->AddButton("Test digits OLD",    "TestDigitsOLD()",   "Create test set of OLD digits");
  pMenu->Show();  
}//TestMenu()
//__________________________________________________________________________________________________
void menu()
{ 
  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  if(ReadAlice()){//it's from file, reconstruct
    pMenu->AddButton("hits->sdigits->digits->clusters","MainTrank()","Convert");
    
    pMenu->AddButton("hits->sdigits"    ,"H_SD()" ,"AliRICH::Hits2SDigits");
    pMenu->AddButton("sdigits->digits"  ,"SD_D()" ,"AliRICHDigitizer");
    pMenu->AddButton("digits->clusters" ,"D_C()"  ,"AliRICHClusterFinder");
    pMenu->AddButton("clusters->recos"  ,"C_R()"  ,"AliRICHRecon");

    pMenu->AddButton("Show",            "Show()","Shows the structure of events in files");
    pMenu->AddButton("Display Fast",    "DisplFast()",           "Display Fast");
    pMenu->AddButton("Control Plots",   "ControlPlots()",        "Display some control histograms");
    pMenu->AddButton("OLD specials->sdigits",          "OLD_S_SD()",       "Perform first phase converstion");
    
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
    pMenu->AddButton("Geo GUI", "new G3GeometryGUI;","Create instance of G4GeometryGUI"); 
  }
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",         "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
  pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
  pMenu->AddButton("Quit",            ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//menu()
//__________________________________________________________________________________________________
void DebugOFF(){  Info("DebugOFF","");  a->SetDebug(0);  r->SetDebug(0);  AliLoader::SetDebug(0);}
void DebugON() {  Info("DebugON","");   a->SetDebug(1);  r->SetDebug(1);  AliLoader::SetDebug(1);}
