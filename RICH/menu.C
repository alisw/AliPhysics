
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
  OLD_S_SD(); SD_D();   AliRICHClusterFinder *z=new AliRICHClusterFinder(r); z->Exec();//delete z;  
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
}//H_SD()
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
        if(id!=kBad&&r->Param()->IsOverTh(chamber,x,y,q))
           r->AddDigit(chamber,x,y,q,combiPid,tid);
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
void TestResponse()
{
  TCanvas *pC=new TCanvas("c","Amplification test",900,800);
  pC->Divide(1,2);
  pC->cd(1);
  TF1 *pF1=new TF1("f1","9e-6*pow(x,4)+2e-7*pow(x,3)-0.0316*pow(x,2)-3e-4*x+25.367",-70,70);
  pF1->Draw();
  
  pC->cd(2);
  
  const Int_t nPoints=8;
  THStack *pStack=new THStack("stack","photons");
  TLegend *pLeg=new TLegend(0.6,0.2,0.9,0.5,"legend");    
  TH1F *apH[nPoints];
  
  Double_t starty=AliRICHParam::DeadZone()/2;
  Double_t deltay=AliRICHParam::SectorSizeY()/nPoints;
  
  for(int i=0;i<nPoints;i++){
    apH[i]=new TH1F(Form("h%i",i),"Qdc for Photon;QDC;Counts",1000,0,1000); apH[i]->SetLineColor(i);
    pStack->Add(apH[i]);                 
    pLeg->AddEntry(apH[i],Form("@(0,%5.2f->%5.2f)",starty+i*deltay,starty+i*deltay-AliRICHParam::SectorSizeY()/2));
  }
        
  
  TVector3 x3(0,0,0);  
  Int_t sector=10;
//  AliRICHParam::ResetWireSag();
  for(Int_t i=0;i<10000;i++){//events loop
    for(int j=0;j<nPoints;j++){
      x3.SetY(starty-j*deltay);
      apH[j]->Fill(AliRICHParam::Loc2TotQdc(x3,400e-9,500000,sector));
    }
  }
  pStack->Draw("nostack");
  pLeg->Draw();
}//TestResponse()
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
//            r->AddSDigit(pHit->C(),padx,pady,r->Param()->Local2PadQdc(localX3,padx,pady),pHit->GetTrack());
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
  pMenu->AddButton("Test segmentation",  "TestSeg()",         "Test AliRICHParam segmentation methods");
  pMenu->AddButton("Test transform",     "TestTransform()",   "Test ALiRICHChamber methods");
  pMenu->AddButton("Test response",      "TestResponse()",    "Test AliRICHParam response methods");
  pMenu->AddButton("Test sdigits",       "TestSD()",          "Create test set of sdigits");
  pMenu->AddButton("Test digits OLD",    "TestDigitsOLD()",   "Create test set of OLD digits");
  pMenu->AddButton("Test clusters",      "TestC()",           "Create test set of clusters");
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
       
  if(ReadAlice()){//it's from file, reconstruct
    pMenu->AddButton("hits->sdigits->digits->clusters","MainTrank()","Convert");
    
    pMenu->AddButton("hits->sdigits",    "H_SD()",       "Perform first phase converstion");
    pMenu->AddButton("sdigits->digits",  "SD_D()",       "Perform first phase converstion");
    pMenu->AddButton("digits->clusters", "D_C()",        "Perform first phase converstion");

    pMenu->AddButton("Show","Show()","Shows the structure of events in files");
    pMenu->AddButton("OLD specials->sdigits",          "OLD_S_SD()",       "Perform first phase converstion");
    
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
  }
  pMenu->AddButton("Geo submenu",     "GeoMenu()",            "Shows geomentry submenu");
  pMenu->AddButton("Geo GUI", "new G3GeometryGUI;","Create instance of G4GeometryGUI"); 
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",         "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Display Fast",    "DisplFast()",           "Display Fast");
  pMenu->AddButton("Control Plots",   "ControlPlots()",        "Display some control histograms");
  pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
  pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
  pMenu->AddButton("Quit",            ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//menu()
//__________________________________________________________________________________________________
void DebugOFF(){  Info("DebugOFF","");  a->SetDebug(0);  r->SetDebug(0);  AliLoader::SetDebug(0);}
void DebugON() {  Info("DebugON","");   a->SetDebug(1);  r->SetDebug(1);  AliLoader::SetDebug(1);}
