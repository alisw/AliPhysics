//globals for easy manual manipulations
AliRun *a;       
AliRunLoader *al;
AliLoader *rl;
AliRICH *r;
AliStack *s;

void ph(Int_t event=0)  {r->PrintHits(event);}    //utility print hits for 'event' event
void ps(Int_t event=0)  {r->PrintSDigits(event);} //utility print sdigits
void pd(Int_t event=0)  {r->PrintDigits(event);}  //utility print digits
void pc(Int_t event=0)  {r->PrintClusters(event);}//utility print clusters
void pt(Int_t event=0)  {r->PrintTracks(event);}  //utility print tracks
Int_t ne(Int_t event=0) {r->Nparticles(kElectron,event);}   //utility number of electrons
Int_t npi0(Int_t event=0) {r->Nparticles(kPi0,event);}   //utility number of electrons
Int_t npip(Int_t event=0) {r->Nparticles(kPiPlus,event);}   //utility number of electrons
Int_t npim(Int_t event=0) {r->Nparticles(kPiMinus,event);}   //utility number of electrons
Int_t nk0(Int_t event=0) {r->Nparticles(kK0,event);}   //utility number of electrons
Int_t nkp(Int_t event=0) {r->Nparticles(kKPlus,event);}   //utility number of electrons
Int_t nkm(Int_t event=0) {r->Nparticles(kKMinus,event);}   //utility number of electrons
//__________________________________________________________________________________________________
void pp(int tid)
{
  if(!al) return;
  al->LoadHeader();  al->LoadKinematics();
  
  if(tid<0||tid>=al->Stack()->GetNtrack())
    cout<<"Valid tid number is 0-"<<al->Stack()->GetNtrack()-1<<" for this event.\n";
  else
    PrintParticleInfo(tid);
  
  al->UnloadKinematics();  al->UnloadHeader();
}
//__________________________________________________________________________________________________
void PrintParticleInfo(int tid)
{
// Prints particle info for a given TID
  TParticle *p=al->Stack()->Particle(tid);
  cout<<p->GetName()<<"("<<tid<<")";
  if(p->GetMother(0)!=-1){cout<<" from "; PrintParticleInfo(p->GetMother(0));}
  else                   {cout<<endl;} 
}    
//__________________________________________________________________________________________________
Int_t prim(Int_t tid)
{
// Provides mother TID for the given TID
  al->LoadHeader();  al->LoadKinematics();
  
  if(tid<0||tid>=al->Stack()->GetNtrack())
    cout<<"Valid tid number is 0-"<<al->Stack()->GetNtrack()-1<<" for this event.\n";
  else
    while(1){
      TParticle *p=al->Stack()->Particle(tid);
      if(p->GetMother(0)==-1) break;
      tid=p->GetMother(0);
    }
  
  al->UnloadKinematics();  al->UnloadHeader();
  return tid;
}

//__________________________________________________________________________________________________
void Show()
{  
//  CreateHists();
  Info("","\n\n\n");
//load all trees  
  al->LoadHeader(); 
    al->LoadKinematics();  
      Bool_t isHits=!rl->LoadHits();  
        Bool_t isSdigits=!rl->LoadSDigits();  
          Bool_t isDigits=!rl->LoadDigits();//loaders
            Bool_t isClusters=!rl->LoadRecPoints();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    Int_t iNparticles=al->Stack()->GetNtrack();
    Int_t iNprims=al->Stack()->GetNprimary();
    
    Int_t iKPlusCounter=0,iKMinusCounter=0;    
    for(Int_t iParticleN=0;iParticleN<iNparticles;iParticleN++){//stack loop
      TParticle *pPart=al->Stack()->Particle(iParticleN);
      switch(pPart->GetPdgCode()){
        case kKPlus: iKPlusCounter++; break;
        case kKMinus:iKMinusCounter++; break;
      }
    }//stack loop
    
    Info("Show-STA","Evt %i->   %i particles %i primaries  %i K+ %i K-",
                     iEventN,   iNparticles,    iNprims,      iKPlusCounter,         iKMinusCounter);
    
    Int_t iHitsCounter=0;
    Int_t iNentries=rl->TreeH()->GetEntries();
    for(Int_t iEntryN=0;iEntryN<iNentries;iEntryN++){//TreeH loop
      rl->TreeH()->GetEntry(iEntryN);//get current entry (prim)          
      
      for(Int_t iHitN=0;iHitN<r->Hits()->GetEntries();iHitN++){//hits loop
        iHitsCounter++;
        AliRICHhit *pHit = (AliRICHhit*)r->Hits()->At(iHitN);//get current hit
        TParticle *pPart=al->Stack()->Particle(pHit->GetTrack());//get stack particle which produced the current hit
        FillContribs(pPart->GetPdgCode(),pHit->C(),kFALSE);
      }//hits loop
      
      if(iEntryN<7) Info("Show","Evt %i-> prim %4i has %4i hits from %s (,%7.2f,%7.2f)",
                  iEventN,iEntryN, r->Hits()->GetEntries(), pPart->GetName(), pPart->Theta()*TMath::RadToDeg(),pPart->Phi()*TMath::RadToDeg());
    }//TreeH loop
    Info("Show-HIT","Evt %i->   %i particles %i primaries  %i entries in TreeH %i hits",
                     iEventN,   iNparticles,    iNprims,      iNentries,         iHitsCounter);
    FillContribs(pPart->GetPdgCode(),pHit->C(),kTRUE);
    
    if(isSdigits){
      rl->TreeS()->GetEntry(0);
      Info("Show-SDI","Evt %i contains %5i sdigits",iEventN,r->SDigits()->GetEntries());
    }
    if(isDigits){
      rl->TreeD()->GetEntry(0);
      for(int i=1;i<=7;i++)
        Info("Show-DIG","Evt %i chamber %i contains %5i digits",
                                 iEventN,   i,           r->Digits(i)->GetEntries());
    }else
        Info("Show-DIG","There is no digits for this event");
    if(isClusters){
      rl->TreeR()->GetEntry(0);
      for(int i=1;i<=7;i++)
        Info("Show-CLU","Evt %i chamber %i contains %5i clusters",
                                 iEventN,   i,           r->Clusters(i)->GetEntries());
    }
    cout<<endl;
  }//events loop
//unload all trees    
  rl->UnloadHits();  
    if(isSdigits) rl->UnloadSDigits(); 
      if(isDigits) rl->UnloadDigits(); 
        if(isClusters) rl->UnloadRecPoints();
          al->UnloadHeader();
            al->UnloadKinematics();
  ShowHists();            
}//void Show()
//__________________________________________________________________________________________________



Bool_t ReadAlice()
{
  Info("ReadAlice","Tring to read ALICE from SIMULATED FILE.");
  if(gAlice) delete gAlice;      
  if(!(al=AliRunLoader::Open("galice.root","AlicE","update"))){
    gSystem->Exec("rm -rf *.root *.dat");
    Error("menu.C::ReadAlice","galice.root broken, removing all this garbage then init new one");
    new AliRun("gAlice","Alice experiment system");
    AliLog::SetModuleDebugLevel("RICH",1);
    gAlice->Init("Config.C");
    r=(AliRICH*)gAlice->GetDetector("RICH");
    a=gAlice; //for manual convinience
    return kFALSE;
  }
  al->LoadgAlice();//
  if(!gAlice) Fatal("menu.C::ReadAlice","No gAlice in file");
  a=al->GetAliRun();//provides pointer to AliRun object
  a->SetDebug(0);    
//RICH      
  if(!(r=(AliRICH*)gAlice->GetDetector("RICH"))) Warning("RICH/menu.C::ReadAlice","No RICH in file");
  if(!(rl=al->GetLoader("RICHLoader")))          Warning("RICH/menu.C::ReadAlice","No RICH loader in file");        
        
  Info("ReadAlice","Run contains %i event(s)",gAlice->GetEventsPerRun());      
  return kTRUE;
}
//__________________________________________________________________________________________________
void TestMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH test");
  pMenu->AddButton("Test segmentation"  ,"rp->TestSeg()"  ,"Test AliRICHParam segmentation methods");
  pMenu->AddButton("Test response"      ,"rp->TestResp()" ,"Test AliRICHParam response methods");
  pMenu->AddButton("Test transformation","rp->TestTrans()","Test AliRICHParam transformation methods");
  pMenu->AddButton("Test opticals"      ,".x Opticals.h"  ,"Test optical properties");
  pMenu->Show();  
}//TestMenu()
//__________________________________________________________________________________________________
void menu()
{ 
  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  if(ReadAlice()){//it's from file, show some info
    pMenu->AddButton("Show",            "Show()",             "Shows the structure of events in files");
    pMenu->AddButton("Display Fast",    "AliRICHDisplFast *d = new AliRICHDisplFast(); d->Exec();",        "Display Fast");
    pMenu->AddButton("Control Plots",   "r->ControlPlots()","Create some control histograms");
    
  }else{//it's aliroot, simulate
    pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
    pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
    pMenu->AddButton("Geo GUI",     "GeomGui()",       "Shows geometry"); 
    pMenu->AddButton("Read RAW",    "ReadRaw()",       "Read a list of digits from test beam file"); 
  }
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",         "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Quit",            ".q",                    "Close session");
  pMenu->Show();
}//menu()
//__________________________________________________________________________________________________
void DebugOFF(){  Info("DebugOFF","");  AliLog::SetGlobalDebugLevel(0);}
void DebugON() {  Info("DebugON","");   AliLog::SetGlobalDebugLevel(AliLog::kDebug);}
//__________________________________________________________________________________________________
TObjArray * CreateHists(Double_t pcut=0.9)
{
  TH2F *pPosH2    =new TH2F("pos"   ,"Pos mixture",5,0,5, 9,0,9); pPosH2->SetStats(0);
  TH2F *pNegH2    =new TH2F("neg"   ,"Neg mixture",5,0,5, 9,0,9); pNegH2->SetStats(0);
  TH2F *pPosCutH2 =new TH2F("poscut",Form("Pos mixture with P>%4.2f GeV",pcut),5,0,5, 9,0,9); pPosCutH2->SetStats(0);
  TH2F *pNegCutH2 =new TH2F("negcut",Form("Neg mixture with P>%4.2f GeV",pcut),5,0,5, 9,0,9); pNegCutH2->SetStats(0);
  pPosH2->GetXaxis()->SetBinLabel(1,"e^{+}");        pNegH2->GetXaxis()->SetBinLabel(1,"e^{-}");  
  pPosH2->GetXaxis()->SetBinLabel(2,"#mu^{+}");      pNegH2->GetXaxis()->SetBinLabel(2,"#mu^{-}");
  pPosH2->GetXaxis()->SetBinLabel(3,"#pi^{+}");      pNegH2->GetXaxis()->SetBinLabel(3,"#pi^{-}");
  pPosH2->GetXaxis()->SetBinLabel(4,"K^{+}");        pNegH2->GetXaxis()->SetBinLabel(4,"K^{-}");  
  pPosH2->GetXaxis()->SetBinLabel(5,"p^{+}");        pNegH2->GetXaxis()->SetBinLabel(5,"p^{-}");  
  
  pPosCutH2->GetXaxis()->SetBinLabel(1,"e^{+}");     pNegCutH2->GetXaxis()->SetBinLabel(1,"e^{-}");  
  pPosCutH2->GetXaxis()->SetBinLabel(2,"#mu^{+}");   pNegCutH2->GetXaxis()->SetBinLabel(2,"#mu^{-}");
  pPosCutH2->GetXaxis()->SetBinLabel(3,"#pi^{+}");   pNegCutH2->GetXaxis()->SetBinLabel(3,"#pi^{-}");
  pPosCutH2->GetXaxis()->SetBinLabel(4,"K^{+}");     pNegCutH2->GetXaxis()->SetBinLabel(4,"K^{-}");  
  pPosCutH2->GetXaxis()->SetBinLabel(5,"p^{+}");     pNegCutH2->GetXaxis()->SetBinLabel(5,"p^{-}");  
  
  pPosH2->GetYaxis()->SetBinLabel(1,"ch1");          pNegH2->GetYaxis()->SetBinLabel(1,"ch1");  
  pPosH2->GetYaxis()->SetBinLabel(2,"ch2");          pNegH2->GetYaxis()->SetBinLabel(2,"ch2");  
  pPosH2->GetYaxis()->SetBinLabel(3,"ch3");          pNegH2->GetYaxis()->SetBinLabel(3,"ch3");  
  pPosH2->GetYaxis()->SetBinLabel(4,"ch4");          pNegH2->GetYaxis()->SetBinLabel(4,"ch4");  
  pPosH2->GetYaxis()->SetBinLabel(5,"ch5");          pNegH2->GetYaxis()->SetBinLabel(5,"ch5");  
  pPosH2->GetYaxis()->SetBinLabel(6,"ch6");          pNegH2->GetYaxis()->SetBinLabel(6,"ch6");  
  pPosH2->GetYaxis()->SetBinLabel(7,"ch7");          pNegH2->GetYaxis()->SetBinLabel(7,"ch7");  
  pPosH2->GetYaxis()->SetBinLabel(8,"prim");         pNegH2->GetYaxis()->SetBinLabel(8,"prim");  
  pPosH2->GetYaxis()->SetBinLabel(9,"tot");          pNegH2->GetYaxis()->SetBinLabel(9,"tot");  

  pPosCutH2->GetYaxis()->SetBinLabel(1,"ch1");          pNegCutH2->GetYaxis()->SetBinLabel(1,"ch1");  
  pPosCutH2->GetYaxis()->SetBinLabel(2,"ch2");          pNegCutH2->GetYaxis()->SetBinLabel(2,"ch2");  
  pPosCutH2->GetYaxis()->SetBinLabel(3,"ch3");          pNegCutH2->GetYaxis()->SetBinLabel(3,"ch3");  
  pPosCutH2->GetYaxis()->SetBinLabel(4,"ch4");          pNegCutH2->GetYaxis()->SetBinLabel(4,"ch4");  
  pPosCutH2->GetYaxis()->SetBinLabel(5,"ch5");          pNegCutH2->GetYaxis()->SetBinLabel(5,"ch5");  
  pPosCutH2->GetYaxis()->SetBinLabel(6,"ch6");          pNegCutH2->GetYaxis()->SetBinLabel(6,"ch6");  
  pPosCutH2->GetYaxis()->SetBinLabel(7,"ch7");          pNegCutH2->GetYaxis()->SetBinLabel(7,"ch7");  
  pPosCutH2->GetYaxis()->SetBinLabel(8,"prim");         pNegCutH2->GetYaxis()->SetBinLabel(8,"prim");  
  pPosCutH2->GetYaxis()->SetBinLabel(9,"tot");          pNegCutH2->GetYaxis()->SetBinLabel(9,"tot");  
  TObjArray *pOA=new TObjArray;
  pOA->Add(pPosH2);pOA->Add(pNegH2);pOA->Add(pPosCutH2);pOA->Add(pNegCutH2);  
  return pOA;
}//ParticleContribs()
//__________________________________________________________________________________________________

void ShowHists()
{
  pC=new TCanvas("c1","Particle composition");
  pC->Divide(2,2);
  pC->cd(1);  pos->Draw("text"); gPad->SetGrid();
  pC->cd(2);  neg->Draw("text"); gPad->SetGrid();
  pC->cd(3);  poscut->Draw("text"); gPad->SetGrid();
  pC->cd(4);  negcut->Draw("text"); gPad->SetGrid();
}


//__________________________________________________________________________________________________
void CheckPR()
{
//Pattern recognition wirh Stack particles
  TFile *pFile = new TFile("$(HOME)/RPR.root","RECREATE","RICH Pattern Recognition");
  TNtupleD *hn = new TNtupleD("hn","ntuple","Pmod:Charge:TrackTheta:TrackPhi:TrackX:TrackY:MinX:MinY:ChargeMIP:ThetaCerenkov:NPhotons:MipIndex:Chamber:Particle");
  printf("\n\n");
  printf("Pattern Recognition done for event %5i",0);
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf);
  for(Int_t iEvtN=0;iEvtN<R()->GetLoader()->GetRunLoader()->GetNumberOfEvents();iEvtN++) {
    al->GetEvent(iEvtN);
    AliRICHTracker *tr = new AliRICHTracker();
    tr->RecWithStack(hn);
//    Info("CheckPR","Pattern Recognition done for event %i \b",iEvtN);
    printf("\b\b\b\b\b%5i",iEvtN+1);
  }
  printf("\n\n");
  pFile->Write();pFile->Close();
}

//__________________________________________________________________________________________________
void GeomGui()
{
  if(gGeoManager){ 
    gGeoManager->GetTopVolume()->Draw(); 
    AliRICHParam::ShowAxis();
  }else 
    new G3GeometryGUI;
}  

