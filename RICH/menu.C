//globals for easy manual manipulations
AliRun *a;       
AliRunLoader *al;
AliLoader *rl;
AliRICH *r;
AliStack *s;

AliRICH * R()    {return r;}
void ph(Int_t event=0) {R()->PrintHits(event);}    //utility print hits for 'event' event
void ps(Int_t event=0) {R()->PrintSDigits(event);} //utility print sdigits
void pd(Int_t event=0) {R()->PrintDigits(event);}  //utility print digits
void pc(Int_t event=0) {R()->PrintClusters(event);}//utility print clusters
void pt(Int_t event=0) {R()->PrintTracks(event);}  //utility print tracks
//void st(Int_t event=0) {R()->PrintStack(event);} //utiltity print stack
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
  R()->GetLoader()->GetRunLoader()->LoadHeader();  R()->GetLoader()->GetRunLoader()->LoadKinematics();
  
  if(tid<0||tid>=R()->GetLoader()->GetRunLoader()->Stack()->GetNtrack())
    cout<<"Valid tid number is 0-"<<R()->GetLoader()->GetRunLoader()->Stack()->GetNtrack()-1<<" for this event.\n";
  else
    while(1){
      TParticle *p=R()->GetLoader()->GetRunLoader()->GetAliRun()->Stack()->Particle(tid);
      if(p->GetMother(0)==-1) break;
      tid=p->GetMother(0);
    }
  
  R()->GetLoader()->GetRunLoader()->UnloadKinematics();  R()->GetLoader()->GetRunLoader()->UnloadHeader();
  return tid;
}

//__________________________________________________________________________________________________
void Show()
{  
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
                  iEventN,iEntryN, r->Hits()->GetEntries(), pPart->GetName(), pPart->Theta()*r2d,pPart->Phi()*r2d);
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
  
  TVector counters=r->Counters();
  
  counters(9)=counters(0); for(Int_t i=1;i<=8;i++) counters(9)-=counters(i);
  counters.Print();
  cout<<endl;
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
    pMenu->AddButton("Control Plots",   "R()->ControlPlots()","Create some control histograms");
    
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
void ReadRaw()
{
// First interation for TB raw data reader  
  char *sRawFileName="beam/run1822_4sigma.dat";
  FILE *pRawFile=fopen(sRawFileName,"r");
  if(!pRawFile){Error("ReadRaw","Something wrong with raw data file %s",sRawFileName); return;}
  
  Int_t iNpads=0,q=0,x=0,y=0;
  
  while(fscanf(pRawFile,"%i",&iNpads)){
    Info("ReadRaw","%i pads fired",iNpads);
    for(Int_t i=0;i<iNpads;i++)
      fscanf(pRawFile,"%i %i %i",&q,&x,&y);
  }
  fclose(pRawFile);
}
//__________________________________________________________________________________________________
void FillContribs(Int_t part,Int_t C, Bool_t print)
{
  static Int_t contrib[10][7][2];
  static Bool_t kEnter=kFALSE;

  C--; // chamber numbering from 1 to 7
  if(!kEnter) {
    for(Int_t i=0;i<10;i++) {
      for(Int_t j=0;j<7;j++) {
        for(Int_t k=0;k<2;k++) contrib[i][j][k]=0;
      }
    }
    kEnter=kTRUE;
  }

  if(print) {
    for(Int_t k=0;k<2;k++) {cout << "----------------Positives  to RICH ---------------" << endl;
                            cout << " chamber    1    2    3     4     5     6    7    " << endl;
                            cout << " -------------------------------------------------" << endl;
      for(Int_t i=0;i<4;i++) {
        if(i==0) cout << " e";
        if(i==1) cout << "pi";
        if(i==2) cout << " K";
        if(i==3) cout << " p";
        if(k==0) cout << "+         ";
        if(k==1) cout << "-         ";
        for(Int_t j=0;j<7;j++) {
          cout << contrib[i][j][k] << "    ";
        }
          cout << endl;
      }
    }
  }

  // +ves
  if(part==kPositron)    contrib[0][C][0]++;
  if(part==kPiPlus)      contrib[1][C][0]++;
  if(part==kKPlus)       contrib[2][C][0]++;
  if(part==kProton)      contrib[3][C][0]++;

  // -ves
  if(part==kElectron)    contrib[0][C][1]++;
  if(part==kPiMinus)     contrib[1][C][1]++;
  if(part==kKMinus)      contrib[2][C][1]++;
  if(part==kProtonBar)   contrib[3][C][1]++;
}
//__________________________________________________________________________________________________
void ParticleContribs()
{
  Bool_t isHits    =!rl->LoadHits();
  if(!isHits){Error("Exec","No hits. No contribs to calculate.");return;}
  
  TH1F *pDistH1 = new TH1F("distH1","length of charge track in amp gap;cm",100,0.,5.);

  al->LoadHeader();  al->LoadKinematics();



  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events Loop
    al->GetEvent(iEventN);
    
    AliStack *pStack=al->Stack();
    Info("","Event %i-> Totaly %i primaries %i",iEventN,pStack->GetNtrack(),pStack->GetNprimary());
    
    
    for(Int_t i=0;i<rl->TreeH()->GetEntries();i++){//prims loop
      rl->TreeH()->GetEntry(i);
      for(Int_t j=0;j<r->Hits()->GetEntries();j++){//hits loop for a given primary
        AliRICHhit *pHit = (AliRICHhit*)r->Hits()->At(j);
        TParticle *pParticle = al->Stack()->Particle(pHit->GetTrack());
//          if(pParticle->P()>1) FillContribs(pParticle->GetPdgCode(),pHit->C(),kFALSE);
        if(pParticle->GetPDG()->Charge()!=0) pDistH1->Fill(pHit->Length());
      }//hits loop
    }//prims loop
  }//events loop
  FillContribs(0,0,kTRUE);
  pDistH1->Draw();
}//ParticleContribs()
//__________________________________________________________________________________________________
void CheckPR()
{
//Pattern recognition wirh Stack particles
  TFile *pFile = new TFile("$(HOME)/RPR.root","RECREATE","RICH Pattern Recognition");
  TNtupleD *hn = new TNtupleD("hn","ntuple","Pmod:Charge:TrackTheta:TrackPhi:TrackX:TrackY:MinX:MinY:ChargeMIP:ThetaCerenkov:NPhotons:MipIndex");
  printf("\n\n");
    printf("Pattern Recognition done for event %5i",0);
  for(Int_t iEvtN=0;iEvtN<R()->GetLoader()->GetRunLoader()->GetNumberOfEvents();iEvtN++) {
    R()->GetLoader()->GetRunLoader()->GetEvent(iEvtN);
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

