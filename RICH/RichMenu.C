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
Int_t nem(Int_t event=0)  {AliRICHDisplFast::Nparticles(kElectron  ,event,al);} //utility number of electrons
Int_t nep(Int_t event=0)  {AliRICHDisplFast::Nparticles(kPositron  ,event,al);} //utility number of positrons
Int_t nmup(Int_t event=0) {AliRICHDisplFast::Nparticles(kMuonPlus  ,event,al);} //utility number of positive muons
Int_t nmum(Int_t event=0) {AliRICHDisplFast::Nparticles(kMuonMinus ,event,al);} //utility number of negative muons
Int_t npi0(Int_t event=0) {AliRICHDisplFast::Nparticles(kPi0       ,event,al);} //utility number of electrons
Int_t npip(Int_t event=0) {AliRICHDisplFast::Nparticles(kPiPlus    ,event,al);} //utility number of electrons
Int_t npim(Int_t event=0) {AliRICHDisplFast::Nparticles(kPiMinus   ,event,al);} //utility number of electrons
Int_t nk0(Int_t event=0)  {AliRICHDisplFast::Nparticles(kK0        ,event,al);} //utility number of electrons
Int_t nkp(Int_t event=0)  {AliRICHDisplFast::Nparticles(kKPlus     ,event,al);} //utility number of electrons
Int_t nkm(Int_t event=0)  {AliRICHDisplFast::Nparticles(kKMinus    ,event,al);} //utility number of electrons
Int_t npp(Int_t event=0)  {AliRICHDisplFast::Nparticles(kProton    ,event,al);} //utility number of protons
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
  CreateHists();
  Info("","\n\n\n");
//load all trees  
  al->LoadHeader(); 
    al->LoadKinematics();  
      Bool_t isHits=!rl->LoadHits();  
        Bool_t isSdigits=!rl->LoadSDigits();  
          Bool_t isDigits=!rl->LoadDigits();//loaders
            Bool_t isClusters=!rl->LoadRecPoints();
  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    
    
    Info("Show-STA","Evt %i->   %i particles %i primaries  %i e- %i e+",
                     iEventN,   iNparticles,    iNprims,      iElectronCounter,      iPositronCounter);
    Info("Show-STA","Evt %i->   %i particles %i primaries  %i p+ %i p-",
                     iEventN,   iNparticles,    iNprims,      iPiPlusCounter,        iPiMinusCounter);
    Info("Show-STA","Evt %i->   %i particles %i primaries  %i K+ %i K-",
                     iEventN,   iNparticles,    iNprims,      iKPlusCounter,         iKMinusCounter);
    Info("Show-STA","Evt %i->   %i particles %i primaries  %i p  %i pbar",
                     iEventN,   iNparticles,    iNprims,      iProtonCounter,        iProtonBarCounter);
    
    Int_t iHitsCounter=0;
    Info("Show-HIT","Evt %i->   %i particles %i primaries  %i entries in TreeH %i hits",
                     iEventN,   iNparticles,    iNprims,      iNentries,         iHitsCounter);
    
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
  if(gAlice){
    delete gAlice->GetRunLoader();
    delete gAlice;
  }      
  if(!(al=AliRunLoader::Open())){//if not possible to read from galice.root, then create the new session
    gSystem->Exec("rm -rf *.root *.dat");
    Error("menu.C::ReadAlice","galice.root broken, removing all this garbage then init new one");
    new AliRun("gAlice","Alice experiment system");
    AliLog::SetModuleDebugLevel("RICH",1);
    gAlice->Init("Config.C");
    r=(AliRICH*)gAlice->GetDetector("RICH");
    a=gAlice; //for manual convinience
    return kFALSE;
  }
  al->LoadgAlice();//before this gAlice is 0;
  if(!gAlice) Fatal("menu.C::ReadAlice","No gAlice in file");
  a=al->GetAliRun();//provides pointer to AliRun object
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
void RichMenu()
{ 
  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  if(ReadAlice()){//it's from file, show some info
    pMenu->AddButton("Show",            "Show()",             "Shows the structure of events in files");
    pMenu->AddButton("Display Fast",    "AliRICHDisplFast *d = new AliRICHDisplFast(); d->Exec();",        "Display Fast");
    pMenu->AddButton("Control Plots",   "ControlPlots()",     "Create some control histograms");
    pMenu->AddButton("Recon with stack","r->CheckPR()",     "Create RSR.root with ntuple hn");    
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
void GeomGui()
{
  if(gGeoManager){ 
    gGeoManager->GetTopVolume()->Draw(); 
    AliRICHParam::ShowAxis();
  }else 
    new G3GeometryGUI;
}  
//__________________________________________________________________________________________________
void ControlPlots()
{
  r->ControlPlots();
  new TBrowser;
}
//__________________________________________________________________________________________________
void ParticleComposition()
{
  TH2F *pFluxH2    =new TH2F("flux","Charged flux for central Hijing event with Vertex<470 P>5MeV for electrons and positrons, P>1GeV for others",10,-5,5, 10,0,10); pFluxH2->SetStats(0);
  pFluxH2->GetXaxis()->SetBinLabel(1,Form("p^{-}>%dGeV"   ,cutPproton));        
  pFluxH2->GetXaxis()->SetBinLabel(2,Form("K^{-}>%dGeV"   ,cutPkaonminus));        
  pFluxH2->GetXaxis()->SetBinLabel(3,Form("#pi^{-}>%dGeV" , ));      
  pFluxH2->GetXaxis()->SetBinLabel(4,Form("#mu^{-}>%dGeV" , ));      
  pFluxH2->GetXaxis()->SetBinLabel(5,Form("e^{+}>%dGeV"   ,));        
  
  pFluxH2->GetXaxis()->SetBinLabel(6,Form("e^{-}>%dGeV"));        
  pFluxH2->GetXaxis()->SetBinLabel(7,Form("#mu^{+}>%dGeV"));      
  pFluxH2->GetXaxis()->SetBinLabel(8,Form("#pi^{+}>%dGeV"));      
  pFluxH2->GetXaxis()->SetBinLabel(9,Form("K^{+}>%dGeV"));        
  pFluxH2->GetXaxis()->SetBinLabel(10,Form("p^{+}>%dGeV"));        
  
  pFluxH2->GetYaxis()->SetBinLabel(1,"sum");  
  pFluxH2->GetYaxis()->SetBinLabel(2,"ch1");  
  pFluxH2->GetYaxis()->SetBinLabel(3,"ch2");  
  pFluxH2->GetYaxis()->SetBinLabel(4,"ch3");  
  pFluxH2->GetYaxis()->SetBinLabel(5,"ch4");  
  pFluxH2->GetYaxis()->SetBinLabel(6,"ch5");  
  pFluxH2->GetYaxis()->SetBinLabel(7,"ch6");  
  pFluxH2->GetYaxis()->SetBinLabel(8,"ch7");  
  pFluxH2->GetYaxis()->SetBinLabel(9,"prim"); 
  pFluxH2->GetYaxis()->SetBinLabel(10,"tot");  

  TH1F *pElecP=new TH1F("Pelec","Electrons made hit in RICH;p [GeV]",1000,-10,10); 
  TH1F *pMuonP=new TH1F("Pmuon","Muons made hit in RICH;p [GeV]"    ,1000,-10,10); 
  TH1F *pPionP=new TH1F("Ppion","Pions made hit in RICH;p [GeV]"    ,1000,-10,10); 
  TH1F *pKaonP=new TH1F("Pkaon","Kaon made hit in RICH;p [GeV]"     ,1000,-10,10); 
  TH1F *pProtP=new TH1F("Pprot","Protons made hit in RICH;p [GeV]"  ,1000,-10,10); 
  
  al->LoadHeader(); 
  al->LoadKinematics();  
  Int_t iNparticles=al->Stack()->GetNtrack();
  Int_t iNprims=al->Stack()->GetNprimary();
  
  for(Int_t iParticleN=0;iParticleN<iNparticles;iParticleN++){//stack loop
    TParticle *pPart=al->Stack()->Particle(iParticleN);

    if(iParticleN%10000==0) Info("Show-STA"," %i particles read",iParticleN);
    
    switch(pPart->GetPdgCode()){
      case kPositron:  pFluxH2->Fill(-1,9);  break;
      case kElectron:  pFluxH2->Fill( 0,9);  break;
      
      case kMuonMinus: pFluxH2->Fill(-2,9);  break;
      case kMuonPlus:  pFluxH2->Fill( 1,9);  break;
      
      case kPiMinus:   pFluxH2->Fill(-3,9);  break;
      case kPiPlus:    pFluxH2->Fill( 2,9);  break;
      
      case kKMinus:    pFluxH2->Fill(-4,9);  break;
      case kKPlus:     pFluxH2->Fill( 3,9);  break;
      
      case kProtonBar: pFluxH2->Fill(-5,9);  break;
      case kProton:    pFluxH2->Fill( 4,9);  break;            
    }//switch
  }//stack loop

  
  
    
  rl->LoadHits(); 
    
        
  for(Int_t iEntryN=0;iEntryN < rl->TreeH()->GetEntries();iEntryN++){//TreeH loop
    rl->TreeH()->GetEntry(iEntryN);//get current entry (prim)                
    for(Int_t iHitN=0;iHitN<r->Hits()->GetEntries();iHitN++){//hits loop
      AliRICHhit *pHit = (AliRICHhit*)r->Hits()->At(iHitN);//get current hit
      TParticle  *pPart=al->Stack()->Particle(pHit->GetTrack());//get stack particle which produced the current hit
      
      if(TMath::Sqrt(pPart->Vx()*pPart->Vx()+pPart->Vy()*pPart->Vy()+pPart->Vz()*pPart->Vz()) > 470) continue; //cut on vertex position
      if(pPart->GetPdgCode()==kElectron && pPart->P()<0.005)                                         continue; //cut on electron momentum 5MeV
      if(pPart->GetPdgCode()==kPositron && pPart->P()<0.005)                                         continue; //cut on electron momentum 5MeV
      if(pPart->GetPdgCode()!=kElectron && pPart->P()<1)                                             continue; //cut on others momentum 1GeV
      
      switch(pPart->GetPdgCode()){
        case kPositron : pElecP->Fill(-pPart->P()); pFluxH2->Fill(-1,pHit->C());break;
        case kElectron : pElecP->Fill( pPart->P()); pFluxH2->Fill( 0,pHit->C());break;
         
        case kMuonPlus : pMuonP->Fill( pPart->P()); pFluxH2->Fill(-2,pHit->C());break;
        case kMuonMinus: pMuonP->Fill(-pPart->P()); pFluxH2->Fill( 1,pHit->C());break;
                     
        case kPiMinus  : pPionP->Fill(-pPart->P()); pFluxH2->Fill(-3,pHit->C());break;
        case kPiPlus   : pPionP->Fill( pPart->P()); pFluxH2->Fill( 2,pHit->C());break;
           
        case kKPlus    : pKaonP->Fill( pPart->P()); pFluxH2->Fill(-4,pHit->C());break;
        case kKMinus   : pKaonP->Fill(-pPart->P()); pFluxH2->Fill( 3,pHit->C());break;
           
        case kProton   : pProtP->Fill( pPart->P()); pFluxH2->Fill(-5,pHit->C());break;
        case kProtonBar: pProtP->Fill(-pPart->P()); pFluxH2->Fill( 4,pHit->C());break;
      }
    }//hits loop      
  }//TreeH loop
                        
  rl->UnloadHits();  
  al->UnloadHeader(); 
  al->UnloadKinematics();  
  
  
  for(Int_t i=1;i<=pFluxH2->GetNbinsX();i++){
    Stat_t sum=0;
    for(Int_t j=2;j<=8;j++)    sum+=pFluxH2->GetBinContent(i,j);    
    pFluxH2->SetBinContent(i,1,sum);
  }
  
  TCanvas *pC1=new TCanvas("canvas1",Form("Event Nprims=%i",iNprims),1000,900);
  pFluxH2->Draw("text");  gPad->SetGrid();
  
  new TCanvas("celec","",200,100); pElecP->Draw();
  new TCanvas("cmuon","",200,100); pMuonP->Draw();
  new TCanvas("cpion","",200,100); pPionP->Draw();
  new TCanvas("ckaon","",200,100); pKaonP->Draw();
  new TCanvas("cprot","",200,100); pProtP->Draw();
  
}









