#if !defined( __CINT__) || defined(__MAKECINT__)
#include <AliRun.h>
#include <AliStack.h>
#include <AliLoader.h>
#include <AliRunLoader.h>

#include "AliRICH.h"
#include "AliRICHDisplFast.h"
#endif

//globals for easy manual manipulations
AliRun *a;    AliStack *s;  AliRunLoader *al; 
AliRICH   *r; AliLoader    *rl,*vl;


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
  if(!p->IsPrimary()){cout<<" from "; PrintParticleInfo(p->GetFirstMother());}
  else                   {cout<<endl;} 
}    
//__________________________________________________________________________________________________
Int_t mother(Int_t tid)
{
// Who is the mother of given track TID?
  al->LoadHeader();  al->LoadKinematics();
  
  if(tid<0||tid>=al->Stack()->GetNtrack())
    cout<<"Valid tid number is 0-"<<al->Stack()->GetNtrack()-1<<" for this event.\n";
  else
    while(1){
      TParticle *p=al->Stack()->Particle(tid);
      if(p->IsPrimary()) break;
      tid=p->GetFirstMother();
    }
  
  al->UnloadKinematics();  al->UnloadHeader();
  return tid;
}
//__________________________________________________________________________________________________

Bool_t AliceRead()
{
  Info("ReadAlice","Tring to read ALICE from SIMULATED FILE...");
  if(gAlice){
    delete gAlice->GetRunLoader();
    delete gAlice;
  }      
  
  if(gSystem->Exec("ls galice.root>/dev/null")==256){//there is no galice.root in current directory
    AliceNew();
    RichGet();
    return kFALSE; //new session
  }else{
    if(!(al=AliRunLoader::Open())){//if not possible to read from galice.root, then remove grabage and reinvoke AliceRead()
      gSystem->Exec("rm -rf *.root *.dat");
      AliceRead();
    }
    al->LoadgAlice();//before this gAlice is 0;
    if(!gAlice) Fatal("menu.C::ReadAlice","No gAlice in file");
    a=al->GetAliRun();//provides pointer to AliRun object
    Info("AliceRead","Run contains %i event(s)",a->GetEventsPerRun());      
    RichGet();
    return kTRUE;   //old session opened from file
  }        
}//AliceRead()
//__________________________________________________________________________________________________
void AliceNew()
{
  Info("AliceNew","Init new session");
  new AliRun("gAlice","Alice experiment system");  gAlice->Init(); a=gAlice; al=gAlice->GetRunLoader();
}//AliceNew()    
//__________________________________________________________________________________________________
void RichGet()
{
  if(!(r=r()))     Warning("RICH/menu.C::ReadAlice","No RICH in file");
  if(!(rl=rl()))   Warning("RICH/menu.C::ReadAlice","No RICH loader in file");        
}

//__________________________________________________________________________________________________
void MenuRich()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH");
  pMenu->AddButton("Display single chambers"          ,"r->Display();"  , "Display Fast");
  pMenu->AddButton("Display ALL chambers"          ,"r->DisplayEvent(0,0);"  , "Display Fast");
  pMenu->AddButton("Print hits"       ,"r->HitsPrint();"      ,"????");
  pMenu->AddButton("Print sdigits"    ,"r->SDigitsPrint();"   ,"????");
  pMenu->AddButton("Print digits"     ,"r->DigitsPrint();"    ,"????");
  pMenu->AddButton("Print clusters"   ,"r->ClustersPrint();"  ,"????");  
  pMenu->AddButton("Print occupancy"  ,"r->OccupancyPrint();" ,"????");  
  pMenu->AddButton("Event Summary  "  ,"r->SummaryOfEvent();" ,"????");  
  pMenu->AddButton("Hits plots"       ,"r->ControlPlots()"    ,"????");
  pMenu->AddButton("Recon with stack" ,"r->CheckPR()"                                           , "Create RSR.root with ntuple hn");    
  pMenu->Show();  
}//TestMenu()
//__________________________________________________________________________________________________
void RichMenu()
{   
  TControlBar *pMenu = new TControlBar("vertical","MAIN");
       
  if(AliceRead()){//it's from file, show some info
    if(r) pMenu->AddButton("RICH submenu"    , "MenuRich()"               , "Show RICH submenu"       );
  }else{//it's aliroot, simulate
    pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
    pMenu->AddButton("Debug OFF",    "DebugOFF();",  "Switch debug on-off");   
    pMenu->AddButton("Run",          "a()->Run(1)",  "Process!");
  }
  pMenu->AddButton("Test segmentation"  ,"rp->TestSeg()"  ,"Test AliRICHParam segmentation methods"     );
  pMenu->AddButton("Test response"      ,"rp->TestResp()" ,"Test AliRICHParam response methods"         );
  pMenu->AddButton("Test transformation","rp->TestTrans()","Test AliRICHParam transformation methods"   );
  pMenu->AddButton("Test opticals"      ,".x Opticals.h"  ,"Test optical properties"                    );
  pMenu->AddButton("Geo GUI"            ,"GeomGui()"      ,"Shows geometry"                             ); 
  pMenu->AddButton("Browser"            ,"new TBrowser;"  ,"Start ROOT TBrowser"                        );
  pMenu->AddButton("Quit"               ,".q"             ,"Close session"                              );
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
    AliRICHParam::DrawAxis();
  }else 
    new G3GeometryGUI;
}  

AliRun    *a() {return al->GetAliRun();}                         //provides pointer to main AliRun object (aka gAlice)
AliRICH   *r() {return (AliRICH*)  a()->GetDetector("RICH");}    //provides pointer to RICH detector
AliLoader *rl(){return             al->GetLoader("RICHLoader");}

void rt(Int_t event=0)    {r->PrintTracks  (event);}                                                       //utility print tracks
Int_t nem(Int_t event=0)  {AliRICH::Nparticles(kElectron  ,event,al);} //utility number of electrons
Int_t nep(Int_t event=0)  {AliRICH::Nparticles(kPositron  ,event,al);} //utility number of positrons
Int_t nmup(Int_t event=0) {AliRICH::Nparticles(kMuonPlus  ,event,al);} //utility number of positive muons
Int_t nmum(Int_t event=0) {AliRICH::Nparticles(kMuonMinus ,event,al);} //utility number of negative muons
Int_t npi0(Int_t event=0) {AliRICH::Nparticles(kPi0       ,event,al);} //utility number of neutral pions 
Int_t npip(Int_t event=0) {AliRICH::Nparticles(kPiPlus    ,event,al);} //utility number of positive pions
Int_t npim(Int_t event=0) {AliRICH::Nparticles(kPiMinus   ,event,al);} //utility number of negative pions
Int_t nk0(Int_t event=0)  {AliRICH::Nparticles(kK0        ,event,al);} //utility number of neutral kaons
Int_t nkp(Int_t event=0)  {AliRICH::Nparticles(kKPlus     ,event,al);} //utility number of positive kaons
Int_t nkm(Int_t event=0)  {AliRICH::Nparticles(kKMinus    ,event,al);} //utility number of negative kaons
Int_t npp(Int_t event=0)  {AliRICH::Nparticles(kProton    ,event,al);} //utility number of protons
Int_t npm(Int_t event=0)  {AliRICH::Nparticles(kProtonBar ,event,al);} //utility number of antiprotons
