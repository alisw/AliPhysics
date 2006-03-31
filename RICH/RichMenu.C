AliRun *a;    AliStack *s;  AliRunLoader *al; //globals for easy manual manipulations
AliRICH   *r; AliLoader    *rl;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  if(!(r=(AliRICH*)al->GetAliRun()->GetDetector("RICH")))     Warning("RICH/menu.C::ReadAlice","No RICH in file");
  if(!(rl=al->GetDetectorLoader("RICH")))           Warning("RICH/menu.C::ReadAlice","No RICH loader in file");        
}

//__________________________________________________________________________________________________
void RichMenu()
{   
  TControlBar *pMenu = new TControlBar("vertical","MAIN");
       
  if(AliceRead()){//it's from file, show some info
    pMenu->AddButton("Display single chambers"         ,"r->Display();"  , "Display Fast");
    pMenu->AddButton("Display ALL chambers"            ,"r->DisplayEvent(0,0);"  , "Display Fast");
    pMenu->AddButton("Recon with stack"                ,"AliRICHReconstructor::CheckPR(        )","Create RSR.root with ntuple hn");    
    pMenu->AddButton("RichAna no Recon"                ,"AliRICHReconstructor::RichAna(0,0,kFALSE)","Create RichAna.root with ntuple hn without PatRec");    
    pMenu->AddButton("RichAna with Recon"              ,"AliRICHReconstructor::RichAna(0,0,kTRUE )","Create RichAna.root with ntuple hn with PatRec");    
    pMenu->AddButton("HITS Print"                      ,"h();"      ,"To print hits: h()");
    pMenu->AddButton("HITS QA"                         ,"hqa()"     ,"QA plots for hits: hqa()");
    pMenu->AddButton("Print sdigits"                   ,"s();"      ,"To print sdigits: s()");
    pMenu->AddButton("Print digits"                    ,"d();"      ,"To print digits: d()");
    pMenu->AddButton("Clusters print"                  ,"c();"      ,"To print clusters: c()");  
    pMenu->AddButton("Clusters QA"                     ,"cqa();"    ,"Clusters QA: cqa() or AliRICHReconstructor::CluQA(al)");  
    pMenu->AddButton("Print ESD"                       ,"e();"      ,"To print ESD status");  
    pMenu->AddButton("Print Matrix"                    ,"m();"      ,"To print probability matrix");  
    pMenu->AddButton("Print occupancy"                 ,"r->OccupancyPrint(-1);" ,"To print occupancy");  
    pMenu->AddButton("Print event summary  "           ,"r->SummaryOfEvent();"   ,"To print a summary of the event");  
  }else{//it's aliroot, simulate
    pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
    pMenu->AddButton("Debug OFF",    "DebugOFF();",  "Switch debug on-off");   
    pMenu->AddButton("Run",          "a->Run(1)",  "Process!");
  }
  pMenu->AddButton("Test segmentation"  ,"rp->TestSeg()"  ,"Test AliRICHParam segmentation methods"     );
  pMenu->AddButton("Test response"      ,"rp->TestResp()" ,"Test AliRICHParam response methods"         );
  pMenu->AddButton("Test transformation","rp->TestTrans()","Test AliRICHParam transformation methods"   );
  pMenu->AddButton("Geo GUI"            ,"geo();"                                                          ,"Shows geometry"                             ); 
  pMenu->AddButton("Debug ON"           ,"AliLog::SetGlobalDebugLevel(AliLog::kDebug);"                    ,"Switch debug on"                            );   
  pMenu->AddButton("Debug OFF"          ,"AliLog::SetGlobalDebugLevel(0);"                                 ,"Switch debug off"                           );   
  pMenu->AddButton("Browser"            ,"new TBrowser;"                                                   ,"Start ROOT TBrowser"                        );
  pMenu->AddButton("Quit"               ,".q"                                                              ,"Close session"                              );
  pMenu->Show();
}//menu()
//__________________________________________________________________________________________________
void DebugOFF(){  Info("DebugOFF","");  AliLog::SetGlobalDebugLevel(0);}
void DebugON() {  Info("DebugON","");   AliLog::SetGlobalDebugLevel(AliLog::kDebug);}
//void geo (            )   { gGeoManager->SetVisOption(0);gGeoManager->GetTopVolume()->Draw(); AliRICHParam::DrawAxis();}
void geo()                  {  if(!gGeoManager) TGeoManager::Import("geometry.root");gGeoManager->GetTopVolume()->Draw();AliRICHParam::DrawAxis();}
  
void dis (Int_t evt=-1)   {r->Display      (evt);}                //utility display 
void dum (            )   {r->Dump         (   );}                //utility display 

void h   (Int_t evt=0 )   {r->HitPrint  (evt);}   //print hits for requested event
void hqa (            )   {r->HitQA     (   );}   //hits QA plots for all events 

void s   (Int_t evt=0 )   {r->SDigPrint (evt);}   //print sdigits for requested event

void d   (Int_t evt=0 )   {r->DigPrint  (evt);}   //print digits for requested event
void dqa (            )   {AliRICHReconstructor::DigQA     (al );}   //digits QA plots for all events


void c   (Int_t evt=0 )   {r->CluPrint  (evt);}   //print clusters for requested event
void cqa (            )   {AliRICHReconstructor::CluQA     (al );}   //clusters QA plots for all events
void ct  (            )   {AliRICHReconstructor::Test      (   );}   //test clusters by predifined list of digits

void t   (Int_t evt=0          )   {AliRICHParam::Stack(evt);}    
void tid (Int_t tid,Int_t evt=0)   {AliRICHParam::Stack(evt,tid);} 
void e   (                     )   {AliRICHTracker::EsdPrint();}                   
void m   (Double_t probCut=0.7 )   {AliRICHTracker::MatrixPrint(probCut);}                   

Int_t nem (Int_t evt=0) {AliRICHParam::StackCount(kElectron  ,evt);} //utility number of electrons
Int_t nep (Int_t evt=0) {AliRICHParam::StackCount(kPositron  ,evt);} //utility number of positrons
Int_t nmup(Int_t evt=0) {AliRICHParam::StackCount(kMuonPlus  ,evt);} //utility number of positive muons
Int_t nmum(Int_t evt=0) {AliRICHParam::StackCount(kMuonMinus ,evt);} //utility number of negative muons
Int_t npi0(Int_t evt=0) {AliRICHParam::StackCount(kPi0       ,evt);} //utility number of neutral pions 
Int_t npip(Int_t evt=0) {AliRICHParam::StackCount(kPiPlus    ,evt);} //utility number of positive pions
Int_t npim(Int_t evt=0) {AliRICHParam::StackCount(kPiMinus   ,evt);} //utility number of negative pions
Int_t nk0 (Int_t evt=0) {AliRICHParam::StackCount(kK0        ,evt);} //utility number of neutral kaons
Int_t nkp (Int_t evt=0) {AliRICHParam::StackCount(kKPlus     ,evt);} //utility number of positive kaons
Int_t nkm (Int_t evt=0) {AliRICHParam::StackCount(kKMinus    ,evt);} //utility number of negative kaons
Int_t npp (Int_t evt=0) {AliRICHParam::StackCount(kProton    ,evt);} //utility number of protons
Int_t npm (Int_t evt=0) {AliRICHParam::StackCount(kProtonBar ,evt);} //utility number of antiprotons
