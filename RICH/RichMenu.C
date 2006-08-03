AliRun *a;    AliStack *s;  AliRunLoader *al; //globals for easy manual manipulations
AliRICH   *r; AliLoader    *rl;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GeomAlign(Bool_t isAlign)
{
  if(gGeoManager) delete gGeoManager;
  if(AliRICHParam::Instance()) delete AliRICHParam::Instance();
  if(isAlign) TGeoManager::Import("geometry.root");
  else        TGeoManager::Import("misaligned_geometry.root");
  AliRICHParam::Instance();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichMenu()
{   
  TString status="Status: ";
  if(gSystem->IsFileInIncludePath("galice.root")){
    status+="galice.root found";
    al=AliRunLoader::Open();                                                //try to open galice.root from current dir 
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    al->LoadgAlice(); a=al->GetAliRun();                                    //take new AliRun object from galice.root   
    rl=al->GetDetectorLoader("RICH");  r=(AliRICH*)a->GetDetector("RICH");  //get RICH object from galice.root
    
    status+=Form(" with %i event(s)",al->GetNumberOfEvents()); status+=(r)? " with RICH": " without RICH";
  }else  
    status+="No galice.root";
  
  TControlBar *pMenu = new TControlBar("horizontal",status.Data(),0,0);
    pMenu->AddButton("                     ","","");
    pMenu->AddButton(" General  "         ,"General()"  ,"general items which do not depend on any files");
    pMenu->AddButton("                     ","","");
    pMenu->AddButton(" Sim data "        ,"SimData()"  ,"items which expect to have simulated files"    );
    pMenu->AddButton("                     ","","");
    pMenu->AddButton(" Raw data "        ,"RawData()"  ,"items which expect to have raw files"          );
    pMenu->AddButton("                     ","","");
    pMenu->AddButton("Quit"            ,".q"         ,"close session"                                 );
  pMenu->Show();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void General()
{         
  TControlBar *pMenu = new TControlBar("vertical","Sim data",60,50);  
    pMenu->AddButton("Debug ON"           ,"don();"                    ,"Switch debug on-off"                        );   
    pMenu->AddButton("Debug OFF"          ,"doff();"                   ,"Switch debug on-off"                        );   
    pMenu->AddButton("Test segmentation"  ,"AliRICHParam::TestSeg();"  ,"Test AliRICHParam segmentation methods"     );
    pMenu->AddButton("Test response"      ,"AliRICHParam::TestResp();" ,"Test AliRICHParam response methods"         );
    pMenu->AddButton("Test transformation","AliRICHParam::TestTrans();","Test AliRICHParam transformation methods"   );
    pMenu->AddButton("Geo GUI"            ,"geo();"                    ,"Shows geometry"                             ); 
    pMenu->AddButton("Browser"            ,"new TBrowser;"             ,"Start ROOT TBrowser"                        );
  pMenu->Show();  
}//menu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimData()
{
  TControlBar *pMenu = new TControlBar("vertical","Sim",190,50);  
    pMenu->AddButton("Display single chambers"         ,"r->Display();"  , "Display Fast");
    pMenu->AddButton("Display ALL chambers"            ,"r->DisplayEvent(0,0);"  , "Display Fast");
    pMenu->AddButton("HITS Print"                      ,"h();"      ,"To print hits: h()");
    pMenu->AddButton("HITS QA"                         ,"hqa()"     ,"QA plots for hits: hqa()");
    pMenu->AddButton("Print sdigits"                   ,"s();"      ,"To print sdigits: s()");
    pMenu->AddButton("Print digits"                    ,"d();"      ,"To print digits: d()");
  pMenu->Show();         
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RawData()
{
  TControlBar *pMenu = new TControlBar("vertical","Raw",350,50);  
    pMenu->AddButton("ESD print"                       ,"ep();"                  ,"To print ESD info: ep()"         );  
    pMenu->AddButton("ESD QA"                          ,"eq();"                  ,"To draw ESD hists: eq()"         );  
    pMenu->AddButton("Clusters print"                  ,"cp();"                  ,"To print clusters: cp()"         );  
    pMenu->AddButton("Clusters QA"                     ,"cq();"                  ,"To draw clusters hists: cq()"    );  
    pMenu->AddButton("Print Matrix"                    ,"mp();"                  ,"To print prob matrix: mp()"      );  
    pMenu->AddButton("Print occupancy"                 ,"r->OccupancyPrint(-1);" ,"To print occupancy"              );  
    pMenu->AddButton("Print event summary  "           ,"r->SummaryOfEvent();"   ,"To print a summary of the event" );  
  pMenu->Show();         
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void doff(){  Printf("DebugOFF");  AliLog::SetGlobalDebugLevel(0);}
void don() {  Printf("DebugON");   AliLog::SetGlobalDebugLevel(AliLog::kDebug);}

void geo()                {  if(!gGeoManager) TGeoManager::Import("geometry.root");gGeoManager->GetTopVolume()->Draw();AliRICHParam::DrawAxis();}
  
void di  (Int_t evt=-1         )   {r->Display      (evt);}                //utility display 
void du  (                     )   {r->Dump         (   );}                //utility display 

void hp  (Int_t evt=0          )   {r->HitPrint  (evt);}   //print hits for requested event
void hq  (                     )   {r->HitQA     (   );}   //hits QA plots for all events 
void sp  (Int_t evt=0          )   {r->SDigPrint (evt);}   //print sdigits for requested event
void sq  (Int_t evt=0          )   {r->SDigPrint (evt);}   //print sdigits for requested event
void dp  (Int_t evt=0          )   {r->DigPrint  (evt);}   //print digits for requested event
void dq  (                     )   {AliRICHReconstructor::DigQA     (al );}   //digits QA plots for all events
void cp  (Int_t evt=0          )   {r->CluPrint  (evt);}                      //print clusters for requested event
void cq  (                     )   {AliRICHReconstructor::CluQA     (al );}   //clusters QA plots for all events
void ep  (                     )   {AliRICHTracker::EsdQA(1);} 
void eq  (                     )   {AliRICHTracker::EsdQA();}                   
void mp  (Double_t probCut=0.7 )   {AliRICHTracker::MatrixPrint(probCut);}                   


void t   (Int_t evt=0          )   {AliRICHParam::Stack(evt);}    
void tid (Int_t tid,Int_t evt=0)   {AliRICHParam::Stack(evt,tid);} 


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
