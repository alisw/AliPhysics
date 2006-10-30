AliRun *a;    AliStack *s;  AliRunLoader *al;   TGeoManager *g; //globals for easy manual manipulations
AliRICH   *r; AliLoader    *rl; AliRICHParam *rp;
Bool_t isGeomType=kFALSE;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GetParam()
{
  isGeomType=!isGeomType;
  if(g) delete g;  if(rp) delete rp; //delete current TGeoManager and AliRICHParam
  if(isGeomType) g=TGeoManager::Import("geometry.root");
  else           g=TGeoManager::Import("misaligned_geometry.root");
  rp=AliRICHParam::Instance();
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
  
  GetParam();
  
  TControlBar *pMenu = new TControlBar("horizontal",status.Data(),0,0);
    pMenu->AddButton("                     ","","");
    pMenu->AddButton(" General  "           ,"General()"  ,"general items which do not depend on any files");
    pMenu->AddButton("                     ",""           ,"");
    pMenu->AddButton(" Sim data "           ,"SimData()"  ,"items which expect to have simulated files"    );
    pMenu->AddButton("                     ",""           ,"");
    pMenu->AddButton(" Raw data "           ,"RawData()"  ,"items which expect to have raw files"          );
    pMenu->AddButton("                     ","print()"    ,"");
    pMenu->AddButton("Test"                 ,"Test()"     ,"all test utilities");
    pMenu->AddButton("                     ","GetParam()" ,"");
    pMenu->AddButton("Quit"                 ,".q"         ,"close session"                                 );
  pMenu->Show();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void General()
{         
  TControlBar *pMenu = new TControlBar("vertical","Sim data",60,50);  
    pMenu->AddButton("Debug ON"           ,"don();"                    ,"Switch debug on-off"                        );   
    pMenu->AddButton("Debug OFF"          ,"doff();"                   ,"Switch debug on-off"                        );   
    pMenu->AddButton("Geo GUI"            ,"geo();"                    ,"Shows geometry"                             ); 
    pMenu->AddButton("Browser"            ,"new TBrowser;"             ,"Start ROOT TBrowser"                        );
  pMenu->Show();  
}//menu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimData()
{
  TControlBar *pMenu = new TControlBar("vertical","Sim",190,50);  
    pMenu->AddButton("Display ALL chambers"            ,"ed();"     , "Display Fast");
    pMenu->AddButton("HITS QA"                         ,"hqa()"     ,"QA plots for hits: hqa()");
    
    pMenu->AddButton("Print hits"                      ,"hp();"      ,"To print hits:     hp(evt)");
    pMenu->AddButton("Print sdigits"                   ,"sp();"      ,"To print sdigits:  sp(evt)");
    pMenu->AddButton("Print digits"                    ,"dp();"      ,"To print digits:   dp(evt)");
    pMenu->AddButton("Print clusters"                  ,"cp();"      ,"To print clusters: cp(evt)");
    
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
}//RawData
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Test()
{         
  TControlBar *pMenu = new TControlBar("vertical","Test",400,50);  
    pMenu->AddButton("Hits->Digits"       ,"thd();"                    ,"test hits->sdigits->digits"                 );   
    pMenu->AddButton("Segmentation"       ,"ts()"                      ,"test segmentation methods"                  );
    pMenu->AddButton("Test response"      ,"AliRICHParam::TestResp();" ,"Test AliRICHParam response methods"         );
    pMenu->AddButton("Test transformation","AliRICHParam::TestTrans();","Test AliRICHParam transformation methods"   );
    pMenu->AddButton("Test Recon"         ,"rec();"                    ,"Test AliRICHRecon"                          );
  pMenu->Show();  
}//Test()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void doff(){  Printf("DebugOFF");  AliLog::SetGlobalDebugLevel(0);}
void don() {  Printf("DebugON");   AliLog::SetGlobalDebugLevel(AliLog::kDebug);}

void geo (                       ) {gGeoManager->GetTopVolume()->Draw("ogl");}
  
void du  (                       ) {r->Dump         (   );}                //utility display 

void hp  (Int_t evt=0            ) {r->HitPrint  (evt);}   //print hits for requested event
void hq  (                       ) {r->HitQA     (   );}   //hits QA plots for all events 
void sp  (Int_t evt=0            ) {r->SdiPrint  (evt);}   //print sdigits for requested event
void sq  (Int_t evt=0            ) {r->SdiPrint  (evt);}   //print sdigits for requested event
void dp  (Int_t evt=0            ) {r->DigPrint  (evt);}   //print digits for requested event
void dq  (                       ) {AliRICHReconstructor::DigQA     (al );}   //digits QA plots for all events
void cp  (Int_t evt=0            ) {r->CluPrint  (evt);                   }   //print clusters for requested event
void cq  (                       ) {AliRICHReconstructor::CluQA     (al );}   //clusters QA plots for all events

void ep  (                       ) {AliRICHTracker::EsdQA(1);            } 
void eq  (                       ) {AliRICHTracker::EsdQA();             }                   
void mp  (Double_t probCut=0.7   ) {AliRICHTracker::MatrixPrint(probCut);}                   


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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tst()
{
  Int_t ch=3;
  TClonesArray hitLst("AliRICHHit");
  TClonesArray sdiLst("AliRICHDigit");
  TObjArray    digLst;                  for(Int_t i=0;i<7;i++) digLst.AddAt(new TClonesArray("AliRICHDigit"),i);
  TClonesArray cluLst("AliRICHCluster");
  
  Int_t iHitCnt=0;  new(hitLst[iHitCnt++]) AliRICHHit(ch,200e-9,kProton,1, 8.40 , 60.14); //c,e,pid,tid,xl,yl
//                    new(hitLst[iHitCnt++]) AliRICHHit(ch,e*1e-9,kProton,2,x,y); //c,e,pid,tid,xl,yl,x,y,z   
  
                                                                     Printf("HIT------HIT---------HIT--------HIT------HIT------HIT");hitLst.Print();Printf("");
  AliRICHv1::Hit2Sdi(&hitLst,&sdiLst);                               Printf("SDI------DIG---------SDI--------SDI------SDI------SDI");sdiLst.Print();Printf("");
  AliRICHDigitizer::Sdi2Dig(&sdiLst,&digLst);                        Printf("DIG------DIG---------DIG--------DIG------DIG------DIG");digLst.Print();Printf("");                   
  AliRICHReconstructor::Dig2Clu((TClonesArray*)digLst[ch],&cluLst);  Printf("CLU------CLU---------CLU--------CLU------CLU------CLU");cluLst.Print();Printf("");                   
  

}


void print()
{
 
  Double_t r2d=TMath::RadToDeg();

  Double_t x=AliRICHDigit::SizeAllX(),y=AliRICHDigit::SizeAllY();
    
  TVector3 c6lt=rp->Lors2Mars(6,0,y);                                              TVector3 c6rt=rp->Lors2Mars(6,x,y);
                                       TVector3 c6ce=rp->Lors2Mars(6,x/2,y/2);
  TVector3 c6lb=rp->Lors2Mars(6,0,0);                                              TVector3 c6rb=rp->Lors2Mars(6,x,0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  TVector3 c5lt=rp->Lors2Mars(5,0,y);                                              TVector3 c5rt=rp->Lors2Mars(5,x,y);
                                       TVector3 c5ce=rp->Lors2Mars(5,x/2,y/2);
  TVector3 c5lb=rp->Lors2Mars(5,0,0);                                              TVector3 c5rb=rp->Lors2Mars(5,x,0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                         
  TVector3 c4lt=rp->Lors2Mars(4,0,y);                                              TVector3 c4rt=rp->Lors2Mars(4,x,y);
                                       TVector3 c4ce=rp->Lors2Mars(4,x/2,y/2);
  TVector3 c4lb=rp->Lors2Mars(4,0,0);                                              TVector3 c4rb=rp->Lors2Mars(4,x,0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TVector3 c3lt=rp->Lors2Mars(3,0,y);                                              TVector3 c3rt=rp->Lors2Mars(3,x,y);
                                       TVector3 c3ce=rp->Lors2Mars(3,x/2,y/2);
  TVector3 c3lb=rp->Lors2Mars(3,0,0);                                              TVector3 c3rb=rp->Lors2Mars(3,x,0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TVector3 c2lt=rp->Lors2Mars(2,0,y);                                              TVector3 c2rt=rp->Lors2Mars(2,x,y);
                                       TVector3 c2ce=rp->Lors2Mars(2,x/2,y/2);
  TVector3 c2lb=rp->Lors2Mars(2,0,0);                                              TVector3 c2rb=rp->Lors2Mars(2,x,0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TVector3 c1lt=rp->Lors2Mars(1,0,y);                                              TVector3 c1rt=rp->Lors2Mars(1,x,y);
                                       TVector3 c1ce=rp->Lors2Mars(1,x/2,y/2);
  TVector3 c1lb=rp->Lors2Mars(1,0,0);                                              TVector3 c1rb=rp->Lors2Mars(1,x,0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TVector3 c0lt=rp->Lors2Mars(0,0,y);                                              TVector3 c0rt=rp->Lors2Mars(0,x,y);
                                       TVector3 c0ce=rp->Lors2Mars(0,x/2,y/2);
  TVector3 c0lb=rp->Lors2Mars(0,0,0);                                              TVector3 c0rb=rp->Lors2Mars(0,x,0);
  
  
  Printf("\n\n\n");                                       
  
  Printf("_______________________________   _______________________________");                                       
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|",c6lt.Mag()          ,c6rt.Mag()        ,c5lt.Mag()        ,c5rt.Mag()                 );
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|",c6lt.Theta()*r2d    ,c6rt.Theta()*r2d  ,c5lt.Theta()*r2d  ,c5rt.Theta()*r2d           );
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|",c6lt.Phi()*r2d      ,c6rt.Phi()*r2d    ,c5lt.Phi()*r2d    ,c5rt.Phi()*r2d             );                                       
  Printf("|                             |   |                             |");                                                                              
  Printf("|         %7.2f             |   |         %7.2f             | Sensitive area  (130.60,126.16)"  ,c6ce.Mag()           ,c5ce.Mag()                );
  Printf("|         %7.2f             |   |         %7.2f             | Lors Center     ( 65.30, 63.08)"  ,c6ce.Theta()*r2d     ,c5ce.Theta()*r2d          );
  Printf("|         %7.2f             |   |         %7.2f             |"                                  ,c6ce.Phi()*r2d       ,c5ce.Phi()*r2d            );
  Printf("|                             |   |                             |");                                                                              
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|",c6lb.Mag()          ,c6rb.Mag()        ,c5lb.Mag()        ,c5rb.Mag()                 );
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|",c6lb.Theta()*r2d    ,c6rb.Theta()*r2d  ,c5lb.Theta()*r2d  ,c5rb.Theta()*r2d           );
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|",c6lb.Phi()*r2d      ,c6rb.Phi()*r2d    ,c5lb.Phi()*r2d    ,c5rb.Phi()*r2d             );                                       
  Printf("-------------------------------   -------------------------------");                                         
  Printf("");
  Printf("_______________________________   _______________________________   _______________________________");                                       
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|   |%7.2f               %7.2f|",c4lt.Mag()      ,c4rt.Mag()      ,c3lt.Mag()      ,c3rt.Mag()      ,c2lt.Mag()       ,c2rt.Mag()      );
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|   |%7.2f               %7.2f|",c4lt.Theta()*r2d,c4rt.Theta()*r2d,c3lt.Theta()*r2d,c3rt.Theta()*r2d,c2lt.Theta()*r2d ,c2rt.Theta()*r2d);
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|   |%7.2f               %7.2f|",c4lt.Phi()*r2d  ,c4rt.Phi()*r2d  ,c3lt.Phi()*r2d  ,c3rt.Phi()*r2d  ,c2lt.Phi()*r2d   ,c2rt.Phi()*r2d  );                                       
  Printf("|                             |   |                             |   |                             |");                                                                              
  Printf("|         %7.2f             |   |         %7.2f             |   |         %7.2f             |"     ,c4ce.Mag()      ,c3ce.Mag()      ,c2ce.Mag());
  Printf("|         %7.2f             |   |         %7.2f             |   |         %7.2f             |"     ,c4ce.Theta()*r2d,c3ce.Theta()*r2d,c2ce.Theta()*r2d);
  Printf("|         %7.2f             |   |         %7.2f             |   |         %7.2f             |"     ,c4ce.Phi()*r2d  ,c3ce.Phi()*r2d  ,c3ce.Phi()*r2d);
  Printf("|                             |   |                             |   |                             |");                                                                              
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|   |%7.2f               %7.2f|",c4lb.Mag()      ,c4rb.Mag()      ,c3lb.Mag()      ,c3rb.Mag()      ,c2lt.Mag()       ,c2rt.Mag()      );
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|   |%7.2f               %7.2f|",c4lb.Theta()*r2d,c4rb.Theta()*r2d,c3lb.Theta()*r2d,c3rb.Theta()*r2d,c2lt.Theta()*r2d ,c2rt.Theta()*r2d);
  Printf("|%7.2f               %7.2f|   |%7.2f               %7.2f|   |%7.2f               %7.2f|",c4lb.Phi()*r2d  ,c4rb.Phi()*r2d  ,c3lb.Phi()*r2d  ,c3rb.Phi()*r2d  ,c2lt.Phi()*r2d   ,c2rt.Phi()*r2d  );                                       
  Printf("-------------------------------   -------------------------------   -------------------------------");                                         
  Printf("");  
  Printf("                                  _______________________________   _______________________________");                                       
  Printf("                                  |%7.2f               %7.2f|   |%7.2f               %7.2f|",c1lt.Mag()      ,c1rt.Mag()      ,c0lt.Mag()      ,c0rt.Mag()      );
  Printf("                                  |%7.2f               %7.2f|   |%7.2f               %7.2f|",c1lt.Theta()*r2d,c1rt.Theta()*r2d,c0lt.Theta()*r2d,c0rt.Theta()*r2d);
  Printf("                                  |%7.2f               %7.2f|   |%7.2f               %7.2f|",c1lt.Phi()*r2d  ,c1rt.Phi()*r2d  ,c0lt.Phi()*r2d  ,c0rt.Phi()*r2d  );                                       
  Printf("                                  |                             |   |                             |");                                                                              
  Printf("                                  |         %7.2f             |   |         %7.2f             |"     ,c1ce.Mag()      ,c0ce.Mag()      );
  Printf("                                  |         %7.2f             |   |         %7.2f             |"     ,c1ce.Theta()*r2d,c0ce.Theta()*r2d);
  Printf("                                  |         %7.2f             |   |         %7.2f             |"     ,c1ce.Phi()*r2d  ,c0ce.Phi()*r2d  );
  Printf("                                  |                             |   |                             |");                                                                              
  Printf("                                  |%7.2f               %7.2f|   |%7.2f               %7.2f|",c1lb.Mag()      ,c1rb.Mag()      ,c0lb.Mag()      ,c0rb.Mag()       );
  Printf("                                  |%7.2f               %7.2f|   |%7.2f               %7.2f|",c1lb.Theta()*r2d,c1rb.Theta()*r2d,c0lb.Theta()*r2d,c0rb.Theta()*r2d );
  Printf("                                  |%7.2f               %7.2f|   |%7.2f               %7.2f|",c1lb.Phi()*r2d  ,c1rb.Phi()*r2d  ,c0lb.Phi()*r2d  ,c0rb.Phi()*r2d   );                                       
  Printf("                                  -------------------------------   -------------------------------");                                       
  
  
  Double_t m[3]; 
  for(int i=0;i<1000;i++){
    Float_t xout=0,xin=gRandom->Rndm()*130.60;
    Float_t yout=0,yin=gRandom->Rndm()*126.16;
    Int_t   c=gRandom->Rndm()*6;
    rp->Lors2Mars(c,xin,yin,m);
    rp->Mars2Lors(c,m,xout,yout);
    if( (xin-xout) != 0) Printf("Problem in X");
    if( (yin-yout) != 0) Printf("Problem in Y");
  }                
}//print()


void rec()
{
  Double_t th=8*TMath::DegToRad();//gRandom->Rndm()*TMath::PiOver4();
  Double_t ph=gRandom->Rndm()*TMath::TwoPi(); 172.51*TMath::DegToRad();
  Double_t  x=gRandom->Rndm()*AliRICHDigit::SizeAllX(); //101.59;
  Double_t  y=gRandom->Rndm()*2*AliRICHDigit::SizePcY();//38.06;
  
  
  Double_t ckovMax=TMath::ACos(1./AliRICHRecon::fkRadIdx),ckovSim;
  Int_t nSim=0;
  while(nSim<3){
    ckovSim=gRandom->Rndm()*ckovMax;//0.6468;
    nSim=20*TMath::Power(TMath::Sin(ckovSim)/TMath::Sin(ckovMax),2); //scale number of photons 
  }
  
  
  TClonesArray *pCluLst=new TClonesArray("AliRICHCluster");
  TPolyMarker  *pMipMap=new TPolyMarker();   pMipMap->SetMarkerStyle(8);  pMipMap->SetMarkerColor(kRed); 
  TPolyMarker  *pPhoMap=new TPolyMarker();   pPhoMap->SetMarkerStyle(4);  pPhoMap->SetMarkerColor(kRed);
  TPolyMarker  *pBkgMap=new TPolyMarker();   pBkgMap->SetMarkerStyle(25); pBkgMap->SetMarkerColor(kRed);
  TPolyLine    *pRing  =new TPolyLine;                                    pRing->SetLineColor(kGreen);     
  TPolyLine    *pOld   =new TPolyLine;                                    pOld->SetLineColor(kBlue);     
  
  Int_t iCluCnt=0; pMipMap->SetPoint(iCluCnt,x,y); new((*pCluLst)[iCluCnt++])   AliRICHCluster(1,x,y,200); //add mip cluster
  
//  for(int j=0;j<30;j++){                                                                                   //add bkg photons  
//    Float_t bkgX=gRandom->Rndm()*AliRICHDigit::SizeAllX();
//    Float_t bkgY=gRandom->Rndm()*AliRICHDigit::SizeAllY();
//    pBkgMap->SetPoint(iCluCnt,bkgX,bkgY); new((*pCluLst)[iCluCnt++]) AliRICHCluster(1,bkgX,bkgY,35);
//  }   

  
  
  
  AliRICHRecon    rec; rec.SetTrack(th,ph,x,y);                                                                
  AliRICHReconOld old; old.SetTrack(th,ph,x,y); 
  
  TVector2 pos;
  for(int i=0;i<nSim;i++){
    rec.TracePhoton(ckovSim,gRandom->Rndm()*2*TMath::Pi(),pos);                                   //add photons 
    if(AliRICHDigit::IsInDead(pos.X(),pos.Y())) continue; 
    pPhoMap->SetPoint(iCluCnt,pos.X(),pos.Y());    new((*pCluLst)[iCluCnt++]) AliRICHCluster(1,pos.X(),pos.Y(),35); 
  }  

  for(int i=0;i<pCluLst->GetEntries();i++){
    AliRICHCluster *pClu=(AliRICHCluster*)pCluLst->At(i);
    Printf("diff phi %f ckov new  %f ckov old %f",rec.FindPhotPhi(pClu->X(),pClu->Y())-old.FindPhotPhi(pClu->X(),pClu->Y()),rec.FindPhotCkov(pClu->X(),pClu->Y()),old.FindPhotCkov(pClu->X(),pClu->Y()));
  }
  
  Int_t nRec=0,nOld=0;
  Double_t ckovRec=rec.CkovAngle(pCluLst,nRec); Double_t err=TMath::Sqrt(rec.CkovSigma2());   
  Double_t ckovOld=old.CkovAngle(pCluLst,nOld);
  
  Printf("---------------- Now reconstructed --------------------");
  
  
  for(int j=0;j<100;j++){rec.TracePhoton(ckovRec,j*0.0628,pos); pRing->SetPoint(j,pos.X(),pos.Y());}  
  for(int j=0;j<100;j++){rec.TracePhoton(ckovOld,j*0.0628,pos); pOld->SetPoint(j,pos.X(),pos.Y());}  
    
  new TCanvas;  AliRICHDigit::DrawPc();  pMipMap->Draw(); pPhoMap->Draw(); pBkgMap->Draw(); pRing->Draw();  pOld->Draw(); 
  
  TLatex txt; txt.SetTextSize(0.03);
  txt.DrawLatex(65,127,Form("#theta=%.4f#pm%.5f with %2i #check{C}"                             ,ckovSim, 0.,nSim             ));
  txt.DrawLatex(65,122,Form("#theta=%.4f#pm%.5f with %2i #check{C} Old=%.4f with %i #check{C}"  ,ckovRec,err,nRec,ckovOld,nOld));
  txt.DrawLatex(0 ,127,Form("#theta=%.2f#circ   #phi=%.2f#circ @(%.2f,%.2f) ",th*TMath::RadToDeg(),ph*TMath::RadToDeg(),x,y));
                   
//  for(int i=0;i<35;i++){
//    Double_t ckov=0.1+i*0.02;
//    Printf("Ckov=%.2f Old=%.3f New=%.3f",ckov,old.FindRingArea(ckov),rec.FindRingArea(ckov));
//  }
  
  
}//rec()

AliRICHRecon rec; AliRICHReconOld old;

void aaa()
{
  rec.SetTrack(10*TMath::DegToRad(),1,30,60); old.SetTrack(10*TMath::DegToRad(),1,30,60);
  TVector2 pos;
  Double_t ckovSim=0.6234;
  rec.TracePhoton(ckovSim,1,pos);
  Printf("ckovSim %f",ckovSim);
  double ckovRec=rec.FindPhotCkov(pos.X(),pos.Y());
  double ckovOld=old.FindPhotCkov(pos.X(),pos.Y());
  Printf("new %f old %f",ckovRec,ckovOld);
}





void ed(Int_t iEvtFrom=-1,Int_t iEvtTo=-1) 
{
// Display digits, reconstructed tracks intersections and RICH rings if available 
  
  
  Int_t iNevt=al->GetNumberOfEvents();
  if(iEvtFrom< 0    ){iEvtFrom=0;iEvtTo=iNevt-1;}
  if(iEvtTo  >=iNevt){           iEvtTo=iNevt-1;}
  
  rl->LoadDigits();
  
  TFile *pEsdFl=TFile::Open("AliESDs.root");     if(!pEsdFl || !pFile->IsOpen()) return;//open AliESDs.root                                                                    
  TTree *pEsdTr=(TTree*) pFile->Get("esdTree");  if(!pEsdTr)                     return;//get ESD tree
                                                                 
  AliESD *pEsd=new AliESD;  pTree->SetBranchAddress("ESD", &pEsd);
  
  
  
  
  TPolyMarker  *pDigMap[7]; //digits map
  TPolyMarker  *pTrkMap[7]; //TRKxPC intersection map
  
  for(Int_t i=0;i<7;i++){
    pDigMap[i]=new TPolyMarker();       pDigMap[i]->SetMarkerStyle(25);   pDigMap[i]->SetMarkerSize(0.5); pDigMap[i]->SetMarkerColor(kGreen); 
    pTrkMap[i]=new TPolyMarker();       pTrkMap[i]->SetMarkerStyle(4);    pTrkMap[i]->SetMarkerSize(0.5); pTrkMap[i]->SetMarkerColor(kRed); 
  }

  TLatex t;
  TCanvas *pC = new TCanvas("RICHDisplay","RICH Display",0,0,1226,900);  pC->Divide(3,3);
  
  for(Int_t iEvt=iEvtFrom;iEvt<=iEvtTo;iEvt++) {                //events loop
    pC->cd(3);  t.DrawText(0.2,0.4,Form("Event %i",iEvt));        //print current event number
        
    pEsdTr->GetEntry(iEvt);                              //get ESD for this event   
    for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//ESD tracks loop
      AliESDtrack *pTrk = pEsd->GetTrack(iTrk);             //
      Float_t th,ph,x,y;
      pTrk->GetRICHtrk(th,ph,x,y);
      
    }//ESD tracks loop
    
    al->GetEvent(iEvt);   rl->TreeD()->GetEntry(0); //get digits list
    for(Int_t iCh=0;iCh<7;iCh++) {//chambers loop
      for(Int_t iDig=0;iDig < r->DigLst(iCh)->GetEntries();iDig++) {      //digits loop
        AliRICHDigit *pDig = (AliRICHDigit*)r->DigLst(iCh)->At(iDig);     
        pDigMap[iCh]->SetPoint(iDig,pDig->LorsX(),pDig->LorsY());
      }                                                             //digits loop

      
      if(iCh==6) pC->cd(1); if(iCh==5) pC->cd(2);
      if(iCh==4) pC->cd(4); if(iCh==3) pC->cd(5); if(iCh==2) pC->cd(6);
                            if(iCh==1) pC->cd(8); if(iCh==0) pC->cd(9);
      
      AliRICHDigit::DrawPc();  pDigMap[iCh]->Draw();
    }//chambers loop
    pC->Update();
    pC->Modified();

    if(iEvt<iEvtTo) {gPad->WaitPrimitive();pC->Clear();}
    
    
    
  }//events loop
  delete pEsf;  pEsdFl->Close();//close AliESDs.root
  rl->UnloadDigits();
}//Display()


void AliRICH::HitQA(Double_t cut,Double_t cutele,Double_t cutR)
{
// Provides a set of control plots intended primarily for charged particle flux analisys
// Arguments: cut (GeV)    - cut on momentum of any charged particles but electrons, 
//            cetele (GeV) - the same for electrons-positrons
//            cutR (cm)    - cut on production vertex radius (cylindrical system)        
  gBenchmark->Start("HitsAna");
  
  Double_t cutPantiproton    =cut;
  Double_t cutPkaonminus     =cut;
  Double_t cutPpionminus     =cut;
  Double_t cutPmuonminus     =cut;
  Double_t cutPpositron      =cutele;
                    
  Double_t cutPelectron      =cutele;
  Double_t cutPmuonplus      =cut;
  Double_t cutPpionplus      =cut;
  Double_t cutPkaonplus      =cut;
  Double_t cutPproton        =cut;
                       

  TH2F *pEleHitRZ    =new TH2F("EleHitRZ"    ,Form("e^{+} e^{-} hit %s;z[cm];R[cm]" ,GetName())     , 400,-300,300 ,400,-500,500);   //R-z plot 0cm<R<550cm -300cm<z<300cm  
  TH2F *pEleHitRP    =new TH2F("EleHitRP"    ,Form("e^{+} e^{-} hit %s;p[GeV];R[cm]",GetName())     ,1000,-1  ,1   ,400,   0,550);   //R-p plot 0cm<R<550cm -1GeV<p<1GeV 
  TH1F *pEleAllP     =new TH1F("EleAllP"     ,     "e^{+} e^{-} all;p[GeV]"                         ,1000,-1  ,1                );  
  TH1F *pEleHitP     =new TH1F("EleHitP"     ,Form("e^{+} e^{-} hit %s;p[GeV]"      ,GetName())     ,1000,-1  ,1                );   
  TH1F *pMuoHitP     =new TH1F("MuoHitP"     ,Form("#mu^{-} #mu^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pPioHitP     =new TH1F("PioHitP"     ,Form("#pi^{-} #pi^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pKaoHitP     =new TH1F("KaoHitP"     ,Form("K^{-} K^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pProHitP     =new TH1F("ProHitP"     ,Form("p^{-} p^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH2F *pFlux        =new TH2F("flux"        ,Form("%s flux with Rvertex<%.1fcm"    ,GetName(),cutR),10  ,-5  ,5   , 10,0    ,10); //special text hist
  TH2F *pVertex      =new TH2F("vertex"      ,Form("%s 2D vertex of RICH hit;x;y"   ,GetName())     ,120 ,0   ,600 ,120,0    ,600); //special text hist
  TH1F *pRho         =new TH1F("rho"         ,Form("%s r of RICH hit"               ,GetName())     ,600 ,0   ,600); //special text hist
  pFlux->SetStats(0);
  pFlux->GetXaxis()->SetBinLabel(1 ,Form("p^{-}>%.3fGeV/c"   ,cutPantiproton));        
  pFlux->GetXaxis()->SetBinLabel(2 ,Form("K^{-}>%.3fGeV/c"   ,cutPkaonminus ));        
  pFlux->GetXaxis()->SetBinLabel(3 ,Form("#pi^{-}>%.3fGeV/c" ,cutPpionminus ));      
  pFlux->GetXaxis()->SetBinLabel(4 ,Form("#mu^{-}>%.3fGeV/c" ,cutPmuonminus ));      
  pFlux->GetXaxis()->SetBinLabel(5 ,Form("e^{+}>%.3fGeV/c"   ,cutPpositron  ));        
  
  pFlux->GetXaxis()->SetBinLabel(6 ,Form("e^{-}>%.3fGeV/c"   ,cutPelectron  ));        
  pFlux->GetXaxis()->SetBinLabel(7 ,Form("#mu^{+}>%.3fGeV/c" ,cutPmuonplus  ));      
  pFlux->GetXaxis()->SetBinLabel(8 ,Form("#pi^{+}>%.3fGeV/c" ,cutPpionplus  ));      
  pFlux->GetXaxis()->SetBinLabel(9 ,Form("K^{+}>%.3fGeV/c"   ,cutPkaonplus  ));        
  pFlux->GetXaxis()->SetBinLabel(10,Form("p^{+}>%.3fGeV/c"   ,cutPproton    ));        
  
  pFlux->GetYaxis()->SetBinLabel(1,"sum");  
  pFlux->GetYaxis()->SetBinLabel(2,"ch1");  
  pFlux->GetYaxis()->SetBinLabel(3,"ch2");  
  pFlux->GetYaxis()->SetBinLabel(4,"ch3");  
  pFlux->GetYaxis()->SetBinLabel(5,"ch4");  
  pFlux->GetYaxis()->SetBinLabel(6,"ch5");  
  pFlux->GetYaxis()->SetBinLabel(7,"ch6");  
  pFlux->GetYaxis()->SetBinLabel(8,"ch7");  
  pFlux->GetYaxis()->SetBinLabel(9,"prim"); 
  pFlux->GetYaxis()->SetBinLabel(10,"tot");  
  
//end of hists definition
   
  Int_t iNevents=fLoader->GetRunLoader()->GetAliRun()->GetEventsPerRun(),iCntPrimParts=0,iCntTotParts=0;
//load all needed trees   
  fLoader->LoadHits(); 
  fLoader->GetRunLoader()->LoadHeader(); 
  fLoader->GetRunLoader()->LoadKinematics();  
  
  for(Int_t iEvtN=0;iEvtN < iNevents;iEvtN++){//events loop
    fLoader->GetRunLoader()->GetEvent(iEvtN);
    AliInfo(Form(" %i event processes",fLoader->GetRunLoader()->GetEventNumber()));
    AliStack *pStack= fLoader->GetRunLoader()->Stack(); 
    
    for(Int_t iParticleN=0;iParticleN<pStack->GetNtrack();iParticleN++){//stack loop
      TParticle *pPart=pStack->Particle(iParticleN);

      if(iParticleN%10000==0) AliInfo(Form(" %i particles read",iParticleN));
    
      switch(pPart->GetPdgCode()){
        case kProtonBar: pFlux->Fill(-4.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-4.5,8); break;
        case kKMinus:    pFlux->Fill(-3.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-3.5,8); break;
        case kPiMinus:   pFlux->Fill(-2.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-2.5,8); break;
        case kMuonMinus: pFlux->Fill(-1.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-1.5,8); break;
        case kPositron:  pFlux->Fill(-0.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-0.5,8); pEleAllP->Fill(-pPart->P()); break;
      
        case kElectron:  pFlux->Fill( 0.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 0.5,8); pEleAllP->Fill( pPart->P()); break;      
        case kMuonPlus:  pFlux->Fill( 1.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 1.5,8); break;      
        case kPiPlus:    pFlux->Fill( 2.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 2.5,8); break;      
        case kKPlus:     pFlux->Fill( 3.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 3.5,8); break;      
        case kProton:    pFlux->Fill( 4.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 4.5,8); break;            
      }//switch
    }//stack loop
//now hits analiser        
    for(Int_t iEntryN=0;iEntryN < fLoader->TreeH()->GetEntries();iEntryN++){//TreeH loop
      fLoader->TreeH()->GetEntry(iEntryN);                                  //get current entry (prim)                
      for(Int_t iHitN=0;iHitN < Hits()->GetEntries();iHitN++){//hits loop
        AliRICHHit *pHit = (AliRICHHit*)Hits()->At(iHitN);            //get current hit
        TParticle  *pPart=pStack->Particle(pHit->GetTrack());      //get stack particle which produced the current hit
        
        if(pPart->GetPDG()->Charge()!=0&&pPart->Rho()>0.1) pVertex->Fill(pPart->Vx(),pPart->Vy()); //safe margin for sec.
        if(pPart->GetPDG()->Charge()!=0) pRho->Fill(pPart->Rho()); //safe margin for sec.
        if(pPart->R()>cutR) continue;                                   //cut on production radius (cylindrical system) 
      
        switch(pPart->GetPdgCode()){
          case kProtonBar: if(pPart->P()>cutPantiproton) {pProHitP->Fill(-pPart->P()); pFlux->Fill(-4.5,pHit->Ch());}break;
          case kKMinus   : if(pPart->P()>cutPkaonminus)  {pKaoHitP->Fill(-pPart->P()); pFlux->Fill(-3.5,pHit->Ch());}break;
          case kPiMinus  : if(pPart->P()>cutPpionminus)  {pPioHitP->Fill(-pPart->P()); pFlux->Fill(-2.5,pHit->Ch());}break;
          case kMuonMinus: if(pPart->P()>cutPmuonminus)  {pMuoHitP->Fill(-pPart->P()); pFlux->Fill(-1.5,pHit->Ch());}break;        
          case kPositron : if(pPart->P()>cutPpositron)   {pEleHitP->Fill(-pPart->P()); pFlux->Fill(-0.5,pHit->Ch()); 
               pEleHitRP->Fill(-pPart->P(),pPart->R());  pEleHitRZ->Fill(pPart->Vz(),pPart->R()); }break;
          
          case kElectron : if(pPart->P()>cutPelectron)   {pEleHitP->Fill( pPart->P()); pFlux->Fill( 0.5,pHit->Ch()); 
               pEleHitRP->Fill( pPart->P(),pPart->R());  pEleHitRZ->Fill(pPart->Vz(),pPart->R()); }break;
          case kMuonPlus : if(pPart->P()>cutPmuonplus)   {pMuoHitP->Fill( pPart->P()); pFlux->Fill( 1.5,pHit->Ch());}break;                     
          case kPiPlus   : if(pPart->P()>cutPpionplus)   {pPioHitP->Fill( pPart->P()); pFlux->Fill( 2.5,pHit->Ch());}break;           
          case kKPlus    : if(pPart->P()>cutPkaonplus)   {pKaoHitP->Fill( pPart->P()); pFlux->Fill( 3.5,pHit->Ch());}break;           
          case kProton   : if(pPart->P()>cutPproton)     {pProHitP->Fill( pPart->P()); pFlux->Fill( 4.5,pHit->Ch());}break;
        }
      }//hits loop      
    }//TreeH loop
    iCntPrimParts +=pStack->GetNprimary();
    iCntTotParts  +=pStack->GetNtrack();
  }//events loop                        
//unload all loaded staff  
  fLoader->UnloadHits();  
  fLoader->GetRunLoader()->UnloadHeader(); 
  fLoader->GetRunLoader()->UnloadKinematics();  
//Calculater some sums
  Stat_t sum=0;
//sum row, sum over rows  
  for(Int_t i=1;i<=pFlux->GetNbinsX();i++){
    sum=0; for(Int_t j=2;j<=8;j++)    sum+=pFlux->GetBinContent(i,j);    
    pFlux->SetBinContent(i,1,sum);
  }
    
//display everything  
  new TCanvas("canvas1",Form("Events %i Nprims=%i Nparticles=%i",iNevents,iCntPrimParts,iCntTotParts),1000,900); pFlux->Draw("text");  gPad->SetGrid();  
//total prims and particles
  TLatex latex; latex.SetTextSize(0.02);
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,10);    latex.DrawLatex(5.1,9.5,Form("%.0f",sum));
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,9);     latex.DrawLatex(5.1,8.5,Form("%.0f",sum));
  for(Int_t iCh=0;iCh<7;iCh++) {
    sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,iCh+2);latex.DrawLatex(5.1,iCh+1.5,Form("%.0f",sum));
  }  
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,1);    latex.DrawLatex(5.1,0.5,Form("%.0f",sum));
  
  new TCanvas("cEleAllP"   ,"e" ,200,100); pEleAllP->Draw();
  new TCanvas("cEleHitRP"  ,"e" ,200,100); pEleHitRP->Draw();
  new TCanvas("cEleHitRZ"  ,"e" ,200,100); pEleHitRZ->Draw();
  new TCanvas("cEleHitP"   ,"e" ,200,100); pEleHitP->Draw();
  new TCanvas("cMuoHitP"   ,"mu",200,100); pMuoHitP->Draw();
  new TCanvas("cPioHitP"   ,"pi",200,100); pPioHitP->Draw();
  new TCanvas("cKaoHitP"   ,"K" ,200,100); pKaoHitP->Draw();
  new TCanvas("cProHitP"   ,"p" ,200,100); pProHitP->Draw();
  new TCanvas("cVertex"    ,"2d vertex" ,200,100); pVertex->Draw();
  new TCanvas("cRho"    ,"Rho of sec" ,200,100); pRho->Draw();
  
  gBenchmark->Show("HitsPlots");
}//HitQA()





