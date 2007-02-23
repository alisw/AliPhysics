AliRun     *a; AliRunLoader *al;   TGeoManager *g; //globals for easy manual manipulations
AliHMPID   *h; AliLoader    *hl; AliHMPIDParam *hp;
Bool_t isGeomType=kFALSE;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GetParam()
{
  isGeomType=!isGeomType;
  if(g) delete g;  if(hp) delete hp; //delete current TGeoManager and AliHMPIDParam
  if(isGeomType) g=TGeoManager::Import("geometry.root");
  else           g=TGeoManager::Import("misaligned_geometry.root");
  hp=AliHMPIDParam::Instance();
}//GetParam()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hmenu()
{   
  TString status="Status: ";
  if(gSystem->IsFileInIncludePath("galice.root")){
    status+="galice.root found";
    al=AliRunLoader::Open();                                                //try to open galice.root from current dir 
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    al->LoadgAlice(); a=al->GetAliRun();                                    //take new AliRun object from galice.root   
    hl=al->GetDetectorLoader("HMPID");  h=(AliHMPID*)a->GetDetector("HMPID");  //get HMPID object from galice.root
    
    status+=Form(" with %i event(s)",al->GetNumberOfEvents()); status+=(h)? " with HMPID": " without HMPID";
  }else  
    status+="No galice.root";
  
  GetParam();
  
  TControlBar *pMenu = new TControlBar("horizontal",status.Data(),0,0);
    pMenu->AddButton("                     ","","");
    pMenu->AddButton("       General       ","General()"  ,"general items which do not depend on any files");
    pMenu->AddButton("                     ",""           ,"");
    pMenu->AddButton("       Sim data      ","SimData()"  ,"items which expect to have simulated files"    );
    pMenu->AddButton("                     ",""           ,"");
    pMenu->AddButton("       Raw data      ","RawData()"  ,"items which expect to have raw files"          );
    pMenu->AddButton("                     ","print()"    ,"");
    pMenu->AddButton("         Test        ","Test()"     ,"all test utilities");
    pMenu->AddButton("                     ","GetParam()" ,"");
    pMenu->AddButton("         Quit        ",".q"         ,"close session"                                 );
  pMenu->Show();
}//Menu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void General()
{         
  TControlBar *pMenu = new TControlBar("vertical","General purpose",100,50);  
    pMenu->AddButton("Debug ON","don();"                    ,"Switch debug on-off"                        );   
    pMenu->AddButton("Debug OFF","doff();"                   ,"Switch debug on-off"                        );   
    pMenu->AddButton("Geo GUI","geo();"                    ,"Shows geometry"                             ); 
    pMenu->AddButton("Browser","new TBrowser;"             ,"Start ROOT TBrowser"                        );
  pMenu->Show();  
}//General()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimData()
{
  TControlBar *pMenu = new TControlBar("vertical","Sim data",310,50);  
    pMenu->AddButton("Display ","hed();"    ,"Display Fast");
    pMenu->AddButton("HITS QA"           ,"hqa()"     ,"QA plots for hits: hqa()");
    pMenu->AddButton("Print stack"       ,"stack();"  ,"To print hits:     hp(evt)");
    pMenu->AddButton("Print hits"        ,"hp();"     ,"To print hits:     hp(evt)");
    pMenu->AddButton("Print sdigits"     ,"sp();"     ,"To print sdigits:  sp(evt)");
    pMenu->AddButton("Print digits"      ,"dp();"     ,"To print digits:   dp(evt)");
    pMenu->AddButton("Print clusters"    ,"cp();"     ,"To print clusters: cp(evt)");
  pMenu->Show();         
}//SimData()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RawData()
{
  TControlBar *pMenu = new TControlBar("vertical","Raw data",580,50);  
    pMenu->AddButton("ESD print"                       ,"ep();"                  ,"To print ESD info: ep()"         );  
    pMenu->AddButton("ESD QA"                          ,"eq();"                  ,"To draw ESD hists: eq()"         );  
    pMenu->AddButton("Clusters print"                  ,"cp();"                  ,"To print clusters: cp()"         );  
    pMenu->AddButton("Clusters QA"                     ,"cq();"                  ,"To draw clusters hists: cq()"    );  
    pMenu->AddButton("Print Matrix"                    ,"mp();"                  ,"To print prob matrix: mp()"      );  
    pMenu->AddButton("Print occupancy"                 ,"r->OccupancyPrint(-1);" ,"To print occupancy"              );  
    pMenu->AddButton("Print event summary  "           ,"r->SummaryOfEvent();"   ,"To print a summary of the event" );  
  pMenu->Show();         
}//RawData()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Test()
{         
  TControlBar *pMenu = new TControlBar("vertical","Test",820,50);  
    pMenu->AddButton("TEST Display "      ,"sed();"    ,"Display Fast");
    pMenu->AddButton("Hits->Digits"       ,"thd();"                    ,"test hits->sdigits->digits"                 );   
    pMenu->AddButton("Segmentation"       ,"ts()"                      ,"test segmentation methods"                  );
    pMenu->AddButton("Test response"      ,"AliHMPIDParam::TestResp();","Test AliHMPIDParam response methods"         );
    pMenu->AddButton("Print map"          ,"PrintMap();"               ,"Test AliHMPIDParam transformation methods"   );
    pMenu->AddButton("Test Recon"         ,"rec();"                    ,"Test AliHMPIDRecon"                          );
  pMenu->Show();  
}//Test()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void doff(){  Printf("DebugOFF");  AliLog::SetGlobalDebugLevel(0);}
void don() {  Printf("DebugON");   AliLog::SetGlobalDebugLevel(AliLog::kDebug);}

void geo (                       ) {gGeoManager->GetTopVolume()->Draw("ogl");}
  
void du  (                       ) {h->Dump         (   );}                //utility display 

void hp  (Int_t evt=0            ) {h->HitPrint  (evt);}   //print hits for requested event
void hq  (                       ) {h->HitQA     (   );}   //hits QA plots for all events 
void sp  (Int_t evt=0            ) {h->SdiPrint  (evt);}   //print sdigits for requested event
void sq  (Int_t evt=0            ) {h->SdiPrint  (evt);}   //print sdigits for requested event
void dp  (Int_t evt=0            ) {h->DigPrint  (evt);}   //print digits for requested event
void dq  (                       ) {AliHMPIDReconstructor::DigQA     (al );}   //digits QA plots for all events
void cp  (Int_t evt=0            ) {h->CluPrint  (evt);                   }   //print clusters for requested event
void cq  (                       ) {AliHMPIDReconstructor::CluQA     (al );}   //clusters QA plots for all events

void ep  (                       ) {AliHMPIDTracker::EsdQA(1);            } 
void eq  (                       ) {AliHMPIDTracker::EsdQA();             }                   
void mp  (Double_t probCut=0.7   ) {AliHMPIDTracker::MatrixPrint(probCut);}                   


void stack(                     )   {AliHMPIDParam::Stack();}    
void tid  (Int_t tid,Int_t evt=0)   {AliHMPIDParam::Stack(evt,tid);} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintMap()
{
 
  Double_t r2d=TMath::RadToDeg();

  Double_t x=AliHMPIDDigit::SizeAllX(),y=AliHMPIDDigit::SizeAllY();
    
  Printf("\n\n\n");                                       
  
  for(int ch=6;ch>=0;ch--){
    AliHMPIDDigit dL(ch,0,1,0,0),dR(ch,0,1,67,0);
    TVector3 lt=rp->Lors2Mars(ch,0,y);                                              TVector3 rt=rp->Lors2Mars(ch,x,y);
                                       TVector3 ce=rp->Lors2Mars(ch,x/2,y/2);
    TVector3 lb=rp->Lors2Mars(ch,0,0);                                              TVector3 rb=rp->Lors2Mars(ch,x,0);
    
    Printf(" ____________________________");                                       
    Printf("|%6.2fcm            %6.2fcm|"         ,lt.Mag()                             , rt.Mag()       );
    Printf("|%6.2fde            %6.2fde|"         ,lt.Theta()*r2d                       , rt.Theta()*r2d );
    Printf("|%6.2fde            %6.2fde|"         ,lt.Phi()*r2d                         , rt.Phi()*r2d   );                                       
    Printf("|                            |"                                                       );
    Printf("|DDL %2i    %7.2fcm   DDL %2i|"       ,dL.DdlIdx()    ,  ce.Mag()           , dR.DdlIdx()    );
    Printf("| 0x%x    %7.2fdeg   0x%x|"           ,dL.DdlId()     ,  ce.Theta()*r2d     , dR.DdlId()     );
    Printf("|          %7.2fdeg        |"                         ,  ce.Phi()*r2d                        );
    Printf("|                            |");                                                                              
    Printf("|%6.2fcm            %6.2fcm|"         ,lb.Mag()                             , rb.Mag()       );
    Printf("|%6.2fde            %6.2fde|"         ,lb.Theta()*r2d                       , rb.Theta()*r2d );
    Printf("|%6.2fde     Ch%i    %6.2fde|"        ,lb.Phi()*r2d   ,  ch                 , rb.Phi()*r2d   );                                       
    Printf(" ----------------------------");                                         
  }
  
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
  
  Int_t ddl,r,d,a,ch,raw,pc,px,py; AliHMPIDDigit dd;
  
  ddl=0;raw=0x2214000;r= 8;d=8;a=20;
  ddl=1;raw=0x2214000;r= 8;d=8;a=20;
  
  
  ddl=2;raw=0x08d6000;r= 2;d=3;a=22;
  ddl=3;raw=0x08d6000;r= 2;d=3;a=22;
  
  
  ddl=6;raw=0x592e000;r=22;d=4;a=46;ch=3;pc=4;px=55;py=5;dd.Raw(ddl,raw); 
  Printf("(ch=%i,pc=%i,x=%2i,y=%2i) ddl=%i raw=0x%h (r=%2i,d=%2i,a=%2i)",
           ch,   pc,  px,   py,     ddl,   raw,      r,    d,    a); dd.Print(); 
  ddl=7;raw=0x592e000;r=22;d=4;a=46;ch=3;pc=1;
  
  
}//PrintMap()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HitQA(Double_t cut,Double_t cutele,Double_t cutR)
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
                       

  TH2F *pEleHitRZ    =new TH2F("EleHitRZ"    ,Form("e^{+} e^{-} hit %s;z[cm];R[cm]" ,GetName())     , 400,-300,300 ,400,-500,500);   //R-z
  TH2F *pEleHitRP    =new TH2F("EleHitRP"    ,Form("e^{+} e^{-} hit %s;p[GeV];R[cm]",GetName())     ,1000,-1  ,1   ,400,   0,550);   //R-p
  TH1F *pEleAllP     =new TH1F("EleAllP"     ,     "e^{+} e^{-} all;p[GeV]"                         ,1000,-1  ,1                );  
  TH1F *pEleHitP     =new TH1F("EleHitP"     ,Form("e^{+} e^{-} hit %s;p[GeV]"      ,GetName())     ,1000,-1  ,1                );   
  TH1F *pMuoHitP     =new TH1F("MuoHitP"     ,Form("#mu^{-} #mu^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pPioHitP     =new TH1F("PioHitP"     ,Form("#pi^{-} #pi^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pKaoHitP     =new TH1F("KaoHitP"     ,Form("K^{-} K^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pProHitP     =new TH1F("ProHitP"     ,Form("p^{-} p^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH2F *pFlux        =new TH2F("flux"        ,Form("%s flux with Rvertex<%.1fcm"    ,GetName(),cutR),10  ,-5  ,5   , 10,0    ,10); //special text hist
  TH2F *pVertex      =new TH2F("vertex"      ,Form("%s 2D vertex of HMPID hit;x;y"   ,GetName())     ,120 ,0   ,600 ,120,0    ,600); //special text hist
  TH1F *pRho         =new TH1F("rho"         ,Form("%s r of HMPID hit"               ,GetName())     ,600 ,0   ,600); //special text hist
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
        AliHMPIDHit *pHit = (AliHMPIDHit*)Hits()->At(iHitN);            //get current hit
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CluQA(AliRunLoader *pAL)
{
// Quality assesment plots for clusters. 
// This methode takes list of digits and form list of clusters again in order to 
// calculate cluster shape and cluster particle mixture    
  AliLoader *pRL=pAL->GetDetectorLoader("HMPID");  AliHMPID *pRich=(AliHMPID*)pAL->GetAliRun()->GetDetector("HMPID");//get pointers for HMPID and HMPID loader
  Int_t iNevt=pAL->GetNumberOfEvents();  if(iNevt==0)             {AliInfoClass("No events");return;}   
                                         if(pRL->LoadDigits())    {AliInfoClass("No digits file");return;}
                                            pAL->LoadHeader();
                                            pAL->LoadKinematics(); 
//                                            AliStack *pStack=pAL->Stack();
  TH1::AddDirectory(kFALSE);
  
        
  TH1F*    pQ=new TH1F("HmpAllQ"  ,"Charge All"           ,4000 ,0  ,4000);// Q hists
  TH1F* pCerQ=new TH1F("HmpCerQ"  ,"Charge Ckov"          ,4000 ,0  ,4000);
  TH1F* pMipQ=new TH1F("HmpMipQ"  ,"Charge MIP"           ,4000 ,0  ,4000);
  
  TH1F*    pS=new TH1F("HmpCluSize"    ,"Cluster size;size"         ,100  ,0  ,100 );// size hists
  TH1F* pCerS=new TH1F("HmpCluCerSize" ,"Ckov size;size"            ,100  ,0  ,100 );
  TH1F* pMipS=new TH1F("HmpCluMipSize" ,"MIP size;size"             ,100  ,0  ,100 );
  
  TH2F*    pM=new TH2F("HmpCluMap"     ,"Cluster map;x [cm];y [cm]" ,1000 ,0  ,AliHMPIDDigit::SizeAllX(),1000,0,AliHMPIDDigit::SizeAllY()); // maps
  TH2F* pMipM=new TH2F("HmpCluMipMap"  ,"MIP map;x [cm];y [cm]"     ,1000 ,0  ,AliHMPIDDigit::SizeAllX(),1000,0,AliHMPIDDigit::SizeAllY());
  TH2F* pCerM=new TH2F("HmpCluCerMap"  ,"Ckov map;x [cm];y [cm]"    ,1000 ,0  ,AliHMPIDDigit::SizeAllX(),1000,0,AliHMPIDDigit::SizeAllY());
 
  
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){
    pAL->GetEvent(iEvt);               
    pRL->TreeD()->GetEntry(0); 
    TClonesArray *pCluLst=new TClonesArray("AliHMPIDCluster");//tmp list of clusters for this event
    
    for(Int_t iCh=0;iCh<7;iCh++) AliHMPIDReconstructor::Dig2Clu(pRich->DigLst(iCh),pCluLst,kFALSE);//cluster finder for all chamber if any digits present
    
    for(Int_t iClu=0;iClu<pCluLst->GetEntriesFast();iClu++){
      AliHMPIDCluster *pClu = (AliHMPIDCluster*)pCluLst->At(iClu);
      Int_t cfm=0; for(Int_t iDig=0;iDig<pClu->Size();iDig++)  cfm+=pClu->Dig(iDig)->Ch(); //collect ckov-fee-mip structure of current cluster ?????
      Int_t iNckov=cfm/1000000;      Int_t iNfee =cfm%1000000/1000;      Int_t iNmip =cfm%1000000%1000; 

                                             pQ   ->Fill(pClu->Q()) ; pS   ->Fill(pClu->Size()) ; pM    ->Fill(pClu->X(),pClu->Y()); //all clusters  
      if(iNckov!=0 && iNfee==0 && iNmip==0) {pCerQ->Fill(pClu->Q()) ; pCerS->Fill(pClu->Size()) ; pCerM ->Fill(pClu->X(),pClu->Y());}//ckov only cluster
      if(iNckov==0 && iNfee==0 && iNmip!=0) {pMipQ->Fill(pClu->Q()) ; pMipS->Fill(pClu->Size()) ; pMipM ->Fill(pClu->X(),pClu->Y());}//mip only cluster
                                       
    }//clusters loop   
    pCluLst->Clear();delete pCluLst;
  }//events loop
  
  pRL->UnloadDigits(); pAL->UnloadKinematics(); pAL->UnloadHeader();
  TCanvas *pC=new TCanvas("RichCluQA",Form("QA for cluster from %i events",iNevt),1000,900); pC->Divide(3,3);
  pC->cd(1);    pM->Draw();          pC->cd(2);    pQ->Draw();       pC->cd(3);    pS->Draw();        
  pC->cd(4); pMipM->Draw();          pC->cd(5); pMipQ->Draw();       pC->cd(6); pMipS->Draw();        
  pC->cd(7); pCerM->Draw();          pC->cd(8); pCerQ->Draw();       pC->cd(9); pCerS->Draw();        
}//CluQA()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OccupancyPrint(Int_t iEvtNreq)
{
//prints occupancy for each chamber in a given event
  Int_t iEvtNmin,iEvtNmax;
  if(iEvtNreq==-1){
    iEvtNmin=0;
    iEvtNmax=gAlice->GetEventsPerRun();
  } else { 
    iEvtNmin=iEvtNreq;iEvtNmax=iEvtNreq+1;
  }
    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;    
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;    
  
//  Info("Occupancy","for event %i",iEvtN);
  if(GetLoader()->LoadHits()) return;
  if(GetLoader()->LoadDigits()) return;

  
  for(Int_t iEvtN=iEvtNmin;iEvtN<iEvtNmax;iEvtN++){    
    Int_t nDigCh[7]={0,0,0,0,0,0,0};  
    Int_t iChHits[7]={0,0,0,0,0,0,0};
    Int_t nPrim[7]={0,0,0,0,0,0,0};
    Int_t nSec[7]={0,0,0,0,0,0,0};
    AliInfo(Form("events processed %i",iEvtN));
    if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
    AliStack *pStack = GetLoader()->GetRunLoader()->Stack();
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);      
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
        AliHMPIDHit *pHit = (AliHMPIDHit*)Hits()->At(iHitN);
        if(pHit->E()>0){
          iChHits[pHit->Ch()]++;
          if(pStack->Particle(pHit->GetTrack())->Rho()<0.01) nPrim[pHit->Ch()]++;else nSec[pHit->Ch()]++;
        }
      }
    }
    
    GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iCh=0;iCh<7;iCh++){
      for(Int_t iDig=0;iDig<DigLst(iCh)->GetEntries();iDig++){
        AliHMPIDDigit *pDig=(AliHMPIDDigit*)DigLst(iCh)->At(iDig);
        nDigCh[pDig->Ch()]++;
      }
      Printf("Occupancy for chamber %i = %4.2f %% and charged prim tracks %i and sec. tracks %i with total %i",
        iCh,Float_t(nDigCh[iCh])*100/AliHMPIDDigit::kPadAll,nPrim[iCh],nSec[iCh],iChHits[iCh]);
    }      
    
    
  }//events loop
  GetLoader()->UnloadHits();
  GetLoader()->UnloadDigits();
  GetLoader()->GetRunLoader()->UnloadHeader();    
  GetLoader()->GetRunLoader()->UnloadKinematics();    
}//OccupancyPrint()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPID::SummaryOfEvent(Int_t iEvtN) const
{
//prints a summary for a given event
  AliInfo(Form("Summary of event %i",iEvtN));
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;
  AliStack *pStack=GetLoader()->GetRunLoader()->Stack();
  
  AliGenEventHeader* pGenHeader =  gAlice->GetHeader()->GenEventHeader();
  if(pGenHeader->InheritsFrom("AliGenHijingEventHeader")) {
    AliInfo(Form(" Hijing event with impact parameter b = %.2f (fm)",((AliGenHijingEventHeader*) pGenHeader)->ImpactParameter()));
  }
  Int_t nChargedPrimaries=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) {
    TParticle *pParticle = pStack->Particle(i);
    if(pParticle->IsPrimary()&&pParticle->GetPDG()->Charge()!=0) nChargedPrimaries++;
    }
  AliInfo(Form("Total number of         primaries %i",pStack->GetNprimary()));
  AliInfo(Form("Total number of charged primaries %i",nChargedPrimaries));
  AliInfo(Form("Total n. of tracks in stack(+sec) %i",pStack->GetNtrack()));
  GetLoader()->GetRunLoader()->UnloadHeader();
  GetLoader()->GetRunLoader()->UnloadKinematics();
}//SummaryOfEvent()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawZoom()
{
// Show info about current cursur position in status bar of the canvas
// Arguments: none
//   Returns: none     
  TCanvas *pC=(TCanvas*)gPad; 
  TRootCanvas *pRC= (TRootCanvas*)pC->GetCanvasImp();
  TGStatusBar *pBar=pRC->GetStatusBar();
  pBar->SetParts(5);
  Float_t x=gPad->AbsPixeltoX(gPad->GetEventX());
  Float_t y=gPad->AbsPixeltoY(gPad->GetEventY());
  AliHMPIDDigit dig;dig.Manual1(1,x,y); UInt_t w32=0; 
  if(IsInDead(x,y))
    pBar->SetText("Out of sensitive area",4);    
  else{
    Int_t ddl=dig.Raw(w32);
    pBar->SetText(Form("(p%i,x%i,y%i) ddl=%i 0x%x (r%i,d%i,a%i) (%.2f,%.2f)",
        dig.Pc(),dig.PadPcX(),dig.PadPcY(),
        ddl,w32,
        dig.Row(),dig.Dilogic(),dig.Addr(),
        dig.LorsX(),dig.LorsY()                            ),4);
  }    
  if(gPad->GetEvent()==1){
    new TCanvas("zoom",Form("Row %i DILOGIC %i",dig.Row(),dig.Dilogic()));  
  }
}//DrawZoom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawPc(Bool_t isFill) 
{ 
// Utility methode draws HMPID chamber PCs on event display.
// Arguments: none
//   Returns: none      
//  y6  ----------  ----------
//      |        |  |        |
//      |    4   |  |    5   |
//  y5  ----------  ----------
//
//  y4  ----------  ----------
//      |        |  |        |
//      |    2   |  |    3   |   view from electronics side
//  y3  ----------  ----------
//          
//  y2  ----------  ----------
//      |        |  |        |
//      |    0   |  |    1   |
//  y1  ----------  ----------
//      x1      x2  x3       x4
  gPad->Range(-10,-10,AliHMPIDDigit::SizeAllX()+5,AliHMPIDDigit::SizeAllY()+5); 
  Float_t x1=0,
          x2=AliHMPIDDigit::SizePcX(),
          x3=AliHMPIDDigit::SizePcX()+AliHMPIDDigit::SizeDead(),   
          x4=AliHMPIDDigit::SizeAllX();
  Float_t y1=0,
          y2=  AliHMPIDDigit::SizePcY(),
          y3=  AliHMPIDDigit::SizePcY()+AliHMPIDDigit::SizeDead(),
          y4=2*AliHMPIDDigit::SizePcY()+AliHMPIDDigit::SizeDead(),
          y5=  AliHMPIDDigit::SizeAllY()-AliHMPIDDigit::SizePcY(),
          y6=  AliHMPIDDigit::SizeAllY();

  Float_t xL[5]={x1,x1,x2,x2,x1}; //clockwise
  Float_t xR[5]={x3,x3,x4,x4,x3};  
  Float_t yD[5]={y1,y2,y2,y1,y1};
  Float_t yC[5]={y3,y4,y4,y3,y3};  
  Float_t yU[5]={y5,y6,y6,y5,y5};
    
  Float_t dX2=0.5*AliHMPIDDigit::SizePadX(),
          dY2=0.5*AliHMPIDDigit::SizePadY() ;
  
  TLatex txt; txt.SetTextSize(0.01);
  Int_t iColLeft=29,iColRight=41;
  TPolyLine *pc=0;  TLine *pL;
  AliHMPIDDigit dig;
  for(Int_t iPc=0;iPc<AliHMPIDDigit::kPcAll;iPc++){
    if(iPc==4) pc=new TPolyLine(5,xL,yU); if(iPc==5) pc=new TPolyLine(5,xR,yU); //draw PCs
    if(iPc==2) pc=new TPolyLine(5,xL,yC); if(iPc==3) pc=new TPolyLine(5,xR,yC);
    if(iPc==0) pc=new TPolyLine(5,xL,yD); if(iPc==1) pc=new TPolyLine(5,xR,yD);
    (iPc%2)? pc->SetFillColor(iColLeft): pc->SetFillColor(iColRight);
    if(!isFill){ pc->Draw();continue;}
    
    pc->Draw("f");
    if(iPc%2) {dig.Manual2(0,iPc,79,25); txt.DrawText(dig.LorsX()+2,dig.LorsY(),Form("PC%i",dig.Pc()));}//print PC#    
    
    txt.SetTextAlign(32);
    for(Int_t iRow=0;iRow<8 ;iRow++){//draw row lines (horizontal)
      dig.Manual2(0,iPc,0,iRow*6);   //set digit to the left-down pad of this row
      if(iPc%2) txt.DrawText(dig.LorsX()-1           ,dig.LorsY(),Form("%i",dig.PadPcY()));                  //print PadY#    
                txt.DrawText(dig.LorsX()-1+(iPc%2)*67,dig.LorsY()+2,Form("r%i",dig.Row()));                  //print Row#    
      pL=new TLine(dig.LorsX()-dX2,dig.LorsY()-dY2,dig.LorsX()+AliHMPIDDigit::SizePcX()-dX2,dig.LorsY()-dY2);//draw horizontal line 
      if(iRow!=0) pL->Draw(); 
    }//row loop  
    
    txt.SetTextAlign(13);
    for(Int_t iDil=0;iDil<10;iDil++){//draw dilogic lines (vertical)
      dig.Manual2(0,iPc,iDil*8,0);       //set this digit to the left-down pad of this dilogic        
                           txt.DrawText(dig.LorsX()  ,dig.LorsY()-1,Form("%i",dig.PadPcX()));                 //print PadX# 
      if(iPc==4 || iPc==5) txt.DrawText(dig.LorsX()+2,dig.LorsY()+42,Form("d%i",dig.Dilogic()));              //print Dilogic#    
      pL=new TLine(dig.LorsX()-dX2,dig.LorsY()-dY2,dig.LorsX()-dX2,dig.LorsY()+AliHMPIDDigit::SizePcY()-dY2); //draw vertical line
      if(iDil!=0)pL->Draw();
    }//dilogic loop        
  }//PC loop      
}//DrawPc()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void hed()
{//event display from files
  static TCanvas *pC=0;
  static Int_t iEvt=0;
  static Int_t iEvtTot=0;
  static TFile *pEsdFl=0;
  static TTree *pEsdTr=0;
  static AliESD *pEsd=0;
  
  if(!pC){
    if(hl==0) {Printf("hed: no HMPID loader");return;}
    Printf("Opening session");
    pEsdFl=TFile::Open("AliESDs.root");     if(!pEsdFl || !pEsdFl->IsOpen()) return;//open AliESDs.root
    pEsdTr=(TTree*) pEsdFl->Get("esdTree"); if(!pEsdTr)                      return;//get ESD tree
    pEsdTr->SetBranchAddress("ESD", &pEsd);
    hl->LoadDigits(); hl->LoadRecPoints();
    iEvtTot=pEsdTr->GetEntries();
    pC=new TCanvas("hed","View from electronics side, IP is behind the picture.",1000,900);  pC->ToggleEventStatus(); pC->Divide(3,3);
    pC->cd(7); TButton *pBtn=new TButton("Next","hed()",0,0,0.2,0.1);   pBtn->Draw(); 
  }
 
  if(iEvt<iEvtTot){
    pEsdTr->GetEntry(iEvt); al->GetEvent(iEvt); hl->TreeD()->GetEntry(0); hl->TreeR()->GetEntry(0);
    TLatex txt;   pC->cd(3);  txt.DrawLatex(0.2,0.2,Form("Event %i Total %i",iEvt,iEvtTot));
    DrawEvt(pC,h->DigLst(),h->CluLst(),pEsd);
    iEvt++;
  }else{
    Printf("Last event");
    pC->Clear();
    delete pC;pC=0x0;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void sed()
{

  static TCanvas *pC1=0;
  
  if(!pC1){
    pC1=new TCanvas("hed","Simulated evets-View from electronics side, IP is behind the picture.",1000,900); pC1->Divide(3,3);
    pC1->cd(7); TButton *pBtn=new TButton("Next","sed()",0,0,0.2,0.1);   pBtn->Draw(); 
  }


  
  AliHMPIDRecon rec;  

  TClonesArray lh("AliHMPIDHit"); 
  TClonesArray ls("AliHMPIDDigit"); 
  TObjArray    ld(7); for(Int_t i=0;i<7;i++) ld.AddAt(new TClonesArray("AliHMPIDDigit"),i);
  TObjArray    lc(7); for(Int_t i=0;i<7;i++) lc.AddAt(new TClonesArray("AliHMPIDCluster"),i);
  AliESD esd;
  
  
  EsdFromStack(&esd);
  HitsFromEsd(&esd,&lh);
  



             AliHMPIDv1::Hit2Sdi(&lh,&ls);                               
      AliHMPIDDigitizer::Sdi2Dig(&ls,&ld);     
  AliHMPIDReconstructor::Dig2Clu(&ld,&lc);
//        AliHMPIDTracker::Recon(&esd,&cl);
  
  DrawEvt(pC1,&ld,&lc,&esd);  
}//SimEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawEvt(TCanvas *pC,TObjArray *pDigLst,TObjArray *pCluLst,AliESD *pEsd)
{//draws all the objects of current event

  AliHMPIDRecon rec;  
  TPolyMarker *pTxC[7];  TPolyMarker *pRin[7]; //intesections and rings
  for(Int_t ch=0;ch<7;ch++){
    pTxC[ch]=new TPolyMarker; pTxC[ch]->SetMarkerStyle(2); pTxC[ch]->SetMarkerColor(kRed); pTxC[ch]->SetMarkerSize(3);
    pRin[ch]=new TPolyMarker; pRin[ch]->SetMarkerStyle(6); pRin[ch]->SetMarkerColor(kMagenta);
  }
  
  
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop to collect cerenkov rings and intersection points
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);
    Int_t ch=pTrk->GetHMPIDcluIdx();
    if(ch<0) continue; //this track does not hit HMPID
    ch/=1000000; 
    Float_t th,ph,xPc,yPc; pTrk->GetHMPIDtrk(xPc,yPc,th,ph);  //get info on current track
    pTxC[ch]->SetNextPoint(xPc,yPc);                          //add new intersection point
    
    Float_t ckov=pTrk->GetHMPIDsignal();  Float_t err=TMath::Sqrt(pTrk->GetHMPIDchi2());
    
    if(ckov>0){
      rec.SetTrack(xPc,yPc,th,ph);
     TVector2 pos;  for(int j=0;j<100;j++){rec.TracePhot(ckov,j*0.0628,pos); pRin[ch]->SetNextPoint(pos.X(),pos.Y());}      
    }
  }//tracks loop
      
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop    
    switch(iCh){
      case 6: pC->cd(1); break; case 5: pC->cd(2); break;
      case 4: pC->cd(4); break; case 3: pC->cd(5); break; case 2: pC->cd(6); break;
                                case 1: pC->cd(8); break; case 0: pC->cd(9); break;
    }
    gPad->SetEditable(kTRUE); gPad->Clear();
    DrawPc(0);
    ((TClonesArray*)pDigLst->At(iCh))->Draw();  //draw digits
    ((TClonesArray*)pCluLst->At(iCh))->Draw();  //draw clusters
                            pTxC[iCh]->Draw();  //draw intersections
                            pRin[iCh]->Draw();  //draw rings
    gPad->SetEditable(kFALSE);
  }//chambers loop
//  TLatex txt; txt.SetTextSize(0.02);
//  txt.DrawLatex(20,-5,Form("#theta=%.4f#pm%.5f with %2i #check{C}"          ,simCkov,simErr,simN));
//  txt.DrawLatex(25,-5,Form("#theta=%.4f#pm%.5f with %2i #check{C}"          ,recCkov,recErr,recN));
//  txt.DrawLatex(0 ,127,Form("#theta=%.2f#circ   #phi=%.2f#circ @(%.2f,%.2f) ",th*TMath::RadToDeg(),ph*TMath::RadToDeg(),radx,rady));
//  Printf("DIG------DIG---------DIG--------DIG------DIG------DIG");pDigLst->Print();Printf("");                   
}//DrawEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void t1(Int_t case=1)
{
  AliHMPIDDigit *d[10]; for(Int_t i=0;i<10;i++) d[i]=new AliHMPIDDigit;
  
  
  Int_t iNdig;
  
  if(case==1){
    iNdig=9;  
  
                                                              d[0]->Manual2(1,2,67,26, 33); 
                                d[1]->Manual2(1,2,66,25,431); d[2]->Manual2(1,2,67,25, 21);
  d[3]->Manual2(1,2,65,24,127); d[4]->Manual2(1,2,66,24, 54); d[5]->Manual2(1,2,67,24,  5);
  d[6]->Manual2(1,2,65,23, 20); d[7]->Manual2(1,2,66,23,  5); d[8]->Manual2(1,2,67,23,  6);
  }else if(case==2){
    iNdig=3;
    d[0]->Manual2(0,0,36,14,  8); 
    d[1]->Manual2(0,0,36,13, 33); d[2]->Manual2(0,0,37,13, 22);
  }
  
  AliHMPIDCluster c;
  for(int i=0;i<iNdig;i++) c.DigAdd(d[i]);  c.Print();
  
  
  TClonesArray *cl=new TClonesArray("AliHMPIDCluster");
  
  c.Solve(cl,kTRUE);
  Printf("");
  
  cl->Print();  
}//t1()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EsdFromStack(AliESD *pEsd)
{
  al->LoadHeader();al->LoadKinematics();
  AliStack *pStk=al->Stack();  
  
  for(Int_t iTrk=0;iTrk<pStk->GetNtrack();iTrk++){//stack loop
    TParticle *pPart=pStk->Particle(iTrk);
    if(pPart->GetPDG()->Charge()==0) continue; //neutral particles are not reconstructed
    if(pPart->GetFirstMother()>0)    continue; //do not consider secondaries
    AliESDtrack trk(pPart);
    pEsd->AddTrack(&trk);
  }//stack loop  
}//EsdFromStack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HitsFromEsd(AliESD *pEsd, TClonesArray *pHitLst)
{
  AliHMPIDRecon rec;
  const Int_t kCerenkov=50000050,kFeedback=50000051;
  Int_t hc=0; TVector2 pos;
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);
    Float_t xRa,yRa;
    Int_t ch=AliHMPIDTracker::IntTrkCha(pTrk,xRa,yRa);
    if(ch<0) continue; //this track does not hit HMPID
    Float_t ckov=0.63;

    Float_t th,ph,xPc,yPc,; pTrk->GetHMPIDtrk(xPc,yPc,th,ph);  
                          new((*pHitLst)[hc++]) AliHMPIDHit(ch,200e-9,kProton  ,iTrk,xPc                ,yPc);                 //mip hit
    for(int i=0;i<4;i++)  new((*pHitLst)[hc++]) AliHMPIDHit(ch,7.5e-9,kFeedback,iTrk,gRandom->Rndm()*130,gRandom->Rndm()*126); //bkg hits 4 per track
    for(int i=0;i<16;i++){
      rec.TracePhot(ckov,gRandom->Rndm()*TMath::TwoPi(),pos);
                          new((*pHitLst)[hc++]) AliHMPIDHit(ch,7.5e-9,kCerenkov,iTrk,pos.X(),pos.Y());}                      //photon hits  
  }//tracks loop    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
