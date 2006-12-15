AliRun *a;     AliRunLoader *al;   TGeoManager *g; //globals for easy manual manipulations
AliHMPID   *r; AliLoader    *rl; AliHMPIDParam *rp;
Bool_t isGeomType=kFALSE;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GetParam()
{
  isGeomType=!isGeomType;
  if(g) delete g;  if(rp) delete rp; //delete current TGeoManager and AliHMPIDParam
  if(isGeomType) g=TGeoManager::Import("geometry.root");
  else           g=TGeoManager::Import("misaligned_geometry.root");
  rp=AliHMPIDParam::Instance();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hmenu()
{   
  TString status="Status: ";
  if(gSystem->IsFileInIncludePath("galice.root")){
    status+="galice.root found";
    al=AliRunLoader::Open();                                                //try to open galice.root from current dir 
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    al->LoadgAlice(); a=al->GetAliRun();                                    //take new AliRun object from galice.root   
    rl=al->GetDetectorLoader("HMPID");  r=(AliHMPID*)a->GetDetector("HMPID");  //get HMPID object from galice.root
    
    status+=Form(" with %i event(s)",al->GetNumberOfEvents()); status+=(r)? " with HMPID": " without HMPID";
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
    pMenu->AddButton("Test response"      ,"AliHMPIDParam::TestResp();","Test AliHMPIDParam response methods"         );
    pMenu->AddButton("Print map"          ,"PrintMap();"               ,"Test AliHMPIDParam transformation methods"   );
    pMenu->AddButton("Test Recon"         ,"rec();"                    ,"Test AliHMPIDRecon"                          );
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
void dq  (                       ) {AliHMPIDReconstructor::DigQA     (al );}   //digits QA plots for all events
void cp  (Int_t evt=0            ) {r->CluPrint  (evt);                   }   //print clusters for requested event
void cq  (                       ) {AliHMPIDReconstructor::CluQA     (al );}   //clusters QA plots for all events

void ep  (                       ) {AliHMPIDTracker::EsdQA(1);            } 
void eq  (                       ) {AliHMPIDTracker::EsdQA();             }                   
void mp  (Double_t probCut=0.7   ) {AliHMPIDTracker::MatrixPrint(probCut);}                   


void t   (Int_t evt=0          )   {AliHMPIDParam::Stack(evt);}    
void tid (Int_t tid,Int_t evt=0)   {AliHMPIDParam::Stack(evt,tid);} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void tst()
{
//create all lists  
  TClonesArray hit("AliHMPIDHit"),sdi("AliHMPIDDigit");
  TObjArray    dig,clu,cog;      for(Int_t i=0;i<7;i++) {dig.AddAt(new TClonesArray("AliHMPIDDigit"),i); 
                                                         clu.AddAt(new TClonesArray("AliHMPIDCluster"),i);
                                                         cog.AddAt(new TClonesArray("AliHMPIDCluster"),i);}
//simulate track  
  TLorentzVector mom; mom.SetPtEtaPhiM(3,0,30*TMath::DegToRad(),0.938);                                                                    
  TLorentzVector vtx; vtx.SetXYZT(0,0,0,0);                                                                 

  AliESD *pESD=new AliESD;  pESD->SetMagneticField(0.2);
  pESD->AddTrack(new AliESDtrack(new TParticle(kProton,0,-1,-1,-1,-1,mom,vtx)));
  pESD->Print();
  
  AliTracker::SetFieldMap(new AliMagF,1);
  AliHMPIDTracker *pTracker=new AliHMPIDTracker;
  
        
  return;                                                                      
  Double_t th=8*TMath::DegToRad();                              //gRandom->Rndm()*TMath::PiOver4();
  Double_t ph=gRandom->Rndm()*TMath::TwoPi(); 
  Double_t radx=gRandom->Rndm()*AliHMPIDDigit::SizeAllX(); 
  Double_t rady=gRandom->Rndm()*AliHMPIDDigit::SizeAllY();
  
  Int_t iHitCnt=0;  
  
  
  
  
      
                                                                      
             AliHMPIDv1::Hit2Sdi(&hitLst,&sdiLst);                               
      AliHMPIDDigitizer::Sdi2Dig(&sdiLst,&digLst);     
  AliHMPIDReconstructor::Dig2Clu(&digLst,&cluLst);   AliHMPIDReconstructor::Dig2Clu(&digLst,&cogLst,0);  
  
  
  Int_t iDigN=0,iCluN=0,iCogN=0; for(Int_t i=0;i<7;i++){ iDigN+=((TClonesArray*)digLst.At(i))->GetEntries();                  
                                                         iCluN+=((TClonesArray*)cluLst.At(i))->GetEntries();
                                                         iCogN+=((TClonesArray*)cogLst.At(i))->GetEntries(); }                   
  
  Printf("SDI------SDI---------SDI--------SDI------SDI------SDI #%i",sdiLst.GetEntries());sdiLst.Print();Printf("");
  Printf("DIG------DIG---------DIG--------DIG------DIG------DIG #%i",iDigN              );digLst.Print();Printf("");                   
  Printf("HIT------HIT---------HIT--------HIT------HIT------HIT #%i",hitLst.GetEntries());hitLst.Print();Printf("");
  Printf("CLU------CLU---------CLU--------CLU------CLU------CLU #%i",iCluN              );cluLst.Print();Printf("");                     
  Printf("COG------COG---------COG--------COG------COG------COG #%i",iCogN              );cogLst.Print();Printf("");                     
}



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


void rec()
{
  
  
  Double_t ckovMax=0.75,ckovSim;
  Int_t nSim=0;
  while(nSim<3){
    ckovSim=gRandom->Rndm()*ckovMax;//0.6468;
    nSim=20*TMath::Power(TMath::Sin(ckovSim)/TMath::Sin(ckovMax),2); //scale number of photons 
  }
  
  
  TClonesArray *pCluLst=new TClonesArray("AliHMPIDCluster");
  TPolyMarker  *pMipMap=new TPolyMarker();   pMipMap->SetMarkerStyle(8);  pMipMap->SetMarkerColor(kRed); 
  TPolyMarker  *pPhoMap=new TPolyMarker();   pPhoMap->SetMarkerStyle(4);  pPhoMap->SetMarkerColor(kRed);
  TPolyMarker  *pBkgMap=new TPolyMarker();   pBkgMap->SetMarkerStyle(25); pBkgMap->SetMarkerColor(kRed);
  TPolyLine    *pRing  =new TPolyLine;                                    pRing->SetLineColor(kGreen);     
  TPolyLine    *pOld   =new TPolyLine;                                    pOld->SetLineColor(kBlue);     
  
  Int_t iCluCnt=0; pMipMap->SetPoint(iCluCnt,x,y); new((*pCluLst)[iCluCnt++])   AliHMPIDCluster(1,x,y,200); //add mip cluster
  
//  for(int j=0;j<30;j++){                                                                                   //add bkg photons  
//    Float_t bkgX=gRandom->Rndm()*AliHMPIDDigit::SizeAllX();
//    Float_t bkgY=gRandom->Rndm()*AliHMPIDDigit::SizeAllY();
//    pBkgMap->SetPoint(iCluCnt,bkgX,bkgY); new((*pCluLst)[iCluCnt++]) AliHMPIDCluster(1,bkgX,bkgY,35);
//  }   

  
  
  
  AliHMPIDRecon    rec; rec.SetTrack(th,ph,x,y);                                                                
  
  TVector2 pos;
  for(int i=0;i<nSim;i++){
    rec.TracePhot(ckovSim,gRandom->Rndm()*2*TMath::Pi(),pos);                                   //add photons 
    if(AliHMPIDDigit::IsInDead(pos.X(),pos.Y())) continue; 
    pPhoMap->SetPoint(iCluCnt,pos.X(),pos.Y());    new((*pCluLst)[iCluCnt++]) AliHMPIDCluster(1,pos.X(),pos.Y(),35); 
  }  

  
  Int_t nRec=0,nOld=0;
  Double_t ckovRec=rec.CkovAngle(pCluLst,nRec); Double_t err=TMath::Sqrt(rec.CkovSigma2());   
  Double_t ckovOld=old.CkovAngle(pCluLst,nOld);
  
  Printf("---------------- Now reconstructed --------------------");
  
  
  for(int j=0;j<100;j++){rec.TracePhot(ckovRec,j*0.0628,pos); pRing->SetPoint(j,pos.X(),pos.Y());}  
  for(int j=0;j<100;j++){rec.TracePhot(ckovOld,j*0.0628,pos); pOld->SetPoint(j,pos.X(),pos.Y());}  
    
  new TCanvas;  AliHMPIDDigit::DrawPc();  pMipMap->Draw(); pPhoMap->Draw(); pBkgMap->Draw(); pRing->Draw();  pOld->Draw(); 
  
  TLatex txt; txt.SetTextSize(0.03);
  txt.DrawLatex(65,127,Form("#theta=%.4f#pm%.5f with %2i #check{C}"                             ,ckovSim, 0.,nSim             ));
  txt.DrawLatex(65,122,Form("#theta=%.4f#pm%.5f with %2i #check{C} Old=%.4f with %i #check{C}"  ,ckovRec,err,nRec,ckovOld,nOld));
  txt.DrawLatex(0 ,127,Form("#theta=%.2f#circ   #phi=%.2f#circ @(%.2f,%.2f) ",th*TMath::RadToDeg(),ph*TMath::RadToDeg(),x,y));
                   
//  for(int i=0;i<35;i++){
//    Double_t ckov=0.1+i*0.02;
//    Printf("Ckov=%.2f Old=%.3f New=%.3f",ckov,old.FindRingArea(ckov),rec.FindRingArea(ckov));
//  }
  
  
}//rec()


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
                       

  TH2F *pEleHitRZ    =new TH2F("EleHitRZ"    ,Form("e^{+} e^{-} hit %s;z[cm];R[cm]" ,GetName())     , 400,-300,300 ,400,-500,500);   //R-z plot 0cm<R<550cm -300cm<z<300cm  
  TH2F *pEleHitRP    =new TH2F("EleHitRP"    ,Form("e^{+} e^{-} hit %s;p[GeV];R[cm]",GetName())     ,1000,-1  ,1   ,400,   0,550);   //R-p plot 0cm<R<550cm -1GeV<p<1GeV 
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


void RecWithStack()
{
  al->LoadHeader();al->LoadKinematics();
  AliStack *pStk=al->Stack();  
  
  AliESD *pEsd=new AliESD;
  for(Int_t iTrk=0;iTrk<pStk->GetNtrack();iTrk++){//stack loop
    TParticle *pPart=pStk->Particle(iTrk);
    if(pPart->GetPDG()->Charge()==0) continue; //neutral particles are not reconstructed
    pEsd->AddTrack(new AliESDtrack(pPart));
  }//stack loop
  
  pEsd->Print();
  AliTracker::SetFieldMap(new AliMagF,1);
  AliHMPIDTracker t; 
  rl->LoadRecPoints(); 
  t.LoadClusters(rl->TreeR()); 
  t.PropagateBack(pEsd);
  rl->UnloadRecPoints();
}


void AliHMPIDReconstructor::CluQA(AliRunLoader *pAL)
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
  
        
  TH1F*    pQ=new TH1F("RiAllQ"  ,"Charge All"           ,4000 ,0  ,4000);// Q hists
  TH1F* pCerQ=new TH1F("RiCerQ"  ,"Charge Ckov"          ,4000 ,0  ,4000);
  TH1F* pMipQ=new TH1F("RiMipQ"  ,"Charge MIP"           ,4000 ,0  ,4000);
  
  TH1F*    pS=new TH1F("RichCluSize"    ,"Cluster size;size"         ,100  ,0  ,100 );// size hists
  TH1F* pCerS=new TH1F("RichCluCerSize" ,"Ckov size;size"            ,100  ,0  ,100 );
  TH1F* pMipS=new TH1F("RichCluMipSize" ,"MIP size;size"             ,100  ,0  ,100 );
  
  TH2F*    pM=new TH2F("RichCluMap"     ,"Cluster map;x [cm];y [cm]" ,1000 ,0  ,AliHMPIDDigit::SizePcX(),1000,0,AliHMPIDDigit::SizePcY()); // maps
  TH2F* pMipM=new TH2F("RichCluMipMap"  ,"MIP map;x [cm];y [cm]"     ,1000 ,0  ,AliHMPIDDigit::SizePcX(),1000,0,AliHMPIDDigit::SizePcY());
  TH2F* pCerM=new TH2F("RichCluCerMap"  ,"Ckov map;x [cm];y [cm]"    ,1000 ,0  ,AliHMPIDDigit::SizePcX(),1000,0,AliHMPIDDigit::SizePcY());
 
  
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){
    pAL->GetEvent(iEvt);               
    pRL->TreeD()->GetEntry(0); 
    TClonesArray *pCluLst=new TClonesArray("AliHMPIDCluster");//tmp list of clusters for this event
    
    for(Int_t iCh=0;iCh<7;iCh++) Dig2Clu(pRich->DigLst(iCh),pCluLst,kFALSE);//cluster finder for all chamber if any digits present
    
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



void AliHMPID::OccupancyPrint(Int_t iEvtNreq)
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
}


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
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
