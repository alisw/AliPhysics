AliRun     *a; AliRunLoader *al;   TGeoManager *g; //globals for easy manual manipulations
AliHMPID   *h; AliLoader    *hl; AliHMPIDParam *hp;
Bool_t isGeomType=kFALSE;

Int_t nCurEvt=0;
Int_t nMaxEvt=0;
TControlBar *pMenu=0;
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
    status+="galice.root: ";
    al=AliRunLoader::Open();                                                //try to open galice.root from current dir 
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    al->LoadgAlice(); a=al->GetAliRun();                                    //take new AliRun object from galice.root   
    hl=al->GetDetectorLoader("HMPID");  h=(AliHMPID*)a->GetDetector("HMPID");  //get HMPID object from galice.root
    
    status+=(h)? "HMPID": "PROBLEM PROBLEM PROBLEM- no HMPID";
    nMaxEvt=al->GetNumberOfEvents()-1;
    status+=Form(" Event(s) 0-%i",nMaxEvt); 
  }else  
    status+="PROBLEM PROBLEM PROBLEM no galice.root";
  
  status+=Form(" curent event %i",nCurEvt);
  GetParam();
  pMenu = new TControlBar("horizontal",status.Data(),0,0);
    pMenu->AddButton("                     ","","");
    pMenu->AddButton("       General       ","General()"  ,"general items which do not depend on any files");
    pMenu->AddButton("                     ",""           ,"");
    pMenu->AddButton("       Sim data      ","SimData()"  ,"items which expect to have simulated files"    );
    pMenu->AddButton("                     ",""           ,"");
    pMenu->AddButton("       Raw data      ","RawData()"  ,"items which expect to have raw files"          );
    pMenu->AddButton("                     ","       "    ,"");
    pMenu->AddButton("         Test        ","Test()"     ,"all test utilities");
    pMenu->AddButton("      PREV EVENT     ","PrevEvent()" ,"Set the previous event"             );
    pMenu->AddButton("      NEXT EVENT     ","NextEvent()","Set the next event"                  );
    pMenu->AddButton("         Quit        ",".q"         ,"close session"                       );
  pMenu->Show();
}//Menu()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void General()
{         
  TControlBar *pMenu = new TControlBar("vertical","General purpose",100,50);  
    pMenu->AddButton("Debug ON","don();"                   ,"Switch debug on-off"                        );   
    pMenu->AddButton("Debug OFF","doff();"                 ,"Switch debug on-off"                        );   
    pMenu->AddButton("Geo GUI","geo();"                    ,"Shows geometry"                             ); 
    pMenu->AddButton("Browser","new TBrowser;"             ,"Start ROOT TBrowser"                        );
  pMenu->Show();  
}//General()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimData()
{
  TControlBar *pSim = new TControlBar("vertical","Sim data",310,50);  
    pSim->AddButton("Display ","hed();"    ,"Display Fast");
    pSim->AddButton("HITS QA"           ,"hqa()"     ,"QA plots for hits: hqa()");
    pSim->AddButton("Print stack"       ,"stack();"  ,"To print hits:     hp(evt)");
    pSim->AddButton("Print hits"        ,"hp(nCurEvt);"     ,"To print hits:     hp(evt)");
    pSim->AddButton("Print sdigits"     ,"sp(nCurEvt);"     ,"To print sdigits:  sp(evt)");
    pSim->AddButton("Print digits"      ,"dp(nCurEvt);"     ,"To print digits:   dp(evt)");
    pSim->AddButton("Print clusters"    ,"cp(nCurEvt);"     ,"To print clusters: cp(evt)");
  pSim->Show();         
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
  TControlBar *pTst = new TControlBar("vertical","Test",625,50);  
    pTst->AddButton("TEST Display "      ,"sed();"                    ,"Display Fast");
    pTst->AddButton("Test all"           ,"tst();"                   ,"test hits->sdigits->digits"                 );   
    pTst->AddButton("Segmentation"       ,"ts()"                      ,"test segmentation methods"                  );
    pTst->AddButton("Test response"      ,"AliHMPIDParam::TestResp();","Test AliHMPIDParam response methods"         );
    pTst->AddButton("Print map"          ,"PrintMap();"               ,"Test AliHMPIDParam transformation methods"   );
    pTst->AddButton("Test Recon"         ,"rec();"                    ,"Test AliHMPIDRecon"                          );
  pTst->Show();  
}//Test()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void doff(){  Printf("DebugOFF");  AliLog::SetGlobalDebugLevel(0);}
void don() {  Printf("DebugON");   AliLog::SetGlobalDebugLevel(AliLog::kDebug);}

void geo (                       ) {gGeoManager->GetTopVolume()->Draw("ogl");}
  
void du  (                       ) {h->Dump         (   );}                //utility display 

void PrevEvent()                   {nCurEvt--;if(nCurEvt<0       )nCurEvt=0      ;pMenu->SetTitle(Form("Event(s): 0-%i Current event %i",nMaxEvt,nCurEvt));}
void NextEvent()                   {nCurEvt++;if(nCurEvt>=nMaxEvt)nCurEvt=nMaxEvt;pMenu->SetTitle(Form("Event(s): 0-%i Current event %i",nMaxEvt,nCurEvt));}
void stack(                     )  {AliHMPIDParam::Stack();}    
void tid  (Int_t tid,Int_t evt=0)  {AliHMPIDParam::Stack(evt,tid);} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintMap()
{
 
  Double_t r2d=TMath::RadToDeg();

  Double_t x=AliHMPIDDigit::SizeAllX(),y=AliHMPIDDigit::SizeAllY();
    
  Printf("\n\n\n");                                       
  
  for(int ch=6;ch>=0;ch--){
    AliHMPIDDigit dL,dR; dL.Manual2(ch,2,0 ,24);
                         dR.Manual2(ch,3,79,24);
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
void tst(Int_t nEvts=111,Int_t type=999)
{
  
  TLegend *lQ=new TLegend(0.5,0.5,0.9,0.9);

  TH1F *hQ7  =new TH1F("hQ7"  ,"" ,300,-50,2000);     hQ7  ->SetLineColor(kRed);      lQ->AddEntry(hQ7  ,"Ckov 7 eV");  hQ7->SetStats(0);
  TH1F *hQ200=new TH1F("hQ200","" ,300,-50,2000);     hQ200->SetLineColor(kBlack);    lQ->AddEntry(hQ200,"mip 200 eV");
  TH1F *hQ500=new TH1F("hQ500","" ,300,-50,2000);     hQ500->SetLineColor(kCyan);     lQ->AddEntry(hQ500,"mip 500 eV");
  TH1F *hQ900=new TH1F("hQ900","" ,300,-50,2000);     hQ900->SetLineColor(kGreen);    lQ->AddEntry(hQ900,"mip 900 eV");
  
  TH1F *hCluPerEvt=new TH1F("hCluPerEvt","# clusters per event",11,-0.5,10.5);
  TH1F *hCluChi2  =new TH1F("hChi2","Chi2 ",1000,0,100);
  TH1F *hCluFlg   =new TH1F("hCluFlg","Cluster flag",14,-1.5,12.5);                       hCluFlg->SetFillColor(5);
  TH1F *hCluRawSize= new TH1F("hCluRawSize","Raw cluster size ",100,0,100);
  
  TH2F *pCluMapSi1=new TH2F("cluMapSi1","Size 1 map"       ,1700,-10,160,1700,-10,160);           
  TH2F *pCluMapLo0=new TH2F("cluMNoLo0","Loc Max 0 map"    ,1700,-10,160,1700,-10,160);     
  TH2F *pCluMapLo1=new TH2F("cluMapLo1","Loc Max 1 map"    ,1700,-10,160,1700,-10,160);        
  TH2F *pCluMapUnf=new TH2F("cluMapUnf","Unfolded map"     ,1700,-10,160,1700,-10,160);         
  TH2F *pCluMapEdg=new TH2F("cluMapEdg","On edge map"      ,1700,-10,160,1700,-10,160);         
  TH2F *pCluMapCoG=new TH2F("cluMapCoG","CoG  map"         ,1700,-10,160,1700,-10,160);            
  TH2F *pCluMapEmp=new TH2F("cluMapEmp","undefined-empty"  ,1700,-10,160,1700,-10,160);      

  TH1F *hHitCluDifX = new TH1F("hHitCluDifX" ,";entries;x_{Hit}-x_{Clu} [cm]"   ,2000,-2,2);          hHitCluDifX->Sumw2();    hHitCluDifX->SetFillColor(kYellow);
  TH1F *hHitCluDifY = new TH1F("hHitCluDifY" ,";entries;y_{Hit}-y_{Clu} [cm]"   ,2000,-2,2);          hHitCluDifY->Sumw2();    hHitCluDifY->SetFillColor(kYellow);
  TH2F *hHitCluDifXY= new TH2F("hHitCluDifXY",";x_{Hit}-x_{Clu};y_{Hit}-y_{Clu}",2000,-2,2,2000,-2,2);hHitCluDifXY->Sumw2();  
  TH1F *hHitCluDifQ = new TH1F("hHitCluDifQ" ,";entries;Q_{Hit}-Q_{Clu}"        ,200 ,-100,100);      hHitCluDifQ->Sumw2();    hHitCluDifQ->SetFillColor(kYellow);
 
  Float_t e200=200e-9,e500=500e-9,e900=900e-9,e7=7e-9;//predefined  Eloss
  
  
  AliHMPIDHit hit(0,0,kProton,0,0,0);
  for(int i=0;i<10000;i++){hQ200->Fill(hit.QdcTot(e200));  hQ500->Fill(hit.QdcTot(e500));   hQ900->Fill(hit.QdcTot(e900));  hQ7  ->Fill(hit.QdcTot(e7));}
  
  TClonesArray hits("AliHMPIDHit");  TClonesArray sdigs("AliHMPIDDigit");
  TObjArray digs(7); for(Int_t i=0;i<7;i++) digs.AddAt(new TClonesArray("AliHMPIDDigit"),i);
  TObjArray clus(7); for(Int_t i=0;i<7;i++) clus.AddAt(new TClonesArray("AliHMPIDCluster"),i);
  
    
  for(Int_t iEvt=0;iEvt<nEvts;iEvt++){//events loop
    if(iEvt%500==0)Printf("============> iEvt = %d ",iEvt);
    
    Int_t ch,pid; Float_t e,hitx,hity,hitq;
    Int_t nHits=(type==999)?1:40;
    for(Int_t iHit=0;iHit<nHits;iHit++){//hits loop for the current event
      switch(iHit){
        case 0:  ch=0;pid=kProton;e=e200;hitx=16.0+gRandom->Rndm()*0.8;hity= 16.8+gRandom->Rndm()*0.84;break; //mip ramdomly distributed in one pad in the middle
        case 1:  ch=0;pid=kProton;e=e200;hitx=0.4;hity=0.42;break; //mip in left-hand bottom coner of chamber 0
        case 2:  ch=0;pid=kProton;e=e200;hitx=0.4;hity=30  ;break; //mip on left edge of chamber 0
        case 3:  ch=0;pid=kProton;e=e200;hitx=40; hity=0.42;break; //mip on bottom edge of chamber 0
        default: ch=gRandom->Rndm()*6; pid=(gRandom->Rndm()>0.9)? kProton:kCerenkov;
                  if(pid==kProton) 
                    e=gRandom->Rndm()*900e-9; 
                  else
                    e=5.5e-9+3e-9*gRandom->Rndm();
                  x=gRandom->Rndm()*AliHMPIDDigit::SizeAllX(); y=gRandom->Rndm()*AliHMPIDDigit::SizeAllY();break; //random hit
      }
      new(hits[iHit]) AliHMPIDHit(ch,e,pid,iHit,hitx,hity);                          
      hitq=e;
    }//hits loop
    
    AliHMPIDv1::Hit2Sdi(&hits,&sdigs);
    AliHMPIDDigitizer::Sdi2Dig(&sdigs,&digs);     
    AliHMPIDReconstructor::Dig2Clu(&digs,&clus);
          
    for(Int_t iCh=AliHMPIDDigit::kMinCh;iCh<=AliHMPIDDigit::kMaxCh;iCh++){//chambers loop
      TClonesArray *pDigs=(TClonesArray *)digs.UncheckedAt(iCh);
      TClonesArray *pClus=(TClonesArray *)clus.UncheckedAt(iCh);
        
      hCluPerEvt->Fill(pClus->GetEntriesFast());
      for(Int_t iClu=0;iClu<pClus->GetEntriesFast();iClu++){//clusters loop
        AliHMPIDCluster *pClu=(AliHMPIDCluster*)pClus->UncheckedAt(iClu);
        Float_t clux=pClu->X(); Float_t cluy=pClu->Y(); Float_t cluq=pClu->Q();
        hCluFlg->Fill(pClu->Status());
        hCluChi2->Fill(pClu->Chi2());
        hCluRawSize->Fill(pClu->Size());
	 
        switch(pClu->Status()){
            case AliHMPIDCluster::kSi1:   pCluMapSi1->Fill(clux,cluy); break;
            case AliHMPIDCluster::kLo1:   pCluMapLo1->Fill(clux,cluy); break;
            case AliHMPIDCluster::kUnf:   pCluMapUnf->Fill(clux,cluy); break; 
            case AliHMPIDCluster::kMax:   pCluMapMax->Fill(clux,cluy); break;
            case AliHMPIDCluster::kEdg:   pCluMapEdg->Fill(clux,cluy); break;
            case AliHMPIDCluster::kCoG:   pCluMapCoG->Fill(clux,cluy); break;
            case AliHMPIDCluster::kNoLoc: pCluMapNoLoc->Fill(clux,cluy); break;
            default:     pCluMapEmp->Fill(clux,cluy); break; 		
        }
        
        hHitCluDifX->Fill(hitx-clux); hHitCluDifY->Fill(hity-cluy); hHitCluDifXY->Fill(hitx-clux,hity-cluy); hHitCluDifQ->Fill(hitq-cluq);
        
      }//clusters loop
    }//chambers loop
      
    hits.Delete();  sdigs.Delete();  for(int i = 0;i<7;i++){((TClonesArray*)digs.At(i))->Delete();((TClonesArray*)clus.At(i))->Delete();}
  }//events loop      
      
  gStyle->SetPalette(1);
  TCanvas *pC2=new TCanvas("Digit canvas","Digit canvas",1280,800); pC2->Divide(3,3);
  pC2->cd(1);gPad->SetLogy(1);hHitCluDifX->Draw("hist");
  pC2->cd(2);gPad->SetLogy(1);hHitCluDifY->Draw("hist");
  pC2->cd(3);gPad->SetLogz(1);hHitCluDifXY->Draw("colz");
  pC2->cd(4);gPad->SetLogy(1);hHitCluDifQ->Draw("hist");
  pC2->cd(5);gPad->SetLogy(1);hCluFlg->Draw();
  pC2->cd(6);gPad->SetLogy(1);hCluChi2->Draw();
  pC2->cd(7);                 hCluRawSize->Draw();
  pC2->cd(8);                 hCluPerEvt->Draw("colz");
  pC2->cd(9);                 hQckov->Draw(); hQm200->Draw("same"); hQm500->Draw("same"); hQm900->Draw("same"); lQ->Draw();
  TCanvas *pC1=new TCanvas("ClusterMaps","Cluster maps",1280,800); pC1->Divide(3,3);
  pC1->cd(1);  pCluMapSi1->Draw();  DrawPc(kFALSE);
  pC1->cd(2);  pCluMapLo1->Draw();  DrawPc(kFALSE);
  pC1->cd(3);  pCluMapUnf->Draw();  DrawPc(kFALSE);
  pC1->cd(5);  pCluMapMax->Draw();  DrawPc(kFALSE);
  pC1->cd(6);  pCluMapEdg->Draw(); DrawPc(kFALSE);
  pC1->cd(7);  pCluMapCoG->Draw(); DrawPc(kFALSE);
  pC1->cd(8);  pCluMapNoLoc->Draw();DrawPc(kFALSE);
  pC1->cd(9);  pCluMapEmp->Draw(); DrawPc(kFALSE);
   
  pC1->SaveAs("$HOME/HitMaps.png");  //?????
  pC2->SaveAs("$HOME/HitCluDif.gif");  
  
  Printf("Digits - raw -digits conversion...");  
  
  AliHMPIDDigit d1,d2; Int_t ddl,r,d,a;UInt_t w32;
  for(Int_t iCh=AliHMPIDDigit::kMinCh;iCh<=AliHMPIDDigit::kMaxCh;iCh++)
  for(Int_t iPc=AliHMPIDDigit::kMinPc;iPc<=AliHMPIDDigit::kMaxPc;iPc++)
  for(Int_t iPx=AliHMPIDDigit::kMinPx;iPx<=AliHMPIDDigit::kMaxPx;iPx++)
  for(Int_t iPy=AliHMPIDDigit::kMinPy;iPy<=AliHMPIDDigit::kMaxPy;iPy++){
    d1.Set(iCh,iPc,iPx,iPy,3040);   //set digit               
    d1.Raw(w32,ddl,r,d,a);          //get raw word for this digit 
    d2.Raw(w32,ddl);                //set another digit from that raw word
    if(d1.Compare(&d2)) {d1.Print(); d2.Print(); Printf("");}//compare them
  }
  Printf("OK");
}//tst()


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void hp(Int_t iEvt=0)
{
//Prints a list of HMPID hits for a given event. Default is event number 0.
  Printf("List of HMPID hits for event %i",iEvt);
  if(al->GetEvent(iEvt)) return;    
  if(hl->LoadHits()) return;
  
  Int_t iTotHits=0;
  for(Int_t iPrim=0;iPrim<hl->TreeH()->GetEntries();iPrim++){//prims loop
    hl->TreeH()->GetEntry(iPrim);      
    h->Hits()->Print();
    iTotHits+=h->Hits()->GetEntries();
  }
  hl->UnloadHits();
  Printf("totally %i hits for event %i",iTotHits,iEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void sp(Int_t iEvt=0)
{
//prints a list of HMPID sdigits  for a given event
  Printf("List of HMPID sdigits for event %i",iEvt);
  if(al->GetEvent(iEvt)) return;    
  if(hl->LoadSDigits()) return;
  
  hl->TreeS()->GetEntry(0);
  h->SdiLst()->Print();
  hl->UnloadSDigits();
  Printf("totally %i sdigits for event %i",h->SdiLst()->GetEntries(),iEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dp(Int_t iEvt=0)
{
//prints a list of HMPID digits  for a given event
  Printf("List of HMPID digits for event %i",iEvt);
  if(al->GetEvent(iEvt)) return;    
  if(hl->LoadDigits()) return;
  
  hl->TreeD()->GetEntry(0);
  h->DigLst()->Print();
  Int_t totDigs=0;
  for(Int_t i=0;i<7;i++) {totDigs+=h->DigLst(i)->GetEntries();}
  hl->UnloadDigits();
  Printf("totally %i digits for event %i",totDigs,iEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cp(Int_t iEvt=0)
{//prints a list of HMPID clusters  for a given event
  Printf("List of HMPID clusters for event %i",iEvt);
  if(al->GetEvent(iEvt)) return;    
  if(hl->LoadRecPoints()) return;
  
  hl->TreeR()->GetEntry(0);
  h->CluLst()->Print();
  
  Int_t iCluCnt=0; for(Int_t iCh=0;iCh<7;iCh++) iCluCnt+=h->CluLst(iCh)->GetEntries();
  
  hl->UnloadRecPoints();
  Printf("totally %i clusters for event %i",iCluCnt,iEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ttt()
{
  TClonesArray hits ("AliHMPIDDigit");
  TClonesArray sdigs("AliHMPIDDigit");
  
  AliHMPIDHit hit(0,45e-9,kProton,33,0,0);
  hit.Hit2Sdi(&sdigs);
  sdigs.Print();
}


#include "Hdisp.C"
