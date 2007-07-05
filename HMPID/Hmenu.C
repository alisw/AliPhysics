AliRun     *a; AliRunLoader *al;   TGeoManager *g; //globals for easy manual manipulations
AliHMPID   *h; AliLoader    *hl; AliHMPIDParam *hp;

Int_t nCurEvt=0;
Int_t nMaxEvt=0;
TControlBar *pMenu=0;
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

  AliHMPIDParam::Instance();      // geometry loaded

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

  Double_t x=AliHMPIDParam::SizeAllX(),y=AliHMPIDParam::SizeAllY();
    
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
