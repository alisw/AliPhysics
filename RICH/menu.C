AliRICH * R()    {return r;}
void ph(Int_t event=0) {R()->PrintHits(event);}    //utility print hits for 'event' event
void ps(Int_t event=0) {R()->PrintSDigits(event);} //utility print sdigits
void pd(Int_t event=0) {R()->PrintDigits(event);}  //utility print digits
void pc(Int_t event=0) {R()->PrintClusters(event);}//utility print clusters
void pt(Int_t event=0) {R()->PrintTracks(event);}  //utility print tracks

//__________________________________________________________________________________________________
void pp(int tid)
{
  R()->GetLoader()->GetRunLoader()->LoadHeader();  R()->GetLoader()->GetRunLoader()->LoadKinematics();
  
  if(tid<0||tid>=R()->GetLoader()->GetRunLoader()->Stack()->GetNtrack())
    cout<<"Valid tid number is 0-"<<R()->GetLoader()->GetRunLoader()->Stack()->GetNtrack()-1<<" for this event.\n";
  else
    PrintParticleInfo(tid);
  
  R()->GetLoader()->GetRunLoader()->UnloadKinematics();  R()->GetLoader()->GetRunLoader()->UnloadHeader();
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

Double_t r2d = TMath::RadToDeg();
Double_t d2r = TMath::DegToRad();

void DisplFast(){ AliRICHDisplFast *d = new AliRICHDisplFast();  d->Exec();}  


void C_R()
{
  AliRICHRecon *detect = new AliRICHRecon("RICH patrec algorithm","Reconstruction");
    

  for (int nev=0; nev< a->GetEventsPerRun(); nev++) {    // Event Loop
    al->GetEvent(nev);
    cout <<endl<< "Processing event:" <<nev<<endl;
    detect->StartProcessEvent();
  } // event loop  
  delete detect;
}  
//__________________________________________________________________________________________________
void D_C()
{
  TStopwatch sw;TDatime time;

   AliRICHClusterFinder *z=new AliRICHClusterFinder(r); z->Exec();

   cout << endl;
   cout << "Info in Digits->Clusters: Start time: ";time.Print();
   cout << "Info in Digits->Clusters: Stop  time: ";time.Set();  time.Print();
   cout << "Info in Digits->Clusters: Time  used: ";sw.Print();
}
//__________________________________________________________________________________________________

void Show()
{  
  cout<<endl;
//load all trees  
  al->LoadHeader(); 
    al->LoadKinematics();  
      rl->LoadHits();  
        Bool_t isSdigits=!rl->LoadSDigits();  
          Bool_t isClusters=!rl->LoadRecPoints();
            Bool_t isDigits=!rl->LoadDigits();//loaders
  cout<<endl;  cout<<endl;  
  for(Int_t iEventN=0;iEventN<a->GetEventsPerRun();iEventN++){//events loop
    Int_t iNparticles=a->GetEvent(iEventN);
    Int_t iNprims=al->Stack()->GetNprimary();
    
    Int_t iTotalHits=0;
    for(Int_t iPrimN=0;iPrimN<iNprims;iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);      
      iTotalHits+=r->Hits()->GetEntries();
      TParticle *pPrim=al->Stack()->Particle(iPrimN);
      if(iPrimN<10) Info("Show","Evt %4i prim %4i has %4i hits from %s (,%7.2f,%7.2f)",
                  iEventN,iPrimN, r->Hits()->GetEntries(), pPrim->GetName(), pPrim->Theta()*r2d,pPrim->Phi()*r2d);
    }//prims loop
    Info("Show-HIT","Evt %i total:  %i particles %i primaries %i hits",
                        iEventN,   iNparticles, iNprims,     iTotalHits);
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


AliRun *a;
AliRunLoader *al;
AliLoader *rl,*tl,*il;

AliRICH *r;

Bool_t ReadAlice()
{
  Info("ReadAlice","Tring to read ALICE from SIMULATED FILE.");
  if(gAlice) delete gAlice;      
  if(!(al=AliRunLoader::Open("galice.root","AlicE","update"))){
    gSystem->Exec("rm -rf *.root *.dat");
    Error("menu.C::ReadAlice","galice.root broken, removing all this garbage then init new one");
    new AliRun("gAlice","Alice experiment system");
    gAlice->Init("Config.C");
    r=(AliRICH*)gAlice->GetDetector("RICH");
    return kFALSE;
  }
  al->LoadgAlice();
  if(!gAlice) Fatal("menu.C::ReadAlice","No gAlice in file");
  a=al->GetAliRun();
  a->SetDebug(0);    
//RICH      
  if(!(r=(AliRICH*)gAlice->GetDetector("RICH"))) Warning("RICH/menu.C::ReadAlice","No RICH in file");
  if(!(rl=al->GetLoader("RICHLoader")))          Warning("RICH/menu.C::ReadAlice","No RICH loader in file");        
        
  Info("ReadAlice","Run contains %i event(s)",gAlice->GetEventsPerRun());      
  return kTRUE;
}
//__________________________________________________________________________________________________
void TestResponse()
{
  TCanvas *pC=new TCanvas("c","Amplification test",900,800);
  pC->Divide(1,2);
  
  
  const Int_t nPoints=8;
  THStack *pStackPhot=new THStack("StackPhot","photons");
  THStack *pStackMip =new THStack("StackMip","mips");
  TLegend *pLeg=new TLegend(0.6,0.2,0.9,0.5,"legend");    
  TH1F *apHphot[nPoints];
  TH1F *apHmip[nPoints];
  
  Double_t starty=0;
  Double_t deltay=AliRICHParam::SectorSizeY()/nPoints;
  
  for(int i=0;i<nPoints;i++){
    apHphot[i]=new TH1F(Form("hphot%i",i),"Qdc for Photon;QDC;Counts",500,0,500); apHphot[i]->SetLineColor(i);pStackPhot->Add(apHphot[i]);                 
    apHmip[i] =new TH1F(Form("hmip%i",i),"Qdc for Mip;QDC;Counts",4000,0,4000);   apHmip[i]->SetLineColor(i);pStackMip->Add(apHmip[i]);                 
    
    pLeg->AddEntry(apHphot[i],Form("@(10,%5.2f->%5.2f)",starty+i*deltay,starty+i*deltay-AliRICHParam::SectorSizeY()/2));
  }
        
  
  TVector2 x2(0,0);  
//  AliRICHParam::ResetWireSag();
  for(Int_t i=0;i<10000;i++){//events loop
    for(int j=0;j<nPoints;j++){
      x2.Set(10,starty+j*deltay);
      apHphot[j]->Fill(AliRICHParam::TotQdc(x2,0));
      apHmip[j]->Fill(AliRICHParam::TotQdc(x2,gRandom->Landau(600,150)*1e-9));
    }
  }
  
  pC->cd(1);  pStackMip->Draw("nostack");
  pC->cd(2);  pStackPhot->Draw("nostack"); pLeg->Draw();
}//TestResponse()
//__________________________________________________________________________________________________
void TestSD()
{
  Info("TestSD","Creating test sdigits.");
  TVector3 hit(426.55,246.28,17.21);        
  TVector2 x2=r->C(4)->Glob2Loc(hit);        
  Int_t iTotQdc=r->P()->TotQdc(x2,624e-9);        
  Int_t iPadXmin,iPadXmax,iPadYmin,iPadYmax;
  Int_t padx,pady;
  Int_t sec=r->P()->Loc2Pad(x2,padx,pady);
  r->P()->Loc2Area(x2,iPadXmin,iPadYmin,iPadXmax,iPadYmax);
  Info("TestSD","Initial hit (%7.2f,%7.2f,%7.2f)->(%7.2f,%7.2f)->(%4i,%4i,%4i) gives %i charge",
                          hit.X(),hit.Y(),hit.Z(),x2.X(),x2.Y(),sec,padx,pady,iTotQdc);
  
  cout<<"left-down=("<<iPadXmin<<","<<iPadYmin<<") right-up=("<<iPadXmax<<','<<iPadYmax<<')'<<endl;
  for(Int_t iPadY=iPadYmin;iPadY<=iPadYmax;iPadY++)
    for(Int_t iPadX=iPadXmin;iPadX<=iPadXmax;iPadX++)
       cout<<r->P()->FracQdc(x2,iPadX,iPadY)<<endl;
  Info("TestSD","Stop.");
}//void TestSdigits()
//__________________________________________________________________________________________________
void TestC()
{
  Info("TestC","Creating test clusters.");
  R()->GetLoader()->MakeTree("R");R()->MakeBranch("R");
  
  AliRICHcluster c;
  c.AddDigit(new AliRICHdigit(1,20,21,200,1,2,3));
  c.AddDigit(new AliRICHdigit(1,22,21,250,1,2,3));
  c.CoG();
  
  R()->AddCluster(c);  
  
  R()->GetLoader()->TreeR()->Fill();  R()->GetLoader()->WriteRecPoints("OVERWRITE");
  R()->GetLoader()->UnloadRecPoints();
  R()->ResetClusters();
  
  Info("TestC","Stop.");
}//TestC()
//__________________________________________________________________________________________________
void TestSeg()
{
  AliRICHParam p;
  
  new TCanvas("name","PC segmentation");
  gPad->Range(-20,-20,200,150);
  AliRICHDisplFast::DrawSectors();
  
  TLatex t; t.SetTextSize(0.02);
  t.DrawText(0,140,"View from interaction point");
  t.DrawLatex(p.PcSizeX()+10,120,Form("Pc  %6.2fx%6.2fcm %3ix%3ipads",p.PcSizeX()    ,p.PcSizeY(),    p.NpadsX()   ,p.NpadsY()));
  t.DrawLatex(p.PcSizeX()+10,115,Form("Sec %6.2fx%5.2fcm %3ix%2ipads",p.SectorSizeX(),p.SectorSizeY(),p.NpadsXsec(),p.NpadsYsec()));
  t.DrawLatex(p.PcSizeX()+10,110,Form("Pad %6.2fx%4.2fcm DeadZone %6.2fcm",p.PadSizeX()   ,p.PadSizeY(),p.DeadZone()));
  
  TVector2 v2;
  t.SetTextAlign(12);
  v2=AliRICHParam::Pad2Loc( 1,24);  t.DrawText(v2.X(),v2.Y(),"sec 1");  
  v2=AliRICHParam::Pad2Loc(81,24);  t.DrawText(v2.X(),v2.Y(),"sec 2");  
  v2=AliRICHParam::Pad2Loc( 1,70);  t.DrawText(v2.X(),v2.Y(),"sec 3");  
  v2=AliRICHParam::Pad2Loc(81,70);  t.DrawText(v2.X(),v2.Y(),"sec 4");  
  v2=AliRICHParam::Pad2Loc( 1,120); t.DrawText(v2.X(),v2.Y(),"sec 5");  
  v2=AliRICHParam::Pad2Loc(81,120); t.DrawText(v2.X(),v2.Y(),"sec 6");  
  
//  TGaxis *pAx=new TGaxis(0,0,140,  0,0,140,510,"-="); pAx->SetTitle("x, cm"); pAx->SetTextSize(0.05); pAx->Draw();
//  TGaxis *pAy=new TGaxis(0,0,  0,140,0,140,510,"-="); pAy->SetTitle("y, cm"); pAy->SetTextSize(0.05); pAy->Draw();
   
  
  t.SetTextColor(kBlue);  

  Int_t padx,pady,sec;
  Double_t margin=5;
  Double_t x0=0; Double_t x1=p.SectorSizeX(); Double_t x2=p.SectorSizeX()+p.DeadZone(); Double_t x3=p.PcSizeX();
  Double_t y0=0; Double_t y1=p.SectorSizeY(); Double_t y2=p.SectorSizeY()+p.DeadZone(); 
  Double_t y3=2*p.SectorSizeY()+p.DeadZone(); Double_t y4=p.PcSizeY()-p.SectorSizeY();
  Double_t y5=p.PcSizeY();
   
//pads along x  
  t.SetTextAlign(11); t.DrawText(x0,y5+margin,"1");  
  t.SetTextAlign(31); t.DrawText(x1,y5+margin,"80");
  t.SetTextAlign(11); t.DrawText(x2,y5+margin,"81");   
  t.SetTextAlign(31); t.DrawText(x3,y5+margin,"160");
//pads along y    
  t.SetTextAlign(11); t.DrawText(x3+margin,y0,"1");  
  t.SetTextAlign(13); t.DrawText(x3+margin,y1,"48"); 
  t.SetTextAlign(11); t.DrawText(x3+margin,y2,"49");
  t.SetTextAlign(13); t.DrawText(x3+margin,y3,"96"); 
  t.SetTextAlign(11); t.DrawText(x3+margin,y4,"97");
  t.SetTextAlign(13); t.DrawText(x3+margin,y5,"144");        
  
  t.SetTextColor(kRed);
   
//positions along x   
  t.SetTextAlign(13);t.DrawText(x0,y0-margin,Form("%5.2f",x0));   
  t.SetTextAlign(33);t.DrawText(x1,y0-margin,Form("%5.2f",x1));   
  t.SetTextAlign(13);t.DrawText(x2,y0-margin,Form("%5.2f",x2));   
  t.SetTextAlign(33);t.DrawText(x3,y0-margin,Form("%5.2f",x3));   
//positions along y
  t.SetTextAlign(31);t.DrawText(x0-margin,y0,Form("%5.2f",y0));   
  t.SetTextAlign(33);t.DrawText(x0-margin,y1,Form("%5.2f",y1));   
  t.SetTextAlign(31);t.DrawText(x0-margin,y2,Form("%5.2f",y2));   
  t.SetTextAlign(33);t.DrawText(x0-margin,y3,Form("%5.2f",y3));   
  t.SetTextAlign(31);t.DrawText(x0-margin,y4,Form("%5.2f",y4));   
  t.SetTextAlign(33);t.DrawText(x0-margin,y5,Form("%5.2f",y5));   
//coners      
  t.SetTextSize(0.02);
  TVector pad(2);
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x0,y0)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x0,y1)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x1,y1)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x1,y0)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x0,y2)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x0,y3)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x1,y3)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x1,y2)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
   
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x0,y4)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x0,y5)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x1,y5)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x1,y4)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x2,y4)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x2,y5)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x3,y5)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x3,y4)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x2,y2)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x2,y3)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x3,y3)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x3,y2)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));

  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x2,y0)));t.SetTextAlign(11);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x2,y1)));t.SetTextAlign(13);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x3,y1)));t.SetTextAlign(33);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
  v2=p.Pad2Loc(p.Loc2Pad(TVector2(x3,y0)));t.SetTextAlign(31);t.DrawText(v2.X(),v2.Y(),Form("%5.2f,%5.2f",v2.X(),v2.Y()));
//new canvas
  new TCanvas("trasform","Test LRS-MRS transform");
  
  TView *pView=new TView(1);
  pView->SetRange(-600,-600,-600,600,600,600);
//axis  
  Double_t X[6]={0,0,0,300,0,0};  Double_t Y[6]={0,0,0,0,300,0};  Double_t Z[6]={0,0,0,0,0,300};  
  TPolyLine3D *pXaxis=new TPolyLine3D(2,X);pXaxis->SetLineColor(kRed);   pXaxis->Draw();
  TPolyLine3D *pYaxis=new TPolyLine3D(2,Y);pYaxis->SetLineColor(kGreen); pYaxis->Draw();
  TPolyLine3D *pZaxis=new TPolyLine3D(2,Z);pZaxis->SetLineColor(kBlue);  pZaxis->Draw();  
//Draw PC for all chambers by trasfering Pc plane using Pc2Mrs methode  
  Int_t iNpointsX=50,iNpointsY=50;  
  for(Int_t iChamberN=1;iChamberN<=7;iChamberN++){//chamber loop
    TPolyMarker3D *pChamber=new TPolyMarker3D(iNpointsX*iNpointsY);
    Int_t i=0;
    for(Double_t x=0;x<p.PcSizeX();x+=p.PcSizeX()/iNpointsX)
      for(Double_t y=0;y<p.PcSizeY();y+=p.PcSizeY()/iNpointsY){//step loop
        TVector3 v3=p.C(iChamberN)->Pc2Mrs(TVector2(x,y));//from regular grid of local PC points to MRS presentation
        pChamber->SetPoint(i++,v3.X(),v3.Y(),v3.Z());//Pc plane poing in MRS
      }//step loop
    pChamber->SetMarkerSize(1);
    pChamber->SetMarkerColor(iChamberN);
    pChamber->Draw();  
    t.SetNDC();t.SetTextColor(iChamberN); t.DrawText(0.1,iChamberN*0.1,Form("Chamber %i",iChamberN));    
  }//chamber loop   
//  gPad->GetView()->RotateView(94,45);
}//void TestSeg()
//__________________________________________________________________________________________________
void TestMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH test");
  pMenu->AddButton("Test segmentation",  "TestSeg()",         "Test AliRICHParam segmentation methods");
  pMenu->AddButton("Test response",      "TestResponse()",    "Test AliRICHParam response methods");
  pMenu->AddButton("Test sdigits",       "TestSD()",          "Create test set of sdigits");
  pMenu->AddButton("Test clusters",      "TestC()",           "Create test set of clusters");
  pMenu->Show();  
}//TestMenu()
//__________________________________________________________________________________________________
void menu()
{ 
  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  if(ReadAlice()){//it's from file, reconstruct
    pMenu->AddButton("Show",            "Show()",             "Shows the structure of events in files");
    pMenu->AddButton("Display Fast",    "DisplFast()",        "Display Fast");
    pMenu->AddButton("Control Plots",   "R()->ControlPlots()","Create some control histograms");
    
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
    pMenu->AddButton("Geo GUI", "new G3GeometryGUI;","Create instance of G4GeometryGUI"); 
    pMenu->AddButton("Read RAW","ReadRaw()","Read a list of digits from test beam file"); 
  }
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",         "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
  pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
  pMenu->AddButton("Quit",            ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//menu()
//__________________________________________________________________________________________________
void DebugOFF(){  Info("DebugOFF","");  a->SetDebug(0);  r->SetDebug(0);  AliLoader::SetDebug(0);}
void DebugON() {  Info("DebugON","");   a->SetDebug(1);  r->SetDebug(1);  AliLoader::SetDebug(1);}

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

  TH1F *distH1 = new TH1F("distH1","distH1",100,0.,5.);
  AliRICH *pRich = (AliRICH*)gAlice->GetDetector("RICH");
  Bool_t isHits    =!pRich->GetLoader()->LoadHits();

  pRich->GetLoader()->GetRunLoader()->LoadHeader();  pRich->GetLoader()->GetRunLoader()->LoadKinematics();


  if(!isHits){Error("Exec","No hits. No contribs to calculate.");return;}

  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events Loop
    pRich->GetLoader()->GetRunLoader()->GetEvent(iEventN);
    cout << " event n. " << iEventN << endl;
    Int_t nPrimaries = (Int_t)pRich->GetLoader()->TreeH()->GetEntries();
    TObjArray * Hits = new TObjArray[nPrimaries];

    for(Int_t i=0;i<nPrimaries;i++) {
      pRich->GetLoader()->TreeH()->GetEntry(i);
      Int_t nHits = pRich->Hits()->GetEntries();
      for(Int_t k=0;k<nHits;k++)         Hits[i].Add(pRich->Hits()->At(k));

    }
//deals with hits
      for(Int_t i=0;i<nPrimaries;i++){//prims loop
        pRich->GetLoader()->TreeH()->GetEntry(i);
        Int_t nHits = pRich->Hits()->GetEntries();
        for(Int_t j=0;j<nHits;j++){//hits loop
          AliRICHhit *pHit = (AliRICHhit*)Hits[i].At(j);
          TParticle *pParticle = pRich->GetLoader()->GetRunLoader()->GetAliRun()->Stack()->Particle(pHit->GetTrack());
//          if(pParticle->P()>1) FillContribs(pParticle->GetPdgCode(),pHit->C(),kFALSE);
          FillContribs(pParticle->GetPdgCode(),pHit->C(),kFALSE);
          if(pParticle->GetPDG()->Charge()!=0) distH1->Fill(pHit->Length());
        }//hits loop
      }//prims loop
    }// event loop
    FillContribs(0,0,kTRUE);
    distH1->Draw();

}//ParticleContribs()
