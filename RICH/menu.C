void Show()
{  
//How to get number of events:  
  Info("RICH/menu.C::Show","3 ways to get number of events 1=%i 2=%f 3=%f",
                          gAlice->GetEventsPerRun(),
                          al->TreeE()->GetEntries(),
                          gAlice->TreeE()->GetEntries());
  rl->LoadHits();
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){// loop on events
    Int_t iNparticles=gAlice->GetEvent(iEventN);
    Int_t iNtracks=rl->TreeH()->GetEntries();
    Info("RICH/menu.C::Show","Event %i contains %i particles in total while %i are primary",
                                     iEventN,    iNparticles,                iNtracks);
    for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){ // loop on tracks
      rl->TreeH()->GetEntry(iTrackN);
      al->Stack()->Particle(iTrackN)->Print();
      Info("RICH/menu.C::Show","track %i has %i hits",iTrackN,r->Hits()->GetEntries());
//      for(AliRICHhit *pRichHit=(AliRICHhit*)gpRich->FirstHit(-1);pRichHit;pRichHit=(AliRICHhit*)gpRich->NextHit()){
//        pRichHit->Dump();        
//              TVector3 armV3=pRICH->ToArm(mrsV3);
//              cout<<"After\n";armV3.Dump();
//      }//loop on hits of given track
    }// loop on tracks
  }//loop on events
}//void Show()

void menu(Int_t iNevents=5)// How many events to generate.
{ 
  Info("RICH/menu.C","%i event(s) are requested",iNevents);
   
  TString runString="gAlice->Run(";  runString=runString+iNevents+");";

  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  pMenu->AddButton("Debug ON-OFF",     "Debug();",   "Switch debug on-off");   
  if(GetAlice()){//it's from file, reconstruct
    pMenu->AddButton("Show","Show()","Shows the structure of events in files");
    pMenu->AddButton("Hits2SDigits","r->Hits2SDigits()","Perform first phase converstion");
    pMenu->AddButton("RingViewer","RingViewer()","Show rings with reconstructed info");
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",           runString.Data(),       "Process!");
  }
  pMenu->AddButton("Geo",          "Geo()",                 "Geomentry submenu");
  pMenu->AddButton("Browser",      "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Quit",         ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//void menu(Int_t iNevents)


void Debug()
{
  if(gAlice->GetDebug()){
    Info("RICH/menu.C::Debug","OFF");
    gAlice->SetDebug(0);
    if(r)r->SetDebug(0);
    if(t)t->SetDebug(0);
    AliLoader::SetDebug(0);
  }else{
    Info("RICH/menu.C::Debug","ON");
    gAlice->SetDebug(1);
    if(r)r->SetDebug(1);
    if(t)t->SetDebug(1);
    AliLoader::SetDebug(1);
  }
}//void Debug()

AliRun *a;
AliRunLoader *al;
AliLoader *rl,*tl,*il;

AliRICH *r;
AliTPC *t;
AliITS *i;
AliPHOS *p;

Bool_t GetAlice()
{
  if(gAlice){//it's aliroot
    Info("RICH/menu.C::GetAlice","gAlice!=NULL, IT'S ALIROOT, EXECUTE Init()");
    gAlice->Init("~/my/AliceConfig.C");
    r=(AliRICH*)gAlice->GetDetector("RICH");
    t=(AliTPC*)gAlice->GetDetector("TPC");
    i=(AliITS*)gAlice->GetDetector("ITS");
    p=(AliPHOS*)gAlice->GetDetector("PHOS");    
    return kFALSE;
  }else{//it's root with ALICE libs loaded
    Info("RICH/menu.C::GetAlice","gAlice=0, GETTING IT FROM SIMULATED FILE.");
        
    if(!(al=AliRunLoader::Open("galice.root","AlicE","update"))) Fatal("RICH/menu.C::GetAlice","Can't get AliRunLoader");
    al->LoadgAlice();
    if(!gAlice) Fatal("RICH/menu.C::GetAlice","No gAlice in file");
//    a=al->GetAliRun();    
    al->LoadHeader();//loads events tree
    al->LoadKinematics();//loads the primaries info
//RICH      
    if(!(r=(AliRICH*)gAlice->GetDetector("RICH"))) Warning("RICH/menu.C::GetAlice","No RICH in file");
    if(!(rl=al->GetLoader("RICHLoader")))          Warning("RICH/menu.C::GetAlice","No RICH loader in file");        
//TPC            
    if(!(t=(AliTPC*)gAlice->GetDetector("TPC")))   Warning("RICH/menu.C::GetAlice","No TPC in file");
    if(!(tl=al->GetLoader("TPCLoader")))           Warning("RICH/menu.C::GetAlice","No TPC loader in file");    
        
    Info("RICH/menu.C::GetAlice","Run contains %i event(s)",gAlice->GetEventsPerRun());    
    return kTRUE;
  }       
}//void GetAlice()         

void RingViewer()
{
  gStyle->SetPalette(1);
  TCanvas *view=new TCanvas("Display","ALICE RICH Display",0,0,600,600);
  
  TH2F *pH2=new TH2F("pH2F","RICH DISPLAY",r->Param()->Nx(),0,r->Param()->Nx(),r->Param()->Ny(),0,r->Param()->Ny());
  pH2->SetStats(0);
  pH2->SetMaximum(100);

  Int_t Nevents = gAlice->GetEventsPerRun();
}
//______________________________________________________________________________
void Geo()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH draw");
  pMenu->AddButton("Draw Isometry", "gMC->Gdraw(\"ALIC\", 60,120,0, 10,10, 0.01,0.01)","Draws ALIC volume in isometry");
  pMenu->AddButton("Draw Front XY", "gMC->Gdraw(\"ALIC\", 0,0,0, 10,10, 0.01,0.01)","Draws ALIC volume in XY view");
  pMenu->AddButton("Draw Side YZ",  "gMC->Gdraw(\"ALIC\",90,180, 0, 10,10, 0.01,0.01)","Draws ALIC volume in YZ view");
  pMenu->AddButton("Draw Top XZ",   "gMC->Gdraw(\"ALIC\",90, 90, 0, 10,10, 0.01,0.01)","Draws ALIC volume in XZ view");
  pMenu->AddButton("ALICE Tree", "((TGeant3*)gMC)->Gdtree(\"ALIC\")","Draws ALICE tree");      
  pMenu->AddButton("RICH Tree",  "((TGeant3*)gMC)->Gdtree(\"RICH\")","Draws RICH tree");      
  pMenu->AddButton("Geo test",  "GeoTest()",   "Draw RICH geo as a macro");
  pMenu->AddButton("Print ref", "PrintGeo()",  "Print RICH chambers default position");
  pMenu->AddButton("Print act", "r->Print();", "Print RICH chambers default position");
  pMenu->Show();  
}//void Draw()
//______________________________________________________________________________
void GeoTest()
{

  TBRIK *pAliceBRIK=new TBRIK("aliceBRIK","ALICE mother volume","void",500,500,500);
  TBRIK *pArmBRIK=new TBRIK("armBRIK","RICH arm1","void",pRICH->GetSizeX(),pRICH->GetSizeY(),pRICH->GetSizeZ());
   
   TNode *pAliceNode=new TNode("aliceNode","Mother volume","aliceBRIK");
   pAliceNode->cd();

// ARM 1 LEFT      
   TRotation rot1;
   TVector3  v1(0,pRICH->GetOffset(),0);
   
         
   rot1.RotateX(pRICH->GetYZAngle()*kDegrad);        v1.RotateX(pRICH->GetYZAngle()*kDegrad);
   rot1.RotateZ(pRICH->GetXYAngle()*kDegrad);        v1.RotateZ(pRICH->GetXYAngle()*kDegrad);
   rot1.RotateZ(pRICH->GetRotAngle()*kDegrad);       v1.RotateZ(pRICH->GetRotAngle()*kDegrad);
         
   TRotMatrix *pArm1RotMatrix=new TRotMatrix("rotArm1","rotArm1",rot1.ThetaX()*kRaddeg, rot1.PhiX()*kRaddeg,
                                                                 rot1.ThetaY()*kRaddeg, rot1.PhiY()*kRaddeg,
                                                                 rot1.ThetaZ()*kRaddeg, rot1.PhiZ()*kRaddeg);

   TNode *pArm1Node=new TNode("arm1Node","Left arm","armBRIK",v1.X(),v1.Y(),v1.Z(),"rotArm1");
   arm1Node->SetLineColor(kRed);
   
// ARM 2 LEFT      
   TRotation rot2;
   TVector3  v2(0,pRICH->GetOffset(),0);
   
         
   rot2.RotateX( pRICH->YZAngle()*kDegrad);        v2.RotateX(pRICH->GetYZAngle()*kDegrad);
   rot2.RotateZ(-pRICH->XYAngle()*kDegrad);        v2.RotateZ(-pRICH->GetXYAngle()*kDegrad);
   rot2.RotateZ( pRICH->RotAngle()*kDegrad);        v2.RotateZ(pRICH->GetRotAngle()*kDegrad);
         
   TRotMatrix *pArm2RotMatrix=new TRotMatrix("rotArm2","rotArm2",rot2.ThetaX()*kRaddeg, rot2.PhiX()*kRaddeg,
                                                                 rot2.ThetaY()*kRaddeg, rot2.PhiY()*kRaddeg,
                                                                 rot2.ThetaZ()*kRaddeg, rot2.PhiZ()*kRaddeg);

   TNode *pArm2Node=new TNode("arm2Node","Left arm","armBRIK",v2.X(),v2.Y(),v2.Z(),"rotArm2");
   arm2Node->SetLineColor(kBlue);
   
   aliceNode->Draw();
}//void GeoTest()
//______________________________________________________________________________
void PrintGeo(Float_t rotDeg=0)
{
  AliRICHParam *p=new AliRICHParam;  
  Double_t r=p->Offset();
  Double_t kP=p->AngleXY();
  Double_t kT=p->AngleYZ();
  Double_t kRot;
  
  if(rotDeg==0)
    kRot=p->AngleRot();
  else
    kRot=rotDeg*deg;
        
  cout<<endl;
  Double_t  phi=90*deg+kRot+kP,theta=90*deg+kT;
  Info("   menu for          1","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot+kP,theta=90*deg;
  Info("   menu for          2","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot,theta=90*deg-kT;    
  Info("   menu for          3","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  
  phi=90*deg+kRot,theta=90*deg;
  Info("   menu for          4","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));

  
  phi=90*deg+kRot,theta=90*deg+kT;
  Info("   menu for          5","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                                r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  
  phi=90*deg+kRot-kP,theta=90*deg;
  Info("   menu for          6","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot-kP,theta=90*deg+kT;
  Info("   menu for          7","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2d,  phi*r2d,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));

  delete p;
}//void PrintGeo()
