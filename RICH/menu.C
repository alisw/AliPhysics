void Show()
{  
//How to get number of events:  
  gpRL->LoadHeader();//without this the following will be zero:
  Info("RICH/menu.C::Show","3 ways to get number of events 1=%i 2=%i 3=%i",
                          gAlice->GetEventsPerRun(),
                          gpRL->TreeE()->GetEntries(),
                          gAlice->TreeE()->GetEntries());
//
  return;
  for(Int_t iEventN=0;iEventN<iNevents;iEventN++){// loop on events
      Int_t iNparticles=gAlice->GetEvent(iEventN);
      Int_t iNtracks=gAlice->TreeH()->GetEntries();
      cout<<"Event "<<iEventN<<" contains "<<iNparticles<<" particles and "<<iNtracks<<" tracks\n";

      for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){ // loop on tracks
         gAlice->TreeH()->GetEntry(iTrackN);
        // gAlice->Particle(iTrackN)->Print();
         for(AliRICHhit *pRICHhit=(AliRICHhit*)pRICH->FirstHit(-1);pRICHhit;pRICHhit->(AliRICHhit*)pRICH->NextHit()){
             TVector3 mrsV3(pRICHhit->X(),pRICHhit->Y(),pRICHhit->Z());
             cout<<"Before\n";mrsV3.Dump();
             TVector3 armV3=pRICH->ToArm(mrsV3);
             cout<<"After\n";armV3.Dump();
         }//loop on hits of given track
      }// loop on tracks
   }//loop on events

}//void Show()

void menu(Int_t iNevents=5)// How many events to generate.
{ 
  Info("RICH/menu.C","%i event(s) are requested",iNevents);
   
  TString runString="gAlice->Run(";  runString=runString+iNevents+");";

  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  if(GetAlice()){//it's from file, reconstruct
    pMenu->AddButton("Show","Show()","Shows the structure of events in files");
    pMenu->AddButton("RingViewer","RingViewer()","Show rings with reconstructed info");
  }else{//it's aliroot, simulate
    pMenu->AddButton("Init",         "Init()",                "Normal init");   
    pMenu->AddButton("Run",           runString.Data(),       "Process!");
  }
  pMenu->AddButton("Debug ON",     "gAlice->SetDebug(1)",   "Switch debug on");   
  pMenu->AddButton("Debug OFF",    "gAlice->SetDebug()",    "Switch debug off");     
  pMenu->AddButton("Geo",          "Geo()",                 "Geomentry submenu");
  pMenu->AddButton("Browser",      "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Quit",         ".q",                    "Close session");
  pMenu->Show();
}//void menu(Int_t iNevents)

void Init()
{
  gAlice->Init("~/my/AliceConfig.C");
}

AliRunLoader *rl,*gpRL;
AliLoader *rrl,*gpRichLoader;
AliRICH *r,*gpRich;
AliRun *a;
Bool_t GetAlice()
{
  if(gAlice){//it's aliroot
    Info("RICH/menu.C::GetAlice","gAlice!=NULL, it's aliroot, execute Init()");
    Init();
    return kFALSE;
  }else{//it's root with ALICE libs loaded
    Info("RICH/menu.C::GetAlice","gAlice=0, getting it from simulated file.");
    
    gpRL=AliRunLoader::Open("galice.root","AlicE","update");
    if(!gpRL) Fatal("GetAlice","Can't get AliRunLoader");
    
    gpRL->LoadgAlice();
    if(!gAlice) Fatal("RICH/menu.C::GetAlice","No gAlice in file");
//    gAlice=gpRL->GetAliRun();    
    
    gpRich=(AliRICH*)gAlice->GetDetector("RICH"); 
    if(!gpRich) Fatal("RICH/menu.C::GetAlice","No RICH in file");
    
    gpRichLoader=gpRL->GetLoader("RICHLoader");     
    if(!gpRichLoader) Fatal("RICH/menu.C::GetAlice","No RICH loader in file");    
    
    Info("RICH/menu.C::GetAlice","contains %i event(s)",gAlice->GetEventsPerRun());
    rl=gpRL;rrl=gpRichLoader;r=gpRich;a=gAlice;//for manual manipulation convinience
    return kTRUE;
  }       
}//void GetAlice()         



void Geo()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH draw");
  pMenu->AddButton("Draw Front", "gMC->Gdraw(\"ALIC\", 0,0,0, 10,10, 0.01,0.01)","Draws ALIC volume in XY view");
  pMenu->AddButton("Draw Side",  "gMC->Gdraw(\"ALIC\",90,180, 0, 10,10, 0.01,0.01)","Draws ALIC volume in YZ view");
  pMenu->AddButton("Draw Top",   "gMC->Gdraw(\"ALIC\",90, 90, 0, 10,10, 0.03,0.03)","Draws ALIC volume in XZ view");
  pMenu->AddButton("ALICE Tree", "((TGeant3*)gMC)->Gdtree(\"ALIC\")","Draws ALICE tree");      
  pMenu->AddButton("RICH Tree",  "((TGeant3*)gMC)->Gdtree(\"RICH\")","Draws RICH tree");      
  pMenu->AddButton("Geo test",  "GeoTest()",   "Draw RICH geo as a macro");
  pMenu->AddButton("Print ref", "PrintGeo()",  "Print RICH chambers default position");
  pMenu->AddButton("Print act", "r->Print();", "Print RICH chambers default position");
  pMenu->Show();  
}//void Draw()

void GeoTest()
{

   TBRIK *pAliceBRIK=new TBRIK("aliceBRIK","ALICE mother volume","void",500,500,500);
   AliRICH *pRICH=(AliRICH*)gAlice->GetDetector("RICH");
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
   
         
   rot2.RotateX(pRICH->GetYZAngle()*kDegrad);        v2.RotateX(pRICH->GetYZAngle()*kDegrad);
   rot2.RotateZ(-pRICH->GetXYAngle()*kDegrad);        v2.RotateZ(-pRICH->GetXYAngle()*kDegrad);
   rot2.RotateZ(pRICH->GetRotAngle()*kDegrad);        v2.RotateZ(pRICH->GetRotAngle()*kDegrad);
         
   TRotMatrix *pArm2RotMatrix=new TRotMatrix("rotArm2","rotArm2",rot2.ThetaX()*kRaddeg, rot2.PhiX()*kRaddeg,
                                                                 rot2.ThetaY()*kRaddeg, rot2.PhiY()*kRaddeg,
                                                                 rot2.ThetaZ()*kRaddeg, rot2.PhiZ()*kRaddeg);

   TNode *pArm2Node=new TNode("arm2Node","Left arm","armBRIK",v2.X(),v2.Y(),v2.Z(),"rotArm2");
   arm2Node->SetLineColor(kBlue);
   
   aliceNode->Draw();
}//void GeoTest()

void RingViewer()
{
  gStyle->SetPalette(1);
  TCanvas *view = new TCanvas("Display","ALICE RICH Display",0,0,1200,750);
  
  TH2F *p2F = new TH2F("p2F","RICH DISPLAY",160,0,160,144,0,144);
  p2F->SetStats(0);
  p2F->SetMaximum(100);

  Int_t Nevents = gAlice->GetEventsPerRun();
}

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
  Double_t phi=90*deg+kRot,theta=90*deg-kT;    
  Info("   menu for          0","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot+kP,theta=90*deg;
  Info("   menu for          1","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot,theta=90*deg;
  Info("   menu for          2","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));

  
  phi=90*deg+kRot-kP,theta=90*deg;
  Info("   menu for          3","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));


  phi=90*deg+kRot+kP,theta=90*deg+kT;
  Info("   menu for          4","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  
  phi=90*deg+kRot,theta=90*deg+kT;
  Info("   menu for          5","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                                r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));
  
  phi=90*deg+kRot-kP,theta=90*deg+kT;
  Info("   menu for          6","r=%8.3f theta=%5.1f phi=%5.1f x=%8.3f y=%8.3f z=%8.3f",
                                 r,      theta*r2g,  phi*r2g,  
                                                               r*sin(theta)*cos(phi),
                                                                       r*sin(theta)*sin(phi),
                                                                               r*cos(theta));

  delete p;
}//void PrintGeo()
