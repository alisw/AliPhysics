void SDigits2Digits()
{
  rl->LoadSDigits();
  rl->MakeTree("D");r->MakeBranch("D");
  
  rl->TreeS()->GetEntry(0);  
  r->Sdigits()->Sort();
  
  Int_t chamber,x,y,qdc,tr[3],id;
  Int_t kBad=-101;
  chamber=x=y=qdc=tr[0]=tr[1]=tr[2]=id=kBad;
  Int_t counter=kBad;//how many sdigits for a given pad
  Int_t start=0,end=0;
        
  for(Int_t i=0;i<r->Sdigits()->GetEntries();i++){//sdigits loop (sorted)
    AliRICHdigit *pSdig=(AliRICHdigit*)r->Sdigits()->At(i);
    if(pSdig->Id()==id){//still the same pad
      end++;
    }else{//new pad, add the pevious one
      if(id!=kBad) r->AddDigits(chamber,x,y,qdc,tr[0],tr[1],tr[2]);//cm-xpad-ypad-qdc-tr1-2-3
      chamber=pSdig->C();x=pSdig->X();y=pSdig->Y();qdc=pSdig->Qdc();tr[0]=pSdig->T(0);id=pSdig->Id();
      start=i;
    }
  }//sdigits loop (sorted)
  
  r->AddDigits(chamber,x,y,qdc,tr[0],tr[1],tr[2]);//add the last digit
        
  rl->TreeD()->Fill();  
  rl->WriteDigits("OVERWRITE");
    
  cout<<endl;
  r->Digits(1)->Print();
}

void Show3()
{  
  cout<<endl;
  al->LoadHeader();
  al->LoadKinematics();
  rl->LoadHits();
  
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events loop
    Int_t iNparticles=gAlice->GetEvent(iEventN);
    Int_t iNprims=rl->TreeH()->GetEntries();
    Info("Show3","Event %i contains %i particles in total while %i are primary",
                                     iEventN,    iNparticles,                iNprims);
    for(Int_t iPrimN=0;iPrimN<iNprims;iPrimN++){//prims loop
      rl->TreeH()->GetEntry(iPrimN);
      TParticle *pPrim=al->Stack()->Particle(iPrimN);
      Info("Show3","prim %i has %i hits and %i sdigits from %s (,%7.2f,%7.2f)",
                                        iPrimN,
                                        r->Hits()->GetEntries(),
                                        r->Specials()->GetEntries(),
                                                           pPrim->GetName(),
                                                           pPrim->Theta()*r2d,pPrim->Phi()*r2d);
//      for(AliRICHhit *pHit=r->FirstHit(-1);pHit;pHit=r->NextHit()){
//        pRichHit->Dump();        
//      }//loop on hits of given prim track
    }//prims loop
    if(!rl->LoadDigits()){
      Info("Show3","Event %i contains %i digits",iEventN,r->Digits(1)->GetEntries());
    }
  }//events loop
  rl->UnloadHits();
  rl->UnloadDigits();
  al->UnloadHeader();
  al->UnloadKinematics();
  cout<<endl;
}//void Show()

void menu()// How many events to generate.
{ 
   

  TControlBar *pMenu = new TControlBar("vertical","RICH main");
       
  pMenu->AddButton("Debug ON",     "DebugON();",   "Switch debug on-off");   
  pMenu->AddButton("Debug OFF",    "DebugOFF();",   "Switch debug on-off");   
  if(CheckAlice()){//it's from file, reconstruct
    pMenu->AddButton("Show3","Show3()","Shows the structure of events in files");
    pMenu->AddButton("Hits2SDigits","r->Hits2SDigits()","Perform first phase converstion");
    pMenu->AddButton("SDigits2Digits","SDigits2Digits()","Perform second phase converstion");
    pMenu->AddButton("RingViewer","RingViewer()","Show rings with reconstructed info");
  }else{//it's aliroot, simulate
    pMenu->AddButton("Run",         "a->Run(1)",       "Process!");
  }
  pMenu->AddButton("Geo submenu",          "Geo()",            "Shows geomentry submenu");
  pMenu->AddButton("Test submenu",    "TestMenu()",            "Shows test submenu");
  pMenu->AddButton("Browser",      "new TBrowser;",         "Start ROOT TBrowser");
  pMenu->AddButton("Quit",         ".q",                    "Close session");
  pMenu->Show();
  a=gAlice;//for manual manipulation convinience
}//void menu(Int_t iNevents)


void DebugOFF()
{
  Info("DebugOFF","");
  a->SetDebug(0);
  r->SetDebug(0);
  AliLoader::SetDebug(0);
}//void DebugOFF()

void DebugON()
{
  Info("DebugON","");  
  a->SetDebug(1);
  r->SetDebug(1);
  AliLoader::SetDebug(1);
}//void DebugON()


AliRun *a;
AliRunLoader *al;
AliLoader *rl,*tl,*il;

AliRICH *r;

Bool_t CheckAlice()
{
  if(gAlice){//it's aliroot
    if(gSystem->Exec("ls galice.root")){
      Info("CheckAlice","It's AliRoot, and no galice.root: SIMULATION");
      gAlice->Init("AliceConfig.C");
      r=(AliRICH*)gAlice->GetDetector("RICH");
      return kFALSE;
    }else{//galice.root is present we want to read alice from file
      ReadAlice();
      return kTRUE;
    }
  }else{//it's root with ALICE libs loaded
    ReadAlice();
    return kTRUE;
  }       
}//void CheckAlice()         

void ReadAlice()
{
  Info("ReadAlice","Reading ALICE from SIMULATED FILE.");
  AliLoader::SetDebug(0);
  if(gAlice) delete gAlice;      
  if(!(al=AliRunLoader::Open("galice.root","AlicE","update"))){
    gSystem->Exec("rm -rf *.root *.dat");
    Fatal("ReadAlice","galice.root broken, removing all this garbage");
  }
  al->LoadgAlice();
  if(!gAlice) Fatal("ReadAlice","No gAlice in file");
  a=al->GetAliRun();
  a->SetDebug(0);    
//RICH      
  if(!(r=(AliRICH*)gAlice->GetDetector("RICH"))) Warning("RICH/menu.C::ReadAlice","No RICH in file");
  r->SetDebug(0);
  if(!(rl=al->GetLoader("RICHLoader")))          Warning("RICH/menu.C::ReadAlice","No RICH loader in file");        
        
  Info("ReadAlice","Run contains %i event(s)",gAlice->GetEventsPerRun());      
}
//__________________________________________________________________________________________________
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
  pMenu->AddButton("RICH Isometry", "gMC->Gdraw(\"ALIC\", 60,120,0, 10,10, 0.01,0.01)","Draws ALIC volume in isometry");
  pMenu->AddButton("RICH Front XY", "gMC->Gdraw(\"ALIC\", 0,0,0, 10,10, 0.01,0.01)","Draws ALIC volume in XY view");
  pMenu->AddButton("RICH Side YZ",  "gMC->Gdraw(\"ALIC\",90,180, 0, 10,10, 0.01,0.01)","Draws ALIC volume in YZ view");
  pMenu->AddButton("RICH Top XZ",   "gMC->Gdraw(\"ALIC\",90, 90, 0, 10,10, 0.01,0.01)","Draws ALIC volume in XZ view");
  pMenu->AddButton("Module Isometry","gMC->Gdraw(\"SRIC\", 30,60,0, 10,10, 0.1,0.1)","Draws SRIC volume in isometry");
  pMenu->AddButton("Module Front XY","gMC->Gdraw(\"SRIC\", 0,0,0, 10,10, 0.1,0.1)","Draws SRIC volume in XY view");
  pMenu->AddButton("Module Top XZ", "gMC->Gdraw(\"SRIC\",90, 90, 0, 10,10, 0.1,0.1)","Draws SRIC volume in XZ view");
  pMenu->AddButton("ALICE Tree", "((TGeant3*)gMC)->Gdtree(\"ALIC\")","Draws ALICE tree");      
  pMenu->AddButton("RICH Tree",  "((TGeant3*)gMC)->Gdtree(\"RICH\")","Draws RICH tree");      
  pMenu->AddButton("Geo test",  "GeoTest()",   "Draw RICH geo as a macro");
  pMenu->AddButton("Print ref", "PrintGeo()",  "Print RICH chambers default position");
  pMenu->AddButton("AliRICH::Print", "r->Print();", "Print RICH chambers default position");
  pMenu->AddButton("Test transform","TestTransform()","Test L2G and G2L methods");
  pMenu->AddButton("Geo GUI", "new G3GeometryGUI;","Create instance of G4GeometryGUI"); 
  pMenu->Show();  
}//void Draw()
//__________________________________________________________________________________________________
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
//__________________________________________________________________________________________________
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



void TestMenu()
{
  TControlBar *pMenu = new TControlBar("vertical","RICH test");
  pMenu->AddButton("Test segmentation",  "TestSegmentation()","Test AliRICHParam::L2P() method");
  pMenu->AddButton("Test transform",     "TestTransform()",   "Test ALiRICHChamber::L2G() and G2L methods");
  pMenu->AddButton("Test gain",          "TestGain()",        "Test AliRICHParam::Gain() method");
  pMenu->AddButton("Test MIP charge",    "TestMipCharge()",   "Test AliRICHParam::TotalCharge() method");
  pMenu->AddButton("Test Sdigits",       "TestSdigits()",     "Create test set of sdigits");
  pMenu->AddButton("Test Digits",        "TestDigits()",      "Create test set of digits");
  pMenu->Show();  
}//void Draw()


void TestGain()
{
  AliRICHParam *pParam=new AliRICHParam;
  AliRICHResponseV0 *pRes=new AliRICHResponseV0;
  
  TLegend *pLegend=new TLegend(0.6,0.3,0.85,0.5);  
  TH1F *pH0=new TH1F("pH1","Gain",100,0,600); 
  TH1F *pH10=new TH1F("pH10","Gain",100,0,600);
  TH1F *pH20=new TH1F("pH20","Gain",100,0,600);
  TH1F *pH30=new TH1F("pH30","Gain",100,0,600);
  TH1F *pHold=new TH1F("pHold","Mip Charge",100,0,2000);
  for(int i=0;i<1000;i++){
    pH0 ->Fill(pParam->Gain(0));
    pH10->Fill(pParam->Gain(10));
    pH20->Fill(pParam->Gain(20));
    pH30->Fill(pParam->Gain(30));
    pHold->Fill(pRes->IntPH(30));
  }
  pH0->Draw();
  pH10->Draw("same");
  pH20->Draw("same");
  pH30->Draw("same");
  pHold->Draw("same");
  pLegend->AddEntry(pH0,"y=0");  
  pLegend->AddEntry(pH10,"y=10");pH10->SetLineColor(kRed);
  pLegend->AddEntry(pH20,"y=20");pH20->SetLineColor(kBlue);
  pLegend->AddEntry(pH30,"y=30");pH30->SetLineColor(kGreen);
  pLegend->AddEntry(pHold,"res");pHold->SetLineColor(kMagenta);
  pLegend->Draw();
}//void TestGain()
//__________________________________________________________________________________________________
void TestMipCharge()
{
  AliRICHParam *pParam=new AliRICHParam;
  AliRICHResponseV0 *pRes=new AliRICHResponseV0;
  
  TLegend *pLegend=new TLegend(0.6,0.3,0.85,0.5);  
  TH1F *pH0= new TH1F("pH1", "Mip Charge",100,0,500); 
  TH1F *pH10=new TH1F("pH10","Mip Charge",100,0,500);
  TH1F *pH20=new TH1F("pH20","Mip Charge",100,0,500);
  TH1F *pH30=new TH1F("pH30","Mip Charge",100,0,500);
  TH1F *pHold=new TH1F("pHold","Mip Charge",100,0,500);
  for(int i=0;i<1000;i++){
    pH0 ->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,0));
    pH10->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,10));
    pH20->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,20));
    pH30->Fill(pParam->TotalCharge(kPiPlus,0.5e-9,30));
    pHold->Fill(pRes->IntPH(0.5e-9,-30));
  }
  pH0->Draw();
  pH10->Draw("same");
  pH20->Draw("same");
  pH30->Draw("same");
  pHold->Draw("same");
  pLegend->AddEntry(pH0,"y=0");  
  pLegend->AddEntry(pH10,"y=10");pH10->SetLineColor(kRed);
  pLegend->AddEntry(pH20,"y=20");pH20->SetLineColor(kBlue);
  pLegend->AddEntry(pH30,"y=30");pH30->SetLineColor(kGreen);
  pLegend->AddEntry(pHold,"res");pHold->SetLineColor(kMagenta);
  pLegend->Draw();
}//void TestGain()
//__________________________________________________________________________________________________
void TestDigits()
{
  rl->MakeTree("D");r->MakeBranch("D");
  
  r->AddDigits(1,10,11,25,0);
  
  r->AddDigits(1,25,30,25,10,12);
  r->AddDigits(1,26,30,25,10,12);
  
  rl->TreeD()->Fill();
  rl->WriteDigits("OVERWRITE");
  rl->UnloadDigits();
  r->ResetDigits();
}//void TestDigits()
//__________________________________________________________________________________________________
void TestSdigits()
{
  rl->MakeTree("S");r->MakeBranch("S");
//totally 19 must be trasformd to 6 digits
  r->AddSdigit(1,40,40,10,40); r->AddSdigit(1,40,40,10,41); r->AddSdigit(1,40,40,10,42); r->AddSdigit(1,40,40,10,43);
  r->AddSdigit(1,45,45,10,45); r->AddSdigit(1,45,45,10,45); r->AddSdigit(1,45,45,10,45); r->AddSdigit(1,45,45,10,45);
  
  r->AddSdigit(1,20,20,10,20);  r->AddSdigit(1,20,20,10,21);
  r->AddSdigit(1,25,25,10,25);  r->AddSdigit(1,25,25,10,25);
  
  r->AddSdigit(1,30,30,10,30); r->AddSdigit(1,30,30,10,31); r->AddSdigit(1,30,30,10,32);
  r->AddSdigit(1,35,35,10,35); r->AddSdigit(1,35,35,10,35); r->AddSdigit(1,35,35,10,35);
  

  r->AddSdigit(1,10,10,10,10);
  
//N2
//   r->AddSdigit(2,23,43,8,0);
//   
//   r->AddSdigit(2,22,42,8,0);
//   r->AddSdigit(2,23,42,8,0);
//   
//   r->AddSdigit(2,18,41,8,0);
//   r->AddSdigit(2,25,41,8,0);
//   
//   r->AddSdigit(2,17,40,8,0);
//   
//   r->AddSdigit(2,27,38,8,0);
//   r->AddSdigit(2,28,38,8,0);
//   
//   r->AddSdigit(2,16,37,8,0);
//   r->AddSdigit(2,21,37,8,0);
//   r->AddSdigit(2,22,37,8,0);
//   
//   r->AddSdigit(2,16,36,8,0);
//   r->AddSdigit(2,21,36,8,0);
//   
//   r->AddSdigit(2,22,35,8,0);
//   r->AddSdigit(2,28,35,8,0);
//   
//   r->AddSdigit(2,16,34,8,0);
//   
//   r->AddSdigit(2,18,32,8,0);
//   r->AddSdigit(2,25,32,8,0);
//   r->AddSdigit(2,27,32,8,0);
//   
//   r->AddSdigit(2,19,31,8,0);
//   r->AddSdigit(2,23,31,8,0);
//   r->AddSdigit(2,24,31,8,0);
//   
//   r->AddSdigit(2,19,30,8,0);
//   r->AddSdigit(2,20,30,8,0);
  
  rl->TreeS()->Fill();
  rl->WriteSDigits("OVERWRITE");
  rl->UnloadSDigits();
  cout<<endl;r->Sdigits()->Print();
  r->ResetSDigits();
}//void TestDigits()
//__________________________________________________________________________________________________
void TestSegmentation()
{
  AliRICHParam *p=r->Param(); 
  Float_t dz=p->DeadZone();
  Float_t sx=p->SectorSizeX();Float_t sy=p->SectorSizeY();
  Float_t px=p->PcSizeX();    Float_t py=p->PcSizeY();
  Int_t padx,pady;
  cout<<endl;
  Info("  1-  1","sec=%i padx=%3i pady=%3i",p->L2P(-px/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 48-  1","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2-dz , -py/2        ,padx,pady),padx,pady);
  Info(" 49-  1","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 96-  1","sec=%i padx=%3i pady=%3i",p->L2P( sx/2    , -py/2        ,padx,pady),padx,pady);
  Info(" 97-  1","sec=%i padx=%3i pady=%3i",p->L2P( sx/2+dz , -py/2        ,padx,pady),padx,pady);
  Info("144-  1","sec=%i padx=%3i pady=%3i",p->L2P( px/2    , -py/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1- 80","sec=%i padx=%3i pady=%3i",p->L2P(-px/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 48- 80","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2-dz , -dz/2        ,padx,pady),padx,pady);
  Info(" 49- 80","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 96- 80","sec=%i padx=%3i pady=%3i",p->L2P( sx/2    , -dz/2        ,padx,pady),padx,pady);
  Info(" 97- 80","sec=%i padx=%3i pady=%3i",p->L2P( sx/2+dz , -dz/2        ,padx,pady),padx,pady);
  Info("144- 80","sec=%i padx=%3i pady=%3i",p->L2P( px/2    , -dz/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1- 81","sec=%i padx=%3i pady=%3i",p->L2P(-px/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 48- 81","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2-dz ,  dz/2        ,padx,pady),padx,pady);
  Info(" 49- 81","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 96- 81","sec=%i padx=%3i pady=%3i",p->L2P( sx/2    ,  dz/2        ,padx,pady),padx,pady);
  Info(" 97- 81","sec=%i padx=%3i pady=%3i",p->L2P( sx/2+dz ,  dz/2        ,padx,pady),padx,pady);
  Info("144- 81","sec=%i padx=%3i pady=%3i",p->L2P( px/2    ,  dz/2        ,padx,pady),padx,pady);
  cout<<endl;
  Info("  1-160","sec=%i padx=%3i pady=%3i",p->L2P(-px/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 48-160","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2-dz ,  py/2        ,padx,pady),padx,pady);
  Info(" 49-160","sec=%i padx=%3i pady=%3i",p->L2P(-sx/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 96-160","sec=%i padx=%3i pady=%3i",p->L2P( sx/2    ,  py/2        ,padx,pady),padx,pady);
  Info(" 97-160","sec=%i padx=%3i pady=%3i",p->L2P( sx/2+dz ,  py/2        ,padx,pady),padx,pady);
  Info("144-160","sec=%i padx=%3i pady=%3i",p->L2P( px/2    ,  py/2        ,padx,pady),padx,pady);  
  cout<<endl;
  Info(" 73-160","sec=%i padx=%3i pady=%3i",p->L2P(    0    ,  py/2      ,padx,pady),padx,pady);    
  Info(" 73- 81","sec=%i padx=%3i pady=%3i",p->L2P(    0    ,  dz/2      ,padx,pady),padx,pady);    
  Info(" 73- 80","sec=%i padx=%3i pady=%3i",p->L2P(    0    ,   0        ,padx,pady),padx,pady);    
  Info("144- 81","sec=%i padx=%3i pady=%3i",p->L2P(    0    , -dz/2      ,padx,pady),padx,pady);    
  Info("144- 81","sec=%i padx=%3i pady=%3i",p->L2P(    0    , -py/2      ,padx,pady),padx,pady);    
}//void TestSegmentation()
//__________________________________________________________________________________________________
void TestClusters()
{
  rl->MakeTree("R");r->MakeBranch("R");
  r->AddCluster  
  rl->TreeR()->Fill();
  rl->WriteRecPoints("OVERWRITE");
  rl->UnloadRecPoints();
  r->ResetClusters();
}//void TestClusters()
