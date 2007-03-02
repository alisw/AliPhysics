//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsd(AliESD *pEsd)
{
  TParticle part; TLorentzVector mom;
  for(Int_t iTrk=0;iTrk<25;iTrk++){//stack loop
    part.SetPdgCode(kProton);
    part.SetProductionVertex(0,0,0,0);  
    Double_t eta= -0.4+gRandom->Rndm()*0.8; //rapidity is random [-0.4,+0.4]
    Double_t phi= gRandom->Rndm()*1.4;      //phi is random      [ 0  , 80 ] degrees    
    mom.SetPtEtaPhiM(2,eta,phi,part.GetMass());
    part.SetMomentum(mom);
    AliESDtrack trk(&part);
    pEsd->AddTrack(&trk);
  }//stack loop  
}//EsdFromStack()
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

    Float_t th,ph,xPc,yPc,; pTrk->GetHMPIDtrk(xPc,yPc,th,ph); rec.SetTrack(xRa,yRa,th,ph); 
    
    if(!AliHMPIDDigit::IsInDead(xPc,yPc)) new((*pHitLst)[hc++]) AliHMPIDHit(ch,200e-9,kProton  ,iTrk,xPc,yPc);                 //mip hit
    for(int i=0;i<4;i++)  new((*pHitLst)[hc++]) AliHMPIDHit(ch,7.5e-9,kFeedback,iTrk,gRandom->Rndm()*130,gRandom->Rndm()*126); //bkg hits 4 per track
    for(int i=0;i<16;i++){
      rec.TracePhot(ckov,gRandom->Rndm()*TMath::TwoPi(),pos);
      new((*pHitLst)[hc++]) AliHMPIDHit(ch,7.5e-9,kCerenkov,iTrk,pos.X(),pos.Y());
    }                      //photon hits  
  }//tracks loop    
}//HitsFromEsd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawSeg()
{
  TCanvas *pC=new TCanvas("seg","Segmentation as seen from electronics side");
  DrawPc(1);

  pC->ToggleEventStatus();
  pC->SetEditable(0);
  pC->AddExec("status","DrawStatus()");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawStatus()
{// Show info about current cursur position in status bar of the canvas
//  Printf("event %i",gPad->GetEvent()); return;
  TPad *pad=gPad;
  TCanvas *pC=(TCanvas*)pad; 
  TRootCanvas *pRC= (TRootCanvas*)pC->GetCanvasImp();
  TGStatusBar *pBar=pRC->GetStatusBar();
  pBar->SetParts(5);
  Float_t x=pad->AbsPixeltoX(pad->GetEventX()); Float_t y=pad->AbsPixeltoY(pad->GetEventY());
  if(AliHMPIDDigit::IsInDead(x,y))
    pBar->SetText("Out of sensitive area",4);    
  else{
    Int_t pc,px,py,w32,ddl,r,d,a;  AliHMPIDDigit::Lors2Pad(x,y,pc,px,py); AliHMPIDDigit dig; dig.Set(1,pc,px,py); dig.Raw(w32,ddl,r,d,a);
    pBar->SetText(Form("(pc%i,px%i,py%i) (r%i,d%i,a%i) (%.2f,%.2f)",
                         pc  ,px  ,py,    r  ,d  ,a   ,dig.LorsX(),dig.LorsY()),4);
  }    
//  if(pad->GetEvent()==1){
//    new TCanvas("zoom",Form("Row %i DILOGIC %i",dig.Row(),dig.Dilogic()));  
//  }
}//DrawStatus()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawPc(Bool_t isFill) 
{ 
  gPad->Range(-10,-10,AliHMPIDDigit::SizeAllX()+5,AliHMPIDDigit::SizeAllY()+5); 
  
    
  Float_t dX2=0.5*AliHMPIDDigit::SizePadX(),
          dY2=0.5*AliHMPIDDigit::SizePadY() ;
  
  TLatex txt; txt.SetTextSize(0.01);
  TLine *pL;
  
  AliHMPIDDigit dig;   UInt_t w32; Int_t ddl,r,d,a;    
  
  for(Int_t iPc=AliHMPIDDigit::kMinPc;iPc<=AliHMPIDDigit::kMaxPc;iPc++){
    TBox *pBox=new TBox(AliHMPIDDigit::fMinPcX[iPc],AliHMPIDDigit::fMinPcY[iPc],
                        AliHMPIDDigit::fMaxPcX[iPc],AliHMPIDDigit::fMaxPcY[iPc]);
    
    if(iPc==0||iPc==2||iPc==4) pBox->SetFillColor(29);
    else                       pBox->SetFillColor(41);
    pBox->Draw();
    
    if(!isFill)  continue;
    
//    if(iPc%2) {dig.Set(0,iPc,79,25); txt.DrawText(dig.LorsX()+2,dig.LorsY(),Form("PC%i",dig.Pc()));}//print PC#    

    txt.SetTextAlign(32);
    for(Int_t iRow=0;iRow<8 ;iRow++){//draw row lines (horizontal)
      dig.Set(0,iPc,0,iRow*6); dig.Raw(w32,ddl,r,d,a);  //set digit to the left-down pad of this row
      
      if(iPc%2) txt.DrawText(dig.LorsX()-1           ,dig.LorsY(),Form("%i",dig.PadPcY()));                  //print PadY#    
                txt.DrawText(dig.LorsX()-1+(iPc%2)*67,dig.LorsY()+2,Form("r%i",r));                          //print Row#    
      pL=new TLine(dig.LorsX()-dX2,dig.LorsY()-dY2,dig.LorsX()+AliHMPIDDigit::SizePcX()-dX2,dig.LorsY()-dY2);//draw horizontal line 
      if(iRow!=0) pL->Draw(); 
    }//row loop  
    
    txt.SetTextAlign(13);
    for(Int_t iDil=0;iDil<10;iDil++){//draw dilogic lines (vertical)
      dig.Set(0,iPc,iDil*8,0); dig.Raw(w32,ddl,r,d,a);       //set this digit to the left-down pad of this dilogic        
      
                           txt.DrawText(dig.LorsX()  ,dig.LorsY()-1,Form("%i",dig.PadPcX()));                 //print PadX# 
      if(iPc==4 || iPc==5) txt.DrawText(dig.LorsX()+2,dig.LorsY()+42,Form("d%i",d));              //print Dilogic#    
      pL=new TLine(dig.LorsX()-dX2,dig.LorsY()-dY2,dig.LorsX()-dX2,dig.LorsY()+AliHMPIDDigit::SizePcY()-dY2); //draw vertical line
      if(iDil!=0)pL->Draw();
    }//dilogic loop        
  }//PC loop      
}//DrawPc()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Close_hed()
{
  TCanvas *pC = ((TCanvas*)gROOT->FindObject("hed"));if(!pC) return;
  pC->Close();
  pC=0x0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ReadHits(TClonesArray *pHitLst)
{
  pHitLst->Delete();  Int_t cnt=0;
  for(Int_t iEnt=0;iEnt<hl->TreeH()->GetEntries();iEnt++){       //TreeH loop
    hl->TreeH()->GetEntry(iEnt);                                 //get current entry (prim)                
    for(Int_t iHit=0;iHit<h->Hits()->GetEntries();iHit++){       //hits loop
      AliHMPIDHit *pHit = (AliHMPIDHit*)h->Hits()->At(iHit);     //get current hit        
      new((*pHitLst)[cnt++]) AliHMPIDHit(*pHit);
    }//hits loop for this entry
  }//tree entries loop
  
}//ReadHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void sed()
{

  static TCanvas *pC1=0;
  
  if(!pC1){
    pC1=new TCanvas("hed","Simulated evets-View from electronics side, IP is behind the picture.",1000,900); pC1->Divide(3,3);
    pC1->cd(7); TButton *pBtn=new TButton("Next","sed()",0,0,0.2,0.1);   pBtn->Draw(); 
    pC1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, "","zoom(Int_t,Int_t,Int_t,TObject*)");
  }

  TClonesArray lh("AliHMPIDHit"); 
  TClonesArray ls("AliHMPIDDigit"); 
  TObjArray    ld(7); for(Int_t i=0;i<7;i++) ld.AddAt(new TClonesArray("AliHMPIDDigit"),i);
  TObjArray    lc(7); for(Int_t i=0;i<7;i++) lc.AddAt(new TClonesArray("AliHMPIDCluster"),i);
  AliESD esd;
  
  
  SimEsd(&esd);
  HitsFromEsd(&esd,&lh);
             AliHMPIDv1::Hit2Sdi(&lh,&ls);                               
      AliHMPIDDigitizer::Sdi2Dig(&ls,&ld);     
  AliHMPIDReconstructor::Dig2Clu(&ld,&lc);
//        AliHMPIDTracker::Recon(&esd,&cl);
  
  DrawEvt(pC1,&lh,&ld,&lc,&esd);  
}//SimEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void zoom(Int_t evt, Int_t pixx, Int_t pixy, TObject *obj)
{
  Printf("pad %i",gPad->GetNumber());
  static Int_t lvl=32; //zoom level
  if(evt!=5 && evt!=6) return;
  if(evt==5&&lvl==2)  return; //max zoom in
  if(evt==6&&lvl==32) return; //max zoom out
  
  Float_t x=gPad->AbsPixeltoX(pixx); Float_t y=gPad->AbsPixeltoY(pixy); 
             
  if(evt==5){ lvl=lvl/2; gPad->Range(x-lvl*2,y-lvl*2,x+lvl*2,y+lvl*2);} //zoom in
  else      { lvl=32;    gPad->Range(-10,-10,150,140); } //zoom out 
  ((TCanvas *)gTQSender)->SetTitle(Form("zoom %i",lvl));
  gPad->Modified();
  gPad->Update();                                              
}
void hed()
{//event display from files
  static TCanvas *pC=0;
  static Int_t iEvt=0;
  static Int_t iEvtTot=999;
  static TFile *pEsdFl=0;
  static TTree *pEsdTr=0;
  static AliESD *pEsd=0;
  
  

  
  if(!pC&&iEvt<iEvtTot){
    iEvt=0;
    iEvtTot=999;
    if(hl==0) {Printf("hed: no HMPID loader");return;}
    Printf("Opening session");
    pEsdFl=TFile::Open("AliESDs.root");     if(!pEsdFl || !pEsdFl->IsOpen()) return;//open AliESDs.root
    pEsdTr=(TTree*) pEsdFl->Get("esdTree"); if(!pEsdTr)                      return;//get ESD tree
    pEsdTr->SetBranchAddress("ESD", &pEsd);
    hl->LoadHits(); hl->LoadDigits(); hl->LoadRecPoints();
    iEvtTot=pEsdTr->GetEntries();
    pC=new TCanvas("hed","View from electronics side, IP is behind the picture.",1000,900);  pC->ToggleEventStatus(); pC->Divide(3,3);
    pC->cd(7); TButton *pBtn=new TButton("Next","hed()",0,0,0.2,0.1);   pBtn->Draw();
    pC->cd(7); TButton *pBtn=new TButton("Quit","Close_hed()",0.2,0,0.4,0.1);   pBtn->Draw(); 
    new TGedEditor(pC);
  }
 
  TLatex txt; txt.SetTextSize(0.1);
  TClonesArray hits("AliHMPIDHit");
      
  if(iEvt<iEvtTot){
    pEsdTr->GetEntry(iEvt); al->GetEvent(iEvt); 
    hl->TreeD()->GetEntry(0); hl->TreeR()->GetEntry(0);
    ReadHits(&hits); 
     
    pC->cd(3);  gPad->Clear(); txt.DrawLatex(0.2,0.2,Form("Event %i (total %i)",iEvt,iEvtTot));
    DrawEvt(pC,&hits,h->DigLst(),h->CluLst(),pEsd);
    
    iEvt++;
  }else{
    Printf("--- No more events available...Bye.");
    pC->Close();
    pC=0x0;
    iEvt=0;
    iEvtTot=999;
  }
}//hed()













void DrawEvt(TCanvas *pC,TClonesArray *pHitLst,TObjArray *pDigLst,TObjArray *pCluLst,AliESD *pEsd)
{//draws all the objects of current event in given canvas

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
    for(Int_t iHit=0;iHit<pHitLst->GetEntries();iHit++){
      AliHMPIDHit *pHit=(AliHMPIDHit*)pHitLst->At(iHit);
      if(pHit->Ch()==iCh)        pHit->Draw();  //draw hits
    }
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








































Reve()
{
  gSystem->Load("libMinuit.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libCDB.so");

  gSystem->Load("libGed.so");
  gSystem->Load("libRGL.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");

  gSystem->Load("libReve.so");

  gSystem->Load("libHMPIDbase.so");
  gSystem->Load("libHMPIDsim.so");
  gSystem->Load("libHMPIDrec.so");
  
  TGeoManager::Import("geometry.root");
  
  AliHMPIDParam *pParam=AliHMPIDParam::Instance();
  
  AliRunLoader *pAL=AliRunLoader::Open(); pAL->LoadHits  ("HMPID"); TTree *pHitT=pAL->GetTreeH("HMPID", false);
                                          pAL->LoadDigits("HMPID"); TTree *pDigT=pAL->GetTreeD("HMPID", false); 
                                          
  Reve::PointSet *pHitPnt=new Reve::PointSet("Hits"));
  TPolyMarker3D  *pDigPnt=new TPolyMarker3D;  pDigPnt->SetName("Digits"); pDigPnt->SetMarkerColor(kGreen);iPntCnt=0;
                                          
  TPointSelector ps(pHitT,pHitPnt,"fX:fY:fZ",""); ps.Select();
  
  TClonesArray *pDigLst=new TClonesArray("AliHMPIDDigit"); //this is tmp dig list per chamber
 
  for(Int_t iCh=0;iCh<7;iCh++){
    pDigT->SetBranchAddress(Form("HMPID%i",iCh+1),&pDigLst);
    pDigT->GetEntry(0);
    for(Int_t iDig=0;iDig<pDigLst->GetEntries();iDig++){
      AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigLst->At(iDig);    
      TVector2 lors=pParam->Pad2Loc(pDig->PadX(),pDig->PadY());
      TVector3 mars=pParam->Lors2Mars(iCh,lors.X(),lors.Y());
      pDigPnt->SetPoint(iPntCnt++,mars.X(),mars.Y(),mars.Z());
    }//digits loop for chamber        
  }//chambers loop
  
  if(!gReve) new Reve::RGTopFrame(0,0,0,2);
  gReve->AddGlobalRenderElement(new Reve::RenderElementObjPtr(pDigPnt));
  gReve->AddRenderElement(pHitPnt);
  gReve->Redraw3D();
}



void zzz(Int_t evt,Int_t px,Int_t py,TObject*o)
{
  Printf("evt %i (%i,%i) class %s",evt,px,py,o->IsA()->ClassName());
}

void Gui()
{
  TGMainFrame   *pMF     =new TGMainFrame(gClient->GetRoot(),300,400);//main frame
//1 level widgets: button and 2 horizontal frames  
  TRootEmbeddedCanvas *pDis;
  pMF->AddFrame(pDis=new TRootEmbeddedCanvas("display",pMF,800,650));
  pDis->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, "","zzz(Int_t,Int_t,Int_t,TObject*)");
  
  pMF->Layout();
  pMF->MapSubwindows();
  pMF->Resize(pMF->GetDefaultSize());
  pMF->MapWindow();
}                                                                      
