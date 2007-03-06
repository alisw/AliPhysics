TCanvas *pAll=0;
AliRunLoader *gAL=0; AliLoader *gHL=0; AliESD *gEsd=0; TTree *gEsdTr=0; AliHMPID *gH=0;
Int_t gEvt=0; Int_t gMaxEvt=0;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hdisp()
{//display events from files if any in current directory or simulated events
  pAll=new TCanvas("all","",1000,900); pAll->Divide(3,3,0,0);pAll->ToggleEditor();
  pAll->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,"","DoZoom(Int_t,Int_t,Int_t,TObject*)");
  
  if(gSystem->IsFileInIncludePath("galice.root")){// tries to open session
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    gAL=AliRunLoader::Open();                                                                    //try to open galice.root from current dir 
    gAL->LoadgAlice();                                                                           //take new AliRun object from galice.root   
    gHL=gAL->GetDetectorLoader("HMPID");  gH=(AliHMPID*)gAL->GetAliRun()->GetDetector("HMPID");  //get HMPID object from galice.root
    gMaxEvt=gAL->GetNumberOfEvents()-1;
    gHL->LoadHits(); gHL->LoadDigits(); gHL->LoadRecPoints();

    TFile *pEsdFl=TFile::Open("AliESDs.root"); gEsdTr=(TTree*) pEsdFl->Get("esdTree"); gEsdTr->SetBranchAddress("ESD", &gEsd);
    pAll->cd(7); TButton *pBtn=new TButton("Next","ReadEvt()",0,0,0.2,0.1);   pBtn->Draw();
    ReadEvt();
  }else{
    pAll->cd(7); TButton *pBtn=new TButton("Next","SimEvt()",0,0,0.2,0.1);   pBtn->Draw(); 
    SimEvt();
  }      
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ReadEvt()
{// Read curent event and display it assumes that session is alredy opened
  TClonesArray hits("AliHMPIDHit");
  if(gEvt>gMaxEvt) gEvt=0; if(gEvt<0) gEvt=gMaxEvt;
  
  gEsdTr->GetEntry(gEvt); gAL->GetEvent(gEvt); 
  ReadHits(&hits); gHL->TreeD()->GetEntry(0); gHL->TreeR()->GetEntry(0);
  
  pAll->cd(3);  gPad->Clear(); TLatex txt;txt.DrawLatex(0.2,0.2,Form("Event %i (total %i)",gEvt,gMaxEvt));
  DrawEvt(&hits,gH->DigLst(),gH->CluLst(),gEsd);
  gEvt++;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEvt()
{
  TClonesArray hits("AliHMPIDHit"); 
  TClonesArray sdig("AliHMPIDDigit"); 
  TObjArray    digs(7); for(Int_t i=0;i<7;i++) digs.AddAt(new TClonesArray("AliHMPIDDigit"),i);
  TObjArray    clus(7); for(Int_t i=0;i<7;i++) clus.AddAt(new TClonesArray("AliHMPIDCluster"),i);
  AliESD esd;
    
  SimEsd(&esd);
  SimHits(&esd,&hits);
             AliHMPIDv1::Hit2Sdi(&hits,&sdig);                               
      AliHMPIDDigitizer::Sdi2Dig(&sdig,&digs);     
  AliHMPIDReconstructor::Dig2Clu(&digs,&clus);
        AliHMPIDTracker::Recon(&esd,&clus);
  
  pAll->cd(3);  gPad->Clear(); TLatex txt;txt.DrawLatex(0.2,0.2,Form("Simulated event %i",gEvt));
  DrawEvt(&hits,&digs,&clus,&esd);  
  gEvt++;
}//SimEvt()
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
void SimHits(AliESD *pEsd, TClonesArray *pHits)
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
    
    if(!AliHMPIDDigit::IsInDead(xPc,yPc)) new((*pHits)[hc++]) AliHMPIDHit(ch,200e-9,kProton  ,iTrk,xPc,yPc);                 //mip hit
    for(int i=0;i<4;i++)  new((*pHits)[hc++]) AliHMPIDHit(ch,7.5e-9,kFeedback,iTrk,gRandom->Rndm()*130,gRandom->Rndm()*126); //bkg hits 4 per track
    for(int i=0;i<16;i++){
      rec.TracePhot(ckov,gRandom->Rndm()*TMath::TwoPi(),pos);
      new((*pHits)[hc++]) AliHMPIDHit(ch,7.5e-9,kCerenkov,iTrk,pos.X(),pos.Y());
    }                      //photon hits  
  }//tracks loop    
}//SimHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ReadHits(TClonesArray *pHitLst)
{
  pHitLst->Delete();  Int_t cnt=0;
  for(Int_t iEnt=0;iEnt<gHL->TreeH()->GetEntries();iEnt++){       //TreeH loop
    gHL->TreeH()->GetEntry(iEnt);                                 //get current entry (prim)                
    for(Int_t iHit=0;iHit<gH->Hits()->GetEntries();iHit++){       //hits loop
      AliHMPIDHit *pHit = (AliHMPIDHit*)gH->Hits()->At(iHit);     //get current hit        
      new((*pHitLst)[cnt++]) AliHMPIDHit(*pHit);
    }//hits loop for this entry
  }//tree entries loop
  
}//ReadHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawCh(Int_t iCh) 
{ 
  gPad->Range(-10,-10,AliHMPIDDigit::SizeAllX()+5,AliHMPIDDigit::SizeAllY()+5); 
  TLatex txt; txt.SetTextSize(0.1); txt.DrawLatex(-5,-5,Form("%i",iCh));
  
  for(Int_t iPc=AliHMPIDDigit::kMinPc;iPc<=AliHMPIDDigit::kMaxPc;iPc++){
    TBox *pBox=new TBox(AliHMPIDDigit::fMinPcX[iPc],AliHMPIDDigit::fMinPcY[iPc],
                        AliHMPIDDigit::fMaxPcX[iPc],AliHMPIDDigit::fMaxPcY[iPc]);
    
    pBox->SetFillStyle(0);
    pBox->Draw();
    
  }//PC loop      
}//DrawCh()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawEvt(TClonesArray *pHitLst,TObjArray *pDigLst,TObjArray *pCluLst,AliESD *pEsd)
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
      case 6: pAll->cd(1); break; case 5: pAll->cd(2); break;
      case 4: pAll->cd(4); break; case 3: pAll->cd(5); break; case 2: pAll->cd(6); break;
                                  case 1: pAll->cd(8); break; case 0: pAll->cd(9); break;
    }
    gPad->SetEditable(kTRUE); gPad->Clear(); 
    DrawCh(iCh);
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
}//DrawEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DoZoom(Int_t evt, Int_t px, Int_t py, TObject *obj)
{
  if(evt!=5 && evt!=6) return; //5- zoom in 6-zoom out
  const Int_t minZoom=64;
  const Int_t maxZoom=2;
  static Int_t zoom=minZoom; //zoom level
  if(evt==5&&zoom==maxZoom) return; 
  if(evt==6&&zoom==minZoom) return; 
  
  if(!obj->IsA()->InheritsFrom("TPad")) return;  //current object is not pad
  TPad *pPad=(TPad*)obj;
  if(pPad->GetNumber()==3 || pPad->GetNumber()==7) return; //current pad is wrong

 // Printf("evt=%i (%i,%i) %s",evt,px,py,obj->GetName());
    
  Float_t x=pPad->AbsPixeltoX(px); Float_t y=pPad->AbsPixeltoY(py); 
 
  if(evt==5){ zoom=zoom/2;     pPad->Range(x-zoom*2,y-zoom*2,x+zoom*2,y+zoom*2);} //zoom in
  else      { zoom=zoom*2;     pPad->Range(x-zoom*2,y-zoom*2,x+zoom*2,y+zoom*2);} //zoom out 
  if(zoom==minZoom) pPad->Range(-10,-10,AliHMPIDDigit::SizeAllX()+5,AliHMPIDDigit::SizeAllY()+5);
  ((TCanvas *)gTQSender)->SetTitle(Form("zoom x%i",minZoom/zoom));
  pPad->Modified();
  pPad->Update();                                              
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
