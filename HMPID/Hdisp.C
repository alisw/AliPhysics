TCanvas *pAll=0;
AliRunLoader *gAL=0; AliLoader *gHL=0; AliESD *gEsd=0; TTree *gEsdTr=0; AliHMPID *gH=0;
Int_t gEvt=0; Int_t gMaxEvt=0;
TObjArray *pNmean;

TChain *pCosCh=new TChain("cosmic");                                                               //clm: Define TChain for cosmic
TObjArray *pCosDigAll=0;                                                                           //clm: Define global Digits
TObjArray *pCosCluAll=0;                                                                           //clm: Define global Clusters
Int_t gCosRun=0;                                                                                   //clm: global cosmic event number
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hdisp(Int_t cosRun=44)                                                    //clm: Select cosmic file for display
{//display events from files if any in current directory or simulated events
  pAll=new TCanvas("all","",1300,900); pAll->Divide(3,3,0,0);
//  pAll->ToggleEditor();
  pAll->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,"","DoZoom(Int_t,Int_t,Int_t,TObject*)");

  OpenCalib();  
     if(gSystem->IsFileInIncludePath("galice.root")){// tries to open session
      if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
      gAL=AliRunLoader::Open();                                                                    //try to open galice.root from current dir 
      gAL->LoadgAlice();                                                                           //take new AliRun object from galice.root   
      gHL=gAL->GetDetectorLoader("HMPID");  gH=(AliHMPID*)gAL->GetAliRun()->GetDetector("HMPID");  //get HMPID object from galice.root
      gMaxEvt=gAL->GetNumberOfEvents()-1;
      gHL->LoadHits(); gHL->LoadSDigits(); gHL->LoadDigits(); gHL->LoadRecPoints();

      TFile *pEsdFl=TFile::Open("AliESDs.root"); gEsdTr=(TTree*) pEsdFl->Get("esdTree"); gEsdTr->SetBranchAddress("ESD", &gEsd);
      pAll->cd(7); TButton *pBtn=new TButton("Next","ReadEvt()",0,0,0.2,0.1);   pBtn->Draw();
                   TButton *pHitBtn=new TButton("Print hits","PrintHits()",0,0.2,0.3,0.3);   pHitBtn->Draw();
                   TButton *pSdiBtn=new TButton("Print sdis","PrintSdis()",0,0.4,0.3,0.5);   pSdiBtn->Draw();
                   TButton *pDigBtn=new TButton("Print digs","PrintDigs()",0,0.6,0.3,0.7);   pDigBtn->Draw();
                   TButton *pCluBtn=new TButton("Print clus","PrintClus()",0,0.8,0.3,0.9);   pCluBtn->Draw();  
      ReadEvt();
    }
    else if ( gSystem->IsFileInIncludePath(Form("cosmic%d.root",cosRun))){                          //clm: Check if cosmic file is in the folder
    gCosRun=cosRun;
    pCosCh->Add(Form("cosmic%d.root",gCosRun));                                                      //clm: Add cosmic file to chain
    pCosCh->SetBranchAddress("Digs",&pCosDigAll);                                                   //clm: Set digit branch address
    pCosCh->SetBranchAddress("Clus",&pCosCluAll);                                                   //clm: Set cluster branch address    
    gMaxEvt=pCosCh->GetEntries()-1;                                                                 //clm: Get number of events from the cosmic chain
    pAll->cd(7); TButton *pCosBtn=new TButton("Next Cosmic","ReadCosEvt()",0,0,0.3,0.1);   pCosBtn->Draw();   //clm: define next button
    ReadCosEvt();                                                                                   //clm: Read first cosmic event  
    }            
    else{
      pAll->cd(7); TButton *pBtn=new TButton("Next","SimEvt()",0,0,0.2,0.1);   pBtn->Draw(); 
      SimEvt();
          }      
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ReadCosEvt()
{// Read curent cosmic event and display it assumes that session is alredy opened
  if(gEvt>gMaxEvt) gEvt=0; if(gEvt<0) gEvt=gMaxEvt;                                     //clm: set event limits
  pCosCh->GetEntry(gEvt);                                                               //clm: read event from chain
  DrawCosEvt(pCosDigAll,pCosCluAll);                                                    //clm: draw cosmic event
  gEvt++;                                                                             
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ReadEvt()
{// Read curent event and display it assumes that session is alredy opened
  TClonesArray hits("AliHMPIDHit");
  if(gEvt>gMaxEvt) gEvt=0; if(gEvt<0) gEvt=gMaxEvt;
  
  gEsdTr->GetEntry(gEvt); gAL->GetEvent(gEvt); 
  ReadHits(&hits); gHL->TreeS()->GetEntry(0); gHL->TreeD()->GetEntry(0); gHL->TreeR()->GetEntry(0);
    
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
  AliHMPIDDigit::fSigmas=4;
  AliHMPIDDigitizer::DoNoise(kFALSE);
  SimEsd(&esd);
  SimHits(&esd,&hits);
             AliHMPIDv1::Hit2Sdi(&hits,&sdig);                               
      AliHMPIDDigitizer::Sdi2Dig(&sdig,&digs);     
      AliHMPIDReconstructor::Dig2Clu(&digs,&clus);
      AliHMPIDTracker::Recon(&esd,&clus,pNmean);
  
  pAll->cd(3);  gPad->Clear(); TLatex txt;txt.DrawLatex(0.2,0.2,Form("Simulated event %i",gEvt));
  DrawEvt(&hits,&digs,&clus,&esd);  
  gEvt++;
}//SimEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsd(AliESD *pEsd)
{
  TParticle part; TLorentzVector mom;
  for(Int_t iTrk=0;iTrk<100;iTrk++){//stack loop
    part.SetPdgCode(kProton);
    part.SetProductionVertex(0,0,0,0);
    Double_t eta= -0.2+gRandom->Rndm()*0.4; //rapidity is random [-0.2,+0.2]
    Double_t phi= gRandom->Rndm()*60.*TMath::DegToRad();   //phi is random      [ 0  , 60 ] degrees    
    mom.SetPtEtaPhiM(5,eta,phi,part.GetMass());
    part.SetMomentum(mom);
    AliESDtrack trk(&part);
    pEsd->AddTrack(&trk);
  }//stack loop  
}//EsdFromStack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimHits(AliESD *pEsd, TClonesArray *pHits)
{
  const Int_t kCerenkov=50000050;
  const Int_t kFeedback=50000051;
  
  AliHMPIDRecon rec;
  
  Int_t hc=0; 
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);
    Float_t xRa,yRa;
    Int_t ch=AliHMPIDTracker::IntTrkCha(pTrk,xRa,yRa);
    if(ch<0) continue; //this track does not hit HMPID
    Float_t beta = pTrk->GetP()/(TMath::Sqrt(pTrk->GetP()*pTrk->GetP()+0.938*0.938));
    Float_t ckov=TMath::ACos(1./(beta*1.292));

    Float_t theta,phi,xPc,yPc,; 
    pTrk->GetHMPIDtrk(xPc,yPc,theta,phi); 
    rec.SetTrack(xRa,yRa,theta,phi); 
    
    if(!AliHMPIDDigit::IsInDead(xPc,yPc)) new((*pHits)[hc++]) AliHMPIDHit(ch,200e-9,kProton  ,iTrk,xPc,yPc);                 //mip hit
    Int_t nPhots = (Int_t)(20.*TMath::Power(TMath::Sin(ckov),2)/TMath::Power(TMath::Sin(TMath::ACos(1./1.292)),2));
    for(int i=0;i<nPhots;i++){
      TVector2 pos;
      pos=rec.TracePhot(ckov,gRandom->Rndm()*TMath::TwoPi());
      if(!AliHMPIDDigit::IsInDead(pos.X(),pos.Y())) new((*pHits)[hc++]) AliHMPIDHit(ch,7.5e-9,kCerenkov,iTrk,pos.X(),pos.Y());
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
  if(iCh>=0){TLatex txt; txt.SetTextSize(0.1); txt.DrawLatex(-5,-5,Form("%i",iCh));}
  
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
  TPolyMarker *pTxC[7], *pRin[7]; TMarker *pMip,*pCko,*pFee,*pDig,*pClu;
  pMip=new TMarker; pMip->SetMarkerColor(kRed);  pMip->SetMarkerStyle(kOpenTriangleUp);
  pCko=new TMarker; pCko->SetMarkerColor(kRed);  pCko->SetMarkerStyle(kOpenCircle);
  pFee=new TMarker; pFee->SetMarkerColor(kRed);  pFee->SetMarkerStyle(kOpenDiamond);
  pDig=new TMarker; pDig->SetMarkerColor(kGreen);pDig->SetMarkerStyle(kOpenSquare);
  pClu=new TMarker; pClu->SetMarkerColor(kBlue); pClu->SetMarkerStyle(kStar);
  for(Int_t ch=0;ch<7;ch++){
    pTxC[ch]=new TPolyMarker; pTxC[ch]->SetMarkerStyle(2); pTxC[ch]->SetMarkerColor(kRed); pTxC[ch]->SetMarkerSize(3);
    pRin[ch]=new TPolyMarker; pRin[ch]->SetMarkerStyle(6); pRin[ch]->SetMarkerColor(kMagenta);
  }
  
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop to collect cerenkov rings and intersection points
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);
    Int_t ch=pTrk->GetHMPIDcluIdx();
    if(ch<0) continue; //this track does not hit HMPID
    ch/=1000000; 
    Float_t th,ph,xRad,yRad; pTrk->GetHMPIDtrk(xRad,yRad,th,ph);//get info on current track
//    pTxC[ch]->SetNextPoint(xPc,yPc);                          //add new intersection point  TEMPORARLY DISABLED...no more available in ESD!
    
    Float_t ckov=pTrk->GetHMPIDsignal();  Float_t err=TMath::Sqrt(pTrk->GetHMPIDchi2());
    if(ckov>0){
      Printf("theta %f phi %f ckov %f",th*TMath::RadToDeg(),ph*TMath::RadToDeg(),ckov);
      rec.SetTrack(xRad,yRad,th,ph);
      for(Int_t j=0;j<100;j++){ 
        TVector2 pos;
        pos=rec.TracePhot(ckov,j*0.0628);
       if(!AliHMPIDDigit::IsInDead(pos.X(),pos.Y())) pRin[ch]->SetNextPoint(pos.X(),pos.Y());
      }      
    }
  }//tracks loop
      
  Int_t totHit=pHitLst->GetEntriesFast(),totDig=0,totClu=0,totTxC=0;
  Int_t totCkov=0, totFeed=0,totMip=0;
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop    
    totTxC+=pTxC[iCh]->GetN();
    totDig+=((TClonesArray*)pDigLst->At(iCh))->GetEntriesFast();
    totClu+=((TClonesArray*)pCluLst->At(iCh))->GetEntriesFast();
    
    switch(iCh){
      case 6: pAll->cd(1); break; case 5: pAll->cd(2); break;
      case 4: pAll->cd(4); break; case 3: pAll->cd(5); break; case 2: pAll->cd(6); break;
                                  case 1: pAll->cd(8); break; case 0: pAll->cd(9); break;
    }
    gPad->SetEditable(kTRUE); gPad->Clear(); 
    DrawCh(iCh);
                           
    ((TClonesArray*)pDigLst->At(iCh))->Draw();               //draw digits
    
    totCkov=0;totFeed=0;totMip=0;
    for(Int_t iHit=0;iHit<pHitLst->GetEntriesFast();iHit++){ // Draw hits
      AliHMPIDHit *pHit=(AliHMPIDHit*)pHitLst->At(iHit);
      if(pHit->Ch()==iCh) pHit->Draw();
      if(pHit->Pid()==50000050) totCkov++;
      else if(pHit->Pid()==50000051) totFeed++;
      else totMip++;
    }    
           
    ((TClonesArray*)pCluLst->At(iCh))->Draw();              //draw clusters
                            pTxC[iCh]->Draw();              //draw intersections
                            pRin[iCh]->Draw("CLP");         //draw rings
    gPad->SetEditable(kFALSE);
  }//chambers loop
  
  pAll->cd(3);  gPad->Clear(); TLegend *pLeg=new TLegend(0.2,0.2,0.8,0.8);
                                        pLeg->SetHeader(Form("Event %i Total %i",gEvt,gMaxEvt+1));
                                        pLeg->AddEntry(pTxC[0],Form("TRKxPC %i"     ,totTxC),"p");
                                        pLeg->AddEntry(pMip   ,Form("Mip hits %i"   ,totMip),"p");    
                                        pLeg->AddEntry(pCko   ,Form("Ckov hits %i"  ,totCkov),"p");    
                                        pLeg->AddEntry(pFee   ,Form("Feed hits %i"  ,totFeed),"p");    
                                        pLeg->AddEntry(pDig   ,Form("Digs %i"       ,totDig),"p");    
                                        pLeg->AddEntry(pClu   ,Form("Clus %i"       ,totClu),"p");    
                                        pLeg->Draw();
}//DrawEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawCosEvt(TObjArray *pCosDigLst,TObjArray *pCosCluLst)                              //clm: new method to read and display cosmic events 
{//draws all the objects of current event in given canvas
  TMarker *pDig,*pClu;
  pDig=new TMarker; pDig->SetMarkerColor(kGreen);pDig->SetMarkerStyle(kOpenSquare);
  pClu=new TMarker; pClu->SetMarkerColor(kBlue); pClu->SetMarkerStyle(kStar);
  
      
  Int_t totDig=0,totClu=0;
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop                                            //clm: chambers loop is now not needed since iCh=0 but kept for provision
    totDig+=((TClonesArray*)pCosDigLst->At(iCh))->GetEntriesFast();
    totClu+=((TClonesArray*)pCosCluLst->At(iCh))->GetEntriesFast();
    
    switch(iCh){
      case 6: pAll->cd(1); break; case 5: pAll->cd(2); break;
      case 4: pAll->cd(4); break; case 3: pAll->cd(5); break; case 2: pAll->cd(6); break;
                                  case 1: pAll->cd(8); break; case 0: pAll->cd(9); break;
    }
   gPad->SetEditable(kTRUE); gPad->Clear(); 
   DrawCh(iCh);
  if(gCosRun < 500 && iCh == 6) {                           //clm: hard coded selection since raw data does not contain ch info which is set to 0
  ((TClonesArray*)pCosDigLst->At(0))->Draw();               //draw digits
  ((TClonesArray*)pCosCluLst->At(0))->Draw();              //draw clusters
  }
  if (gCosRun >= 500 && iCh == 5) { 
  ((TClonesArray*)pCosDigLst->At(0))->Draw();               //draw digits
  ((TClonesArray*)pCosCluLst->At(0))->Draw();              //draw clusters
  }
   
   gPad->SetEditable(kFALSE);
  }//chambers loop
  pAll->cd(3);  gPad->Clear(); TLegend *pLeg=new TLegend(0.2,0.2,0.8,0.8);
                                        pLeg->SetHeader(Form("Cosmic Run %i Event %i Total %i",gCosRun,gEvt,gMaxEvt+1));
                                        pLeg->AddEntry(pDig   ,Form("Digs %i"       ,totDig),"p");    
                                        pLeg->AddEntry(pClu   ,Form("Clus %i"       ,totClu),"p");    
                                        pLeg->Draw();
}//DrawCosEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DoZoom(Int_t evt, Int_t px, Int_t py, TObject *obj)
{
  if(evt!=5 && evt!=6) return; //5- zoom in 6-zoom out
  const Int_t minZoom=64;
  const Int_t maxZoom=2;
  static Int_t zoom=minZoom; //zoom level
  if(evt==5&&zoom==maxZoom) return; 
  if(evt==6&&zoom==minZoom) return; 
  
 // if(!obj->IsA()->InheritsFrom("TPad")) return;  //current object is not pad
  TVirtualPad *pPad=gPad->GetSelectedPad();
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
void PrintHits()
{
//Prints a list of HMPID hits for a given event. Default is event number 0.
  Printf("List of HMPID hits for event %i",gEvt);
  
  Int_t iTotHits=0;
  for(Int_t iPrim=0;iPrim<gHL->TreeH()->GetEntries();iPrim++){//prims loop
    gHL->TreeH()->GetEntry(iPrim);      
    gH->Hits()->Print();
    iTotHits+=gH->Hits()->GetEntries();
  }
  Printf("totally %i hits for event %i",iTotHits,gEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintSdis()
{
//prints a list of HMPID sdigits  for a given event
  Printf("List of HMPID sdigits for event %i",gEvt);
  gH->SdiLst()->Print();
  Printf("totally %i sdigits for event %i",gH->SdiLst()->GetEntries(),gEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintDigs()
{
//prints a list of HMPID digits
  Printf("List of HMPID digits for event %i",gEvt);  
  gH->DigLst()->Print();
  Int_t totDigs=0;  for(Int_t i=0;i<7;i++) {totDigs+=gH->DigLst(i)->GetEntries();}
  Printf("totally %i digits for event %i",totDigs,gEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintClus()
{//prints a list of HMPID clusters  for a given event
  Printf("List of HMPID clusters for event %i",gEvt);
  gH->CluLst()->Print();
  
  Int_t iCluCnt=0; for(Int_t iCh=0;iCh<7;iCh++) iCluCnt+=gH->CluLst(iCh)->GetEntries();
  
  Printf("totally %i clusters for event %i",iCluCnt,gEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OpenCalib()
{
  AliCDBManager* pCDB = AliCDBManager::Instance();
  pCDB->SetDefaultStorage("local://$HOME");
  AliCDBEntry *pQthreEnt=pCDB->Get("HMPID/Calib/Qthre",0);
  AliCDBEntry *pNmeanEnt=pCDB->Get("HMPID/Calib/Nmean",0);
  
  if(!pQthreEnt || ! pNmeanEnt) return;
  
  pNmean=(TObjArray*)pNmeanEnt->GetObject(); 
}
