#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TSystem.h>
#include <TFile.h>
#include <TVirtualX.h>
#include <TTree.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TCanvasImp.h>
#include <TStyle.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TPDGCode.h>
#include <TLegend.h>
#include <TPolyMarker.h>
#include <TBox.h>
#include <AliESDEvent.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include "AliHMPIDHit.h"
#include "AliHMPIDv2.h"
#include "AliHMPIDReconstructor.h"
#include "AliHMPIDRecon.h"
#include "AliHMPIDParam.h"
#include "AliHMPIDCluster.h"
#include "AliTracker.h"

#endif

AliHMPIDParam *fParam;
TDatabasePDG *fPdg;

TCanvas *fCanvas=0; Int_t fType=3; Int_t fEvt=-1; Int_t fNevt=0;                      
TCanvasImp *fCanvasImp;
TFile *fHitFile; TTree *fHitTree; TClonesArray *fHitLst; TPolyMarker *fRenMip[7]; TPolyMarker *fRenCko[7]; TPolyMarker *fRenFee[7];
                                  TClonesArray *fSdiLst; 
TFile *fDigFile; TTree *fDigTree; TObjArray    *fDigLst; TBox *fRenDig[7][160*144]; TBox *fBox[7][160*144];   
TFile *fCluFile; TTree *fCluTree; TObjArray    *fCluLst; TPolyMarker *fRenClu[7];
TFile *fEsdFile; TTree *fEsdTree; AliESDEvent  *fEsd;    TPolyMarker *fRenTxC[7]; TPolyMarker *fRenRin[7];  
TFile *fCosFile; TTree *fCosTree;

TButton *fHitMipBok,*fHitCkoBok,*fHitFeeBok,*fDigBok,*fCluBok,*fEsdBok;

TString fStHitMip = "ON";
TString fStHitCko = "ON";
TString fStHitFee = "ON";
TString fStDig    = "ON";
TString fStClu    = "ON";
TString fStEsd    = "ON";

Int_t fTimeArrow = 1;

Int_t fMaxCharge = 200;  //maximum charge to saturate color to red

Int_t fChamN[9]={6,5,-1,4,3,2,-1,1,0};

Int_t fTotPads[7],fTotClus[7];

AliRunLoader *gAL=0; 

Int_t nDigs[7];

enum EObjectType {kHitMip=1,kHitCko,kHitFee,kDigit,kCluster,kTrack,kRing};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CreateContainers()
{//to create all containers
  fHitLst=new TClonesArray("AliHMPIDHit");
  fSdiLst=new TClonesArray("AliHMPIDDigit");
  fDigLst=new TObjArray(7); for(Int_t i=0;i<7;i++) fDigLst->AddAt(new TClonesArray("AliHMPIDDigit"),i);       fDigLst->SetOwner(kTRUE);
  fCluLst=new TObjArray(7); for(Int_t i=0;i<7;i++) fCluLst->AddAt(new TClonesArray("AliHMPIDCluster"),i);     fCluLst->SetOwner(kTRUE); 
  fEsd   =new AliESDEvent;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CreateRenders()
{
  for(Int_t ch=0;ch<7;ch++){
    fRenMip[ch]=new TPolyMarker; fRenMip[ch]->SetMarkerStyle(kOpenTriangleUp);fRenMip[ch]->SetMarkerColor(kBlack)  ;fRenMip[ch]->SetUniqueID(kHitMip);
    fRenCko[ch]=new TPolyMarker; fRenCko[ch]->SetMarkerStyle(kOpenCircle);    fRenCko[ch]->SetMarkerColor(kBlack)  ;fRenCko[ch]->SetUniqueID(kHitCko);
    fRenFee[ch]=new TPolyMarker; fRenFee[ch]->SetMarkerStyle(kOpenDiamond);   fRenFee[ch]->SetMarkerColor(kBlack)  ;fRenFee[ch]->SetUniqueID(kHitFee);
    fRenClu[ch]=new TPolyMarker; fRenClu[ch]->SetMarkerStyle(kStar);          fRenClu[ch]->SetMarkerColor(kMagenta);fRenClu[ch]->SetUniqueID(kCluster);
    fRenTxC[ch]=new TPolyMarker; fRenTxC[ch]->SetMarkerStyle(kPlus);          fRenTxC[ch]->SetMarkerColor(kRed)    ;fRenTxC[ch]->SetUniqueID(kTrack);
                                 fRenTxC[ch]->SetMarkerSize(3);
    fRenRin[ch]=new TPolyMarker; fRenRin[ch]->SetMarkerStyle(kFullDotSmall);  fRenRin[ch]->SetMarkerColor(kMagenta);fRenRin[ch]->SetUniqueID(kRing);

    for(Int_t iDig=0;iDig<160*144;iDig++) {
      fRenDig[ch][iDig] = new TBox;
      fRenDig[ch][iDig]->SetFillStyle(1);
      fBox[ch][iDig] = new TBox;
      fBox[ch][iDig]->SetFillStyle(0);
    }
  }
}//CreateRenders()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ClearRenders()
{
  for(Int_t ch=0;ch<7;ch++){
    fRenTxC[ch]->SetPolyMarker(0); 
    fRenRin[ch]->SetPolyMarker(0); 
    fRenMip[ch]->SetPolyMarker(0); 
    fRenCko[ch]->SetPolyMarker(0); 
    fRenFee[ch]->SetPolyMarker(0); 
    nDigs[ch] = 0;
    fRenClu[ch]->SetPolyMarker(0); 
  }
}//ClearRenders()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintHits()
{
//Prints a list of HMPID hits for a given event. Default is event number 0.
  if(!fHitTree) return;
  Printf("List of HMPID hits for event %i",fEvt);
  Int_t iTot=0;
  for(Int_t iEnt=0;iEnt<fHitTree->GetEntries();iEnt++){//entries loop
    fHitTree->GetEntry(iEnt);      
    fHitLst->Print();
    iTot+=fHitLst->GetEntries();
  }
  Printf("totally %i hits for event %i",iTot,fEvt);
}//PrintHits();
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintSdis()
{//prints a list of HMPID sdigits  for a given event
  Printf("List of HMPID sdigits for event %i",fEvt);
  fSdiLst->Print();
  Printf("totally %i sdigits for event %i",fSdiLst->GetEntries(),fEvt);
}//PrintSdis()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintDigs()
{//prints a list of HMPID digits
  Printf("List of HMPID digits for event %i",fEvt);  
  fDigLst->Print();
  Int_t iTot=0;  for(Int_t i=0;i<7;i++) {iTot+=((TClonesArray*)fDigLst->At(i))->GetEntries();}
  Printf("totally %i digits for event %i",iTot,fEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintClus()
{//prints a list of HMPID clusters  for a given event
  Printf("List of HMPID clusters for event %i",fEvt);
  
  fCluLst->Print();
  Int_t iTot=0; for(Int_t iCh=0;iCh<7;iCh++) iTot+=((TClonesArray*)fCluLst->At(iCh))->GetEntries();
  
  Printf("totally %i clusters for event %i",iTot,fEvt);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrintEsd()
{//prints a list of HMPID Esd  for a given event
  Printf("List of HMPID ESD summary for event %i",fEvt);
  for(Int_t iTrk=0;iTrk<fEsd->GetNumberOfTracks();iTrk++){
    AliESDtrack *pTrk = fEsd->GetTrack(iTrk);
    Float_t x,y;Int_t q,nacc;   pTrk->GetHMPIDmip(x,y,q,nacc);
    Printf("Track %02i with p %7.2f with ThetaCer %5.3f with %i photons",iTrk,pTrk->GetP(),pTrk->GetHMPIDsignal(),nacc);
  }  
}//PrintEsd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawChamber(Int_t iCh) 
{//used by Draw() to Draw() chamber structure
  gPad->Range(-10,-10,AliHMPIDParam::SizeAllX()+5,AliHMPIDParam::SizeAllY()+5);
  if(iCh>=0){TLatex txt; txt.SetTextSize(0.1); txt.DrawLatex(-5,-5,Form("%i",iCh));}
  
  for(Int_t iPc=AliHMPIDParam::kMinPc;iPc<=AliHMPIDParam::kMaxPc;iPc++){
    TBox *pBox=new TBox(AliHMPIDParam::MinPcX(iPc),AliHMPIDParam::MinPcY(iPc),
                        AliHMPIDParam::MaxPcX(iPc),AliHMPIDParam::MaxPcY(iPc));
    pBox->SetFillStyle(0);  pBox->Draw();
  }//PC loop      
}//DrawChamber()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawLegend()
{//used by Draw() to draw legend
  Int_t nTxC=0,nMip=0,nCko=0,nFee=0,nDig=0,nClu=0;
  for(Int_t ch=0;ch<7;ch++){
    nTxC+=fRenTxC[ch]->Size();
    nMip+=fRenMip[ch]->Size();
    nCko+=fRenCko[ch]->Size();
    nFee+=fRenFee[ch]->Size();
    nClu+=fRenClu[ch]->Size();
    nDig+=nDigs[ch];
  }
  TLegend *pLeg=new TLegend(0.2,0.2,0.8,0.8);
  pLeg->SetHeader(Form("Event %i Total %i",fEvt,fNevt));
  pLeg->AddEntry(fRenTxC[0],   Form("TRKxPC %i"     ,nTxC),"p");
  pLeg->AddEntry(fRenMip[0],   Form("Mip hits %i"   ,nMip),"p");    
  pLeg->AddEntry(fRenCko[0],   Form("Ckov hits %i"  ,nCko),"p");    
  pLeg->AddEntry(fRenFee[0],   Form("Feed hits %i"  ,nFee),"p");    
  pLeg->AddEntry(fRenDig[0][0],Form("Digs %i"       ,nDig),"p");    
  pLeg->AddEntry(fRenClu[0],   Form("Clus %i"       ,nClu),"p");    
  pLeg->Draw();
}//DrawLegend()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Draw()
{//draws all the objects of current event in given canvas
  Int_t iPadN[7]={9,8,6,5,4,2,1};
  Int_t sampleCol = fMaxCharge/gStyle->GetNumberOfColors();
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop    
    fCanvas->cd(iPadN[iCh]);
    gPad->SetEditable(kTRUE); gPad->Clear(); DrawChamber(iCh);
    
    if(fDigFile){
      if(fStDig   =="ON") {
        for(Int_t iDig=0;iDig<nDigs[iCh];iDig++) {
          Int_t charge = fRenDig[iCh][iDig]->GetUniqueID();
          Int_t color = charge/sampleCol;if(color>=gStyle->GetNumberOfColors()) color = gStyle->GetNumberOfColors()-1;
          fRenDig[iCh][iDig]->SetFillColor(gStyle->GetColorPalette(color));
          fRenDig[iCh][iDig]->Draw();
          fBox[iCh][iDig]->Draw();
        }
      }
    }
    
    if(fEsdFile){if(fStEsd   =="ON") fRenTxC[iCh]->Draw();}
    
    if(fHitFile){
      if(fStHitMip=="ON") fRenMip[iCh]->Draw();
      if(fStHitFee=="ON") fRenFee[iCh]->Draw();
      if(fStHitCko=="ON") fRenCko[iCh]->Draw();
    }
    
    if(fEsdFile){if(fStEsd   =="ON") fRenRin[iCh]->Draw();}
    if(fCluFile){if(fStClu   =="ON") fRenClu[iCh]->Draw();}
    gPad->SetEditable(kFALSE);
  }//chambers loop
  
  fCanvas->cd(3);  gPad->Clear(); DrawLegend();
  
  fCanvas->cd(7);
  
  if(fHitFile){
    fHitMipBok->SetTitle(Form("Mips  %s",fStHitMip.Data())); fHitMipBok->Modified();
    fHitCkoBok->SetTitle(Form("Ckov  %s",fStHitCko.Data())); fHitCkoBok->Modified();
    fHitFeeBok->SetTitle(Form("Fdbk  %s",fStHitFee.Data())); fHitFeeBok->Modified();
  }
  if(fDigFile){fDigBok->SetTitle(fStDig); fDigBok->Modified();}  
  if(fCluFile){fCluBok->SetTitle(fStClu); fCluBok->Modified();} 
  if(fEsdFile){fEsdBok->SetTitle(fStEsd); fEsdBok->Modified();}
  
}//Draw()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RenderHit(TClonesArray *pHitLst)
{//used by ReadEvent() and SimulateEvent() to render hits to polymarker structures, one per chamber
  for(Int_t iHit=0;iHit<pHitLst->GetEntries();iHit++){       //hits loop
    AliHMPIDHit *pHit = (AliHMPIDHit*)pHitLst->At(iHit); Int_t ch=pHit->Ch(); Float_t x=pHit->LorsX(); Float_t y=pHit->LorsY();    //get current hit        
    switch(pHit->Pid()){
      case 50000050: fRenCko[ch]->SetNextPoint(x,y);break; 
      case 50000051: fRenFee[ch]->SetNextPoint(x,y);break;
      default:       fRenMip[ch]->SetNextPoint(x,y);break;
    }//switch hit PID      
//    Printf("----------ihit %i ch %i chhit %i MIP %i CKO %i FEE %i",iHit,ch,pHit->Pid());
  }//hits loop for this entry
}//RenderHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RenderDig(TObjArray *pDigLst)
{//used by ReadEvent() and SimulateEvent() to render digs to Tbox structures, one per chamber
  for(Int_t iCh=0;iCh<=6;iCh++){                                    //chambers loop   
    TClonesArray *pDigCham=(TClonesArray*)pDigLst->At(iCh);         //get digs list for this chamber
    nDigs[iCh] = 0;
    fTotPads[iCh] = pDigCham->GetEntries();
    for(Int_t iDig=0;iDig<pDigCham->GetEntries();iDig++){            //digs loop
      AliHMPIDDigit *pDig = (AliHMPIDDigit*)pDigCham->At(iDig); Float_t x=pDig->LorsX(); Float_t y=pDig->LorsY();    //get current hit        
      fRenDig[iCh][iDig]->SetX1(x-0.5*AliHMPIDParam::SizePadX());
      fRenDig[iCh][iDig]->SetX2(x+0.5*AliHMPIDParam::SizePadX());
      fRenDig[iCh][iDig]->SetY1(y-0.5*AliHMPIDParam::SizePadY());
      fRenDig[iCh][iDig]->SetY2(y+0.5*AliHMPIDParam::SizePadY());
      fRenDig[iCh][iDig]->SetUniqueID((Int_t)pDig->Q());
      fBox[iCh][iDig]->SetX1(x-0.5*AliHMPIDParam::SizePadX());
      fBox[iCh][iDig]->SetX2(x+0.5*AliHMPIDParam::SizePadX());
      fBox[iCh][iDig]->SetY1(y-0.5*AliHMPIDParam::SizePadY());
      fBox[iCh][iDig]->SetY2(y+0.5*AliHMPIDParam::SizePadY());
      nDigs[iCh]++;
    }//Digit loop
  }//hits loop for this entry
}//RenderHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RenderClu(TObjArray *pClus)
{//used by ReadEvent() and SimulateEvent() to render clusters to polymarker structures, one per chamber
  for(Int_t iCh=0;iCh<=6;iCh++){                                     //chambers loop   
    TClonesArray *pClusCham=(TClonesArray*)pClus->At(iCh);           //get clusters list for this chamber
    fTotClus[iCh] = pClusCham->GetEntries();
    for(Int_t iClu=0;iClu<pClusCham->GetEntries();iClu++){           //clusters loop
      AliHMPIDCluster *pClu = (AliHMPIDCluster*)pClusCham->At(iClu); //get current cluster        
      fRenClu[iCh]->SetNextPoint(pClu->X(),pClu->Y());
//      Printf("RenderClu: ch %i x %f y %f",iCh,pClu->X(),pClu->Y());
    }//cluster loop
  }//chamber loop for this entry
}//RenderClus()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RenderEsd(AliESDEvent *pEsd)
{//used by ReadEvent() or SimulateEvent() to render ESD to polymarker structures for rings and intersections one per chamber
  AliHMPIDRecon rec;
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop to collect cerenkov rings and intersection points
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);    Int_t ch=pTrk->GetHMPIDcluIdx(); //get track and chamber intersected by it
    if(ch<0) continue;                                                          //this track does not intersect any chamber
    Float_t thRa,phRa,xRa,yRa; pTrk->GetHMPIDtrk(xRa,yRa,thRa,phRa);            //get info on current track
    ch/=1000000;                            
    Float_t xPc=0,yPc=0; AliHMPIDTracker::IntTrkCha(pTrk,xPc,yPc);              //find again intersection of track with PC--> it is not stored in ESD!
    fRenTxC[ch]->SetNextPoint(xPc,yPc);                                         //add this intersection point
    Float_t ckov=pTrk->GetHMPIDsignal();                                        //get ckov angle stored for this track  
    if(ckov>0){
      rec.SetTrack(xRa,yRa,thRa,phRa);
      for(Int_t j=0;j<500;j++){
        TVector2 pos; pos=rec.TracePhot(ckov,j*6.28/500.);
       if(!AliHMPIDParam::IsInDead(pos.X(),pos.Y())) fRenRin[ch]->SetNextPoint(pos.X(),pos.Y());
      }      
    }//if ckov is valid
  }//tracks loop  
}//RenEsd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimulateEsd(AliESDEvent *pEsd)
{
  TParticle part; TLorentzVector mom;
  for(Int_t iTrk=0;iTrk<100;iTrk++){//stack loop
    Printf("SimulateEsd: kProton %i",kProton);
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
void SimulateHits(AliESDEvent *pEsd, TClonesArray *pHits)
{//used by SimulateEvent to simulate hits out from provided ESD
  const Int_t kCerenkov=50000050;
  const Int_t kFeedback=50000051;
  
  AliHMPIDRecon rec;
  Float_t eMip=200e-9,ePho=7.5e-9; 
  Int_t hc=0; 
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);
    Float_t xRa,yRa;
    Int_t ch=AliHMPIDTracker::IntTrkCha(pTrk,xRa,yRa);
    if(ch<0) continue; //this track does not hit HMPID
    Float_t beta = pTrk->GetP()/(TMath::Sqrt(pTrk->GetP()*pTrk->GetP()+0.938*0.938));
    Float_t ckov=TMath::ACos(1./(beta*1.292));

    Float_t theta,phi,xPc,yPc,; pTrk->GetHMPIDtrk(xPc,yPc,theta,phi); rec.SetTrack(xRa,yRa,theta,phi); 
    
    if(!AliHMPIDParam::IsInDead(xPc,yPc)) new((*pHits)[hc++]) AliHMPIDHit(ch,eMip,kProton  ,iTrk,xPc,yPc);                 //mip hit
    Int_t nPhots = (Int_t)(20.*TMath::Power(TMath::Sin(ckov),2)/TMath::Power(TMath::Sin(TMath::ACos(1./1.292)),2));
    for(int i=0;i<nPhots;i++){
      TVector2 pos;
      pos=rec.TracePhot(ckov,gRandom->Rndm()*TMath::TwoPi());
      if(!AliHMPIDParam::IsInDead(pos.X(),pos.Y())) new((*pHits)[hc++]) AliHMPIDHit(ch,ePho,kCerenkov,iTrk,pos.X(),pos.Y());
    }                      //photon hits  
    for(int i=0;i<3;i++){//feedback photons
      Float_t x=gRandom->Rndm()*160; Float_t y=gRandom->Rndm()*150;
      if(!AliHMPIDParam::IsInDead(x,y)) new((*pHits)[hc++]) AliHMPIDHit(ch,ePho,kFeedback,iTrk,x,y);                 //feedback hits  
    }//photon hits loop                      
  }//tracks loop    
}//SimulateHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DispTB(TString t0,TString t1,TString t2,TString t3)
{
  fCanvasImp->SetStatusText(t0,0);
  fCanvasImp->SetStatusText(t1,1);
  fCanvasImp->SetStatusText(t2,2);
  fCanvasImp->SetStatusText(t3,3);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DisplayInfo(Int_t px, Int_t py)
{
  if(!gPad) return;
  TVirtualPad *pPad=gPad->GetSelectedPad();
  Int_t padN = pPad->GetNumber()-1;
  if(padN<0 || padN>8) return;
  Int_t ch,pc,padX,padY;
  ch = fChamN[padN];
  if(ch<0) return;
  TString text0,text1,text2,text3;
  
  Float_t x=pPad->AbsPixeltoX(px); Float_t y=pPad->AbsPixeltoY(py); 
  fParam->Lors2Pad(x,y,pc,padX,padY);
  TVector3 xyz=fParam->Lors2Mars(ch,x,y);
  
  if(padX>=0&&padY>=0) {
    text0.Append(Form("Pad(%i,%i)-LORS(%6.2f,%6.2f)-MARS(%7.2f,%7.2f,%7.2f)",padX,padY,x,y,xyz.X(),xyz.Y(),xyz.Z()));
    text1.Append(Form("Module %i Sector %i",ch,pc));
    text2="";
    text3.Append(Form("Pads = %i - Clusters = %i - Multiplicity %5.2f%%",fTotPads[ch],fTotClus[ch],100.*fTotPads[ch]/(144.*160.)));
  }
  
// Find object DIGIT
  TObject *obj = fCanvas->GetSelected();
  TString name = obj->GetName();

  if(name=="") {DispTB(text0,text1,text2,text3);return;}
    
  if(name.Contains("TBox")) {

    TBox *box = (TBox*)obj;
    if(box->GetUniqueID()==0) return; // just black frame hit!
    text2.Append(Form("charge = %i ADC",box->GetUniqueID()));
    
  } else {
//
// Find object HITs CLUSTER TRACK and RING (based on TPolyMarker)
//
    TPolyMarker *b = (TPolyMarker*)obj;

    const Int_t big = 9999;

    // check if point is near one of the points
    Int_t distance = big;
    Int_t index=0;

    Double_t *xPol = b->GetX();
    Double_t *yPol = b->GetY();

    for(Int_t i=0;i<b->Size();i++) {
      Int_t pxp = pPad->XtoAbsPixel(pPad->XtoPad(xPol[i]));
      Int_t pyp = pPad->YtoAbsPixel(pPad->YtoPad(yPol[i]));
      Int_t d   = TMath::Abs(pxp-px) + TMath::Abs(pyp-py);
      if (d < distance) {distance = d; index=i;}
    }

    Int_t type = b->GetUniqueID();
    
    // Case Hit Mip
    
    TString nameHit;
    Int_t iCko[7]={0,0,0,0,0,0,0};
    Int_t iFee[7]={0,0,0,0,0,0,0};
    Int_t iMip[7]={0,0,0,0,0,0,0};
    Int_t iHit;
    Int_t chHit;
    Int_t indHit=0;
    
    switch(type) {
      
      case kHitMip:
      case kHitCko:
      case kHitFee:
        
        for(Int_t iEnt=0;iEnt<fHitTree->GetEntries();iEnt++){    
          fHitTree->GetEntry(iEnt);                              
          for(iHit=0;iHit<fHitLst->GetEntries();iHit++) {       //hits loop
            AliHMPIDHit *pHit = (AliHMPIDHit*)fHitLst->At(iHit);
            chHit = pHit->Ch();
            Int_t pid = pHit->Pid();
            if(pid==50000050) iCko[chHit]++;
            if(pid==50000051) iFee[chHit]++;
            if(pid!=50000050 && pid!=50000051) iMip[chHit]++;
            
            if(type==kHitMip && ch==chHit && index==iMip[chHit]-1) {indHit = iHit;break;}
            if(type==kHitCko && ch==chHit && index==iCko[chHit]-1) {indHit = iHit;break;}
            if(type==kHitFee && ch==chHit && index==iFee[chHit]-1) {indHit = iHit;break;}
          }
        }
        
        fHitTree->GetEntry(6-ch);                              
        AliHMPIDHit *pHit = (AliHMPIDHit*)fHitLst->At(indHit);
        
        Float_t x=pHit->LorsX(); 
        Float_t y=pHit->LorsY();
        Float_t charge = pHit->Q();
        if(pHit->Pid()==50000050) {nameHit = " Cherenkov Photon ";}
        else if(pHit->Pid()==50000051) {nameHit = " Feedback Photon ";}
        else if(fPdg->GetParticle(pHit->Pid())) nameHit = fPdg->GetParticle(pHit->Pid())->GetName();
        text0="";text0.Append(Form("Hit x %6.2f y %6.2f",x,y));
        text2="";text2.Append(Form("Q = %7.2f ADC",charge));
        text3="";text3="Particle: "+ nameHit;
        break;
        
      case kCluster:
        TClonesArray *pClusCham=(TClonesArray*)fCluLst->At(ch);         //get clusters list for this chamber
        AliHMPIDCluster *pClu = (AliHMPIDCluster*)pClusCham->At(index); //get current cluster
        text0="";text0.Append(Form("CLUSTER: x %6.2f y %6.2f",pClu->X(),pClu->Y()));
        text2="";text2.Append(Form("charge = %i ADC",(Int_t)pClu->Q()));
        break;
        
      case kTrack:
        AliESDtrack *pTrk=fEsd->GetTrack(index);
        Float_t thRa,phRa,xRa,yRa; pTrk->GetHMPIDtrk(xRa,yRa,thRa,phRa);
        Float_t xPc,yPc; AliHMPIDTracker::IntTrkCha(pTrk,xPc,yPc);
        text0="";text0.Append(Form("TRACK: x %6.2f y %6.2f at PC plane",xPc,yPc));
        text2="";text2.Append(Form("p = %7.2f GeV/c",pTrk->GetP()));
        Float_t ckov=pTrk->GetHMPIDsignal();                             
        if(ckov>0){
          Float_t x,y;Int_t q,nacc;   pTrk->GetHMPIDmip(x,y,q,nacc);
          text3="";text3.Append(Form("Theta Cherenkov %5.3f with %i photons",pTrk->GetHMPIDsignal(),nacc));
        }      
        break;
      default:
        break;
    }
  }
//Update toolbar status  

  DispTB(text0,text1,text2,text3);  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CloseInfo()
{
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DoZoom(Int_t evt, Int_t px, Int_t py, TObject *)
{
//  Printf(" px %i py%i event %i",px,py,evt);
  if(evt==51) DisplayInfo(px,py);
  if(evt==11)  CloseInfo();
  if(evt!=5 && evt!=6) return; //5- zoom in 6-zoom out
  const Int_t minZoom=64;
  const Int_t maxZoom=2;
  static Int_t zoom[7]={64,64,64,64,64,64,64}; //zoom level
  
 // if(!obj->IsA()->InheritsFrom("TPad")) return;  //current object is not pad
  TVirtualPad *pPad=gPad->GetSelectedPad();
  Int_t padN = pPad->GetNumber()-1;
  if(padN<0 || padN>8) return;
  
  Int_t iCh = fChamN[padN];
  if(iCh < 0) return;
  
  if(evt==5&&zoom[iCh]==maxZoom) return; 
  if(evt==6&&zoom[iCh]==minZoom) return; 
    
  Float_t x=pPad->AbsPixeltoX(px); Float_t y=pPad->AbsPixeltoY(py); 
 
  if(evt==5){ zoom[iCh]=zoom[iCh]/2;     pPad->Range(x-zoom[iCh]*2,y-zoom[iCh]*2,x+zoom[iCh]*2,y+zoom[iCh]*2);} //zoom in
  if(evt==6){ zoom[iCh]=zoom[iCh]*2;     pPad->Range(x-zoom[iCh]*2,y-zoom[iCh]*2,x+zoom[iCh]*2,y+zoom[iCh]*2);} //zoom out 
  if(zoom[iCh]==minZoom) pPad->Range(-10,-10,AliHMPIDParam::SizeAllX()+5,AliHMPIDParam::SizeAllY()+5);
  ((TCanvas *)gTQSender)->SetTitle(Form("zoom x%i",minZoom/zoom[iCh]));
  pPad->Modified();
  pPad->Update();                                              
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ReadEvent()
{//used by NextEvent()  to read curent event and construct all render elements
  if(fNevt && fEvt>=fNevt) fEvt=0; //loop over max event
  if(fNevt && fEvt<0) fEvt=fNevt-1; //loop over max event
  
  if(fHitFile){ 
    if(fHitTree) delete fHitTree;   fHitTree=(TTree*)fHitFile->Get(Form("Event%i/TreeH",fEvt));
    
                          fHitTree->SetBranchAddress("HMPID",&fHitLst);
    for(Int_t iEnt=0;iEnt<fHitTree->GetEntries();iEnt++){    
                          fHitTree->GetEntry(iEnt);                              
      RenderHit(fHitLst);   
    }//prim loop
  }//if hits file

  if(fDigFile){ 
    if(fDigTree) delete fDigTree;   fDigTree=(TTree*)fDigFile->Get(Form("Event%i/TreeD",fEvt));
    for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) fDigTree->SetBranchAddress(Form("HMPID%i",iCh),&(*fDigLst)[iCh]);      
    fDigTree->GetEntry(0);                              
    RenderDig(fDigLst);   
  }//if digs file
    
  if(fCluFile){ 
    if(fCluTree) delete fCluTree; fCluTree=(TTree*)fCluFile->Get(Form("Event%i/TreeR",fEvt));
    for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++) fCluTree->SetBranchAddress(Form("HMPID%i",iCh),&(*fCluLst)[iCh]);      
    fCluTree->GetEntry(0);
    RenderClu(fCluLst);   
  }//if clus file
  
  if(fEsdFile){//if ESD file open
    fEsdTree->GetEntry(fEvt);  
    RenderEsd(fEsd);
  }//if ESD file 
}//ReadEvent()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimulateEvent()
{// used by NextEvent() to simulate all info
  AliCDBManager* pCDB = AliCDBManager::Instance();  pCDB->SetDefaultStorage("local://$HOME"); pCDB->SetRun(0);
  AliCDBEntry *pNmeanEnt=pCDB->Get("HMPID/Calib/Nmean");

  SimulateEsd(fEsd);
  SimulateHits(fEsd,fHitLst);
                 AliHMPIDv2::Hit2Sdi(fHitLst,fSdiLst);                               
          AliHMPIDDigitizer::Sdi2Dig(fSdiLst,fDigLst);     
      AliHMPIDReconstructor::Dig2Clu(fDigLst,fCluLst);
            AliHMPIDTracker::Recon(fEsd,fCluLst,(TObjArray*)pNmeanEnt->GetObject());
            
  RenderHit(fHitLst);
  RenderClu(fCluLst);
  RenderEsd(fEsd);            
}//SimulateEvent()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CheckStatus()
{
  if(fHitFile){
    if(fStHitMip=="OFF") {fHitMipBok->SetFillColor(kRed);} else {fHitMipBok->SetFillColor(18);}
    if(fStHitCko=="OFF") {fHitCkoBok->SetFillColor(kRed);} else {fHitCkoBok->SetFillColor(18);}
    if(fStHitFee=="OFF") {fHitFeeBok->SetFillColor(kRed);} else {fHitFeeBok->SetFillColor(18);}
  }
  if(fDigFile){if(fStDig   =="OFF") {fDigBok->SetFillColor(kRed);} else {fDigBok->SetFillColor(18);}}
  if(fCluFile){if(fStClu   =="OFF") {fCluBok->SetFillColor(kRed);} else {fCluBok->SetFillColor(18);}}
  if(fEsdFile){if(fStEsd   =="OFF") {fEsdBok->SetFillColor(kRed);} else {fEsdBok->SetFillColor(18);}}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GetEvent()
{
  ClearRenders();
  CheckStatus();
  switch(fType){
    case 1: ReadEvent();break;
//    case 2: ReadCosmic(); break;
    case 3: SimulateEvent();     break;
    default:                 return;
  }  
  Draw();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NextEvent()
{
  fTimeArrow = 1;
  fEvt+=fTimeArrow;
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PrevEvent()
{
  fTimeArrow = -1;
  fEvt+=fTimeArrow;
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void End()
{
  gSystem->Exit(0);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hdisp()                                  
{//display events from files if any in current directory or simulated events
  
  fParam=AliHMPIDParam::Instance();                          // first invocation of AliHMPIDParam to initialize geometry...
  fPdg = TDatabasePDG::Instance();                           // first invocation of TDatabasePDG to retrieve particle infos...
  
  CreateContainers();
  CreateRenders();
  
  TString title="Session with";
  if(gSystem->IsFileInIncludePath("HMPID.Hits.root")){// tries to open hits
    fHitFile=TFile::Open("HMPID.Hits.root");       fNevt=fHitFile->GetNkeys(); fType=1; title+=Form(" HITS-%i ",fNevt);
  }
  
  if(gSystem->IsFileInIncludePath("HMPID.Digits.root")){// tries to open clusters
     fDigFile=TFile::Open("HMPID.Digits.root");  fNevt=fDigFile->GetNkeys(); fType=1;  title+=Form(" DIGITS-%i ",fNevt);
   }
  
  if(gSystem->IsFileInIncludePath("HMPID.RecPoints.root")){// tries to open clusters
    fCluFile=TFile::Open("HMPID.RecPoints.root");  fNevt=fCluFile->GetNkeys(); fType=1;  title+=Form(" CLUSTERS-%i ",fNevt);
   }
  
  if(gSystem->IsFileInIncludePath("AliESDs.root")){     
    fEsdFile=TFile::Open("AliESDs.root"); 
    fEsdTree=(TTree*)fEsdFile->Get("esdTree"); 
    fEsd->ReadFromTree(fEsdTree); fEsd->GetStdContent();                       //clm: new ESD schema: see Task Force meeting 20th June, 2007
    fNevt=fEsdTree->GetEntries(); fType=1;  title+=Form(" ESD-%i ",fNevt);
    
    //clm: we need to set the magnetic field  
    if(gSystem->IsFileInIncludePath("galice.root")){
    if(gAlice) delete gAlice;                       
    gAL=AliRunLoader::Open();                       
    gAL->LoadgAlice();           
    AliHMPIDTracker::SetFieldMap(gAL->GetAliRun()->Field(),kTRUE);                         
   }else{
          Printf("=============== NO galice file! Magnetic field for ESD tracking is: %f ===============",AliTracker::GetBz());
     }  
      
  }

  if(gSystem->IsFileInIncludePath("cosmic.root")){                          //clm: Check if cosmic file is in the folder
    fCosFile=TFile::Open("cosmic.root");  fCosTree=(TTree*)fCosFile->Get("cosmic");  fNevt=fCosTree->GetEntries();  fType=2;
    fCosTree->SetBranchAddress("Digs",&fDigLst);   fCosTree->SetBranchAddress("Clus",&fCluLst); 
  }            

  fCanvas=new TCanvas("all","",-1024,768);fCanvas->SetWindowSize(1024,768);
  fCanvasImp = fCanvas->GetCanvasImp();
  fCanvasImp->ShowStatusBar();
  fCanvas->Divide(3,3,0,0);//  pAll->ToggleEditor();
  fCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0,"DoZoom(Int_t,Int_t,Int_t,TObject*)");
  fCanvas->cd(7); 
  gStyle->SetPalette(1,0);
  switch(fType){
    case 1: fCanvas->SetTitle(title.Data());                      
                         TButton *pNxtBtn = new TButton("Next"      ,"NextEvent()"   ,0.0,0.0,0.3,0.1); pNxtBtn->Draw(); 
                         TButton *pPrvBtn = new TButton("Previous"  ,"PrevEvent()"   ,0.3,0.0,0.6,0.1); pPrvBtn->Draw(); 
            if(fHitFile){TButton *pHitBtn = new TButton("Print hits","PrintHits()"   ,0.0,0.2,0.3,0.3); pHitBtn->Draw();}  
            if(fDigFile){TButton *pDigBtn = new TButton("Print digs","PrintDigs()"   ,0.0,0.4,0.3,0.5); pDigBtn->Draw();}  
            if(fCluFile){TButton *pCluBtn = new TButton("Print clus","PrintClus()"   ,0.0,0.6,0.3,0.7); pCluBtn->Draw();} 
            if(fEsdFile){TButton *pEsdBtn = new TButton("Print ESD ","PrintEsd()"    ,0.0,0.8,0.3,0.9); pEsdBtn->Draw();}
            if(fHitFile){      fHitMipBok = new TButton(Form("Mip  %s",fStHitMip.Data()),"SwitchHitMip()",0.3,0.16,0.6,0.22); fHitMipBok->Draw();}
            if(fHitFile){      fHitCkoBok = new TButton(Form("Ckov %s",fStHitCko.Data()),"SwitchHitCko()",0.3,0.22,0.6,0.28); fHitCkoBok->Draw();}  
            if(fHitFile){      fHitFeeBok = new TButton(Form("Fdbk %s",fStHitFee.Data()),"SwitchHitFee()",0.3,0.28,0.6,0.34); fHitFeeBok->Draw();}  
            if(fDigFile){         fDigBok = new TButton(fStDig      ,"SwitchDigs()"  ,0.3,0.4,0.6,0.5); fDigBok->Draw();}  
            if(fCluFile){         fCluBok = new TButton(fStClu      ,"SwitchClus()"  ,0.3,0.6,0.6,0.7); fCluBok->Draw();} 
            if(fEsdFile){         fEsdBok = new TButton(fStEsd      ,"SwitchEsd()"   ,0.3,0.8,0.6,0.9); fEsdBok->Draw();}

            if(fHitFile){TButton *fHitMipOnly = new TButton(" Only ","OnlyHitMip()",0.6,0.16,0.9,0.22); fHitMipOnly->Draw();}  
            if(fHitFile){TButton *fHitCkoOnly = new TButton(" Only ","OnlyHitCko()",0.6,0.22,0.9,0.28); fHitCkoOnly->Draw();}  
            if(fHitFile){TButton *fHitFeeOnly = new TButton(" Only ","OnlyHitFee()",0.6,0.28,0.9,0.34); fHitFeeOnly->Draw();}  
            if(fDigFile){TButton *fDigOnly    = new TButton(" Only ","OnlyDigs()"  ,0.6,0.40,0.9,0.50); fDigOnly->Draw();}  
            if(fCluFile){TButton *fCluOnly    = new TButton(" Only ","OnlyClus()"  ,0.6,0.60,0.9,0.70); fCluOnly->Draw();} 
            if(fEsdFile){TButton *fEsdOnly    = new TButton(" Only ","OnlyEsd()"   ,0.6,0.80,0.9,0.90); fEsdOnly->Draw();}
                         TButton *pAll        = new TButton(" ALL  ","All()"       ,0.5,0.90,0.7,1.00); pAll->Draw();
            
                                     TButton *pEnd    = new TButton("Quit"      ,"End()"         ,0.6,0.0,0.9,0.1); pEnd->Draw();
    break;
    case 2: fCanvas->SetTitle("COSMIC");
                         TButton *pCosBtn = new TButton("Next"      ,"NextEvent()",0.0,0.0,0.3,0.1); pCosBtn->Draw(); 
             
    break;
    case 3: fCanvas->SetTitle("SIMULATION"); 
            TButton *pSimBtn = new TButton("Simulate"  ,"NextEvent()",0.0,0.0,0.3,0.1); pSimBtn->Draw(); 
            TButton *pHitBtn = new TButton("Print hits","PrintHits()",0.0,0.2,0.3,0.3); pHitBtn->Draw();  
            TButton *pDigBtn = new TButton("Print digs","PrintDigs()",0.0,0.4,0.3,0.5); pDigBtn->Draw();
            TButton *pCluBtn = new TButton("Print clus","PrintClus()",0.0,0.6,0.3,0.7); pCluBtn->Draw(); 
            TButton *pEsdBtn = new TButton("Print ESD" ,"PrintEsd()" ,0.0,0.8,0.3,0.9); pEsdBtn->Draw();
            if(fHitFile){      fHitMipBok = new TButton(Form("Mips %s",fStHitMip.Data()),"SwitchHitMip()",0.3,0.16,0.6,0.22); fHitMipBok->Draw();}  
            if(fHitFile){      fHitCkoBok = new TButton(Form("Ckov %s",fStHitCko.Data()),"SwitchHitCko()",0.3,0.22,0.6,0.28); fHitCkoBok->Draw();}  
            if(fHitFile){      fHitFeeBok = new TButton(Form("Fdbk %s",fStHitFee.Data()),"SwitchHitFee()",0.3,0.28,0.6,0.34); fHitFeeBok->Draw();}  
            if(fDigFile){         fDigBok = new TButton(fStDig      ,"SwitchDigs()"  ,0.3,0.4,0.6,0.5); fDigBok->Draw();}  
            if(fCluFile){         fCluBok = new TButton(fStClu      ,"SwitchClus()"  ,0.3,0.6,0.6,0.7); fCluBok->Draw();} 
            if(fEsdFile){         fEsdBok = new TButton(fStEsd      ,"SwitchEsd()"   ,0.3,0.8,0.6,0.9); fEsdBok->Draw();}
            
            if(fHitFile){TButton *fHitMipOnly = new TButton(" Only ","OnlyHitMip()",0.6,0.16,0.9,0.22); fHitMipOnly->Draw();}  
            if(fHitFile){TButton *fHitCkoOnly = new TButton(" Only ","OnlyHitCko()",0.6,0.22,0.9,0.28); fHitCkoOnly->Draw();}  
            if(fHitFile){TButton *fHitFeeOnly = new TButton(" Only ","OnlyHitFee()",0.6,0.28,0.9,0.34); fHitFeeOnly->Draw();}  
            if(fDigFile){TButton *fDigOnly    = new TButton(" Only ","OnlyDigs()"  ,0.6,0.40,0.9,0.50); fDigOnly->Draw();}  
            if(fCluFile){TButton *fCluOnly    = new TButton(" Only ","OnlyClus()"  ,0.6,0.60,0.9,0.70); fCluOnly->Draw();} 
            if(fEsdFile){TButton *fEsdOnly    = new TButton(" Only ","OnlyEsd()"   ,0.6,0.80,0.9,0.90); fEsdOnly->Draw();}
                         TButton *pSimAll     = new TButton(" ALL  ","All()"       ,0.5,0.90,0.7,1.00); pSimAll->Draw();
    break;
  }
    
  NextEvent();           
      
}//Hdisp()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SwitchHitMip()
{
  fStHitMip=(fStHitMip=="OFF")? "ON":"OFF";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SwitchHitCko()
{
  fStHitCko=(fStHitCko=="OFF")? "ON":"OFF";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SwitchHitFee()
{
  fStHitFee=(fStHitFee=="OFF")? "ON":"OFF";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SwitchDigs()
{
  fStDig=(fStDig=="OFF")? "ON":"OFF";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SwitchClus()
{
  fStClu=(fStClu=="OFF")? "ON":"OFF";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SwitchEsd()
{
  fStEsd=(fStEsd=="OFF")? "ON":"OFF";
   GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AllNo()
{
  if(fHitFile) {
    fStHitMip = "OFF";
    fStHitCko = "OFF";
    fStHitFee = "OFF";
  }
  if(fDigFile) fStDig    = "OFF";
  if(fCluFile) fStClu    = "OFF";
  if(fEsdFile) fStEsd    = "OFF";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AllYes()
{
  if(fHitFile) {
    fStHitMip = "ON";
    fStHitCko = "ON";
    fStHitFee = "ON";
  }
  if(fDigFile) fStDig    = "ON";
  if(fCluFile) fStClu    = "ON";
  if(fEsdFile) fStEsd    = "ON";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void All()
{
  AllYes();
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OnlyHitMip()
{
  AllNo();
  fStHitMip = "ON";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OnlyHitCko()
{
  AllNo();
  fStHitCko = "ON";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OnlyHitFee()
{
  AllNo();
  fStHitFee = "ON";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OnlyDigs()
{
  AllNo();
  fStDig    = "ON";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OnlyClus()
{
  AllNo();
  fStClu    = "ON";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OnlyEsd()
{
  AllNo();
  fStEsd    = "ON";
  GetEvent();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
