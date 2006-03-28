//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliRICH.h"
#include "AliRICHParam.h"
#include "AliRICHHelix.h" //ReadESD
#include <TFile.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <AliStack.h>
#include <AliRun.h>
#include <AliMC.h>       //ctor
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliGenHijingEventHeader.h>
#include <AliESD.h>
#include <TH1F.h>        //HitQA()
#include <TH2F.h>        //Display() 
#include <TBenchmark.h>
#include <AliLog.h>
#include <TLatex.h>      //Display()
#include <TCanvas.h>     //Display()
#include <TGraph.h>      //Display()
#include <TStyle.h>      //Display()
#include <TMarker.h>     //Display()
ClassImp(AliRICH)    
//__________________________________________________________________________________________________
// RICH manager class   
//BEGIN_HTML
/*
  <img src="gif/alirich.gif">
*/
//END_HTML
//__________________________________________________________________________________________________
AliRICH::AliRICH():AliDetector(),fSdig(0),fSdigCnt(0),fDig(0),fClu(0) 
{
//Default ctor should not contain any new operators
//AliDetector ctor deals with Hits and Digits  
}//AliRICH::AliRICH()
//__________________________________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title),fSdig(0),fSdigCnt(0),fDig(0),fClu(0)
{
//Named ctor
  AliDebug(1,"Start.");
//AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  HitCreate();          gAlice->GetMCApp()->AddHitList(fHits);
  fNcham=7;
  fCounters.ResizeTo(40); fCounters.Zero();
  AliDebug(1,"Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{
//dtor
  AliDebug(1,"Start.");

  
  if(fHits)      delete fHits;
  if(fSdig)      delete fSdig;
  if(fDigits)    delete fDigits;
  if(fDig)      {fDig->Delete();   delete fDig;}
  if(fClu)      {fClu->Delete();   delete fClu;}
  AliDebug(1,"Stop.");    
}//AliRICH::~AliRICH()
//__________________________________________________________________________________________________
void AliRICH::BuildGeometry() 
{
//Builds a TNode geometry for event display
  AliDebug(1,"Start.");
  
  AliDebug(1,"Stop.");    
}//void AliRICH::BuildGeometry()
//__________________________________________________________________________________________________
void AliRICH::MakeBranch(Option_t* option)
{
//Create Tree branches for the RICH.
  AliDebug(1,Form("Start with option= %s.",option));
    
  const Int_t kBufferSize = 4000;
      
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&TreeH()){//H
    HitCreate();      //branch will be created in AliDetector::MakeBranch
  }//H     
  AliDetector::MakeBranch(option);//this is after cH because we need to guarantee that fHits array is created
      
  if(cS&&fLoader->TreeS()){//S  
    SDigCreate();   MakeBranchInTree(fLoader->TreeS(),"RICH",&fSdig,kBufferSize,0) ;
  }//S
   
  if(cD&&fLoader->TreeD()){//D
    DigCreate();
    for(Int_t i=0;i<fNcham;i++){ 
      MakeBranchInTree(fLoader->TreeD(),Form("%s%d",GetName(),i+1),&((*fDig)[i]),kBufferSize,0);
    }
  }//D
  
  if(cR&&fLoader->TreeR()){//R
    CluCreate();
    for(Int_t i=0;i<fNcham;i++)
      MakeBranchInTree(fLoader->TreeR(),Form("%sClusters%d",GetName(),i+1), &((*fClu)[i]), kBufferSize, 0);    
  }//R
  AliDebug(1,"Stop.");   
}//void AliRICH::MakeBranch(Option_t* option)
//__________________________________________________________________________________________________
void AliRICH::SetTreeAddress()
{
//Set branch address for the Hits and Digits Tree.
  AliDebug(1,"Start.");
      
  TBranch *branch;
    
  if(fLoader->TreeH()){//H
    AliDebug(1,"tree H is requested.");
    HitCreate();//branch map will be in AliDetector::SetTreeAddress    
  }//H
  AliDetector::SetTreeAddress();//this is after TreeH because we need to guarantee that fHits array is created

  if(fLoader->TreeS()){//S
    AliDebug(1,"tree S is requested.");
    branch=fLoader->TreeS()->GetBranch(GetName());        if(branch){SDigCreate();   branch->SetAddress(&fSdig);}
  }//S
    
  if(fLoader->TreeD()){//D    
    AliDebug(1,"tree D is requested.");
    for(int i=0;i<fNcham;i++){   branch=fLoader->TreeD()->GetBranch(Form("%s%d",GetName(),i+1));  if(branch){DigCreate(); branch->SetAddress(&((*fDig)[i]));}}
  }//D
    
  if(fLoader->TreeR()){//R
    AliDebug(1,"tree R is requested.");
    for(int i=0;i<fNcham;i++){   branch=fLoader->TreeR()->GetBranch(Form("%sClusters%d" ,GetName(),i+1));  if(branch){CluCreate(); branch->SetAddress(&((*fClu)[i]));}}
  }//R
  AliDebug(1,"Stop.");
}//void AliRICH::SetTreeAddress()
//__________________________________________________________________________________________________
// AliRICHHit* AliRICH::Hit(Int_t tid)const
// {
// // Search for the first RICH hit belonging to the given tid
//   GetLoader()->LoadHits();
//   for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop      
//     GetLoader()->TreeH()->GetEntry(iPrimN);
//     for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
//       AliRICHHit *pHit=(AliRICHHit*)Hits()->At(iHitN);
//       if(tid==pHit->Track()) {GetLoader()->UnloadHits();return pHit;}
//     }//hits
//   }//prims loop
//   GetLoader()->UnloadHits();
//   return 0;
// }
//__________________________________________________________________________________________________
void AliRICH::HitPrint(Int_t iEvtN)const
{
//Prints a list of RICH hits for a given event. Default is event number 0.
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  AliInfo(Form("List of RICH hits for event %i",iEvtN));
  if(GetLoader()->LoadHits()) return;
  
  Int_t iTotalHits=0;
  for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
    GetLoader()->TreeH()->GetEntry(iPrimN);      
    Hits()->Print();
    iTotalHits+=Hits()->GetEntries();
  }
  GetLoader()->UnloadHits();
  AliInfo(Form("totally %i hits",iTotalHits));
}
//__________________________________________________________________________________________________
void AliRICH::SDigPrint(Int_t iEvtN)const
{
//prints a list of RICH sdigits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  Info("PrintSDigits","List of RICH sdigits for event %i",iEvtN);
  if(GetLoader()->LoadSDigits()) return;
  
  GetLoader()->TreeS()->GetEntry(0);
  fSdig->Print();
  GetLoader()->UnloadSDigits();
  Printf("totally %i sdigits",fSdig->GetEntries());
}
//__________________________________________________________________________________________________
void AliRICH::DigPrint(Int_t iEvtN)const
{
//prints a list of RICH digits  for a given event
  if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
  Printf("List of RICH digits for event %i",iEvtN);
  if(GetLoader()->LoadDigits()) return;
  
  Int_t iTotalDigits=0;
  GetLoader()->TreeD()->GetEntry(0);
  if(!fDig) return;
  for(Int_t iCham=0;iCham<fNcham;iCham++){
    TClonesArray *pDigs=(TClonesArray*)fDig->At(iCham);    iTotalDigits+=pDigs->GetEntries();    pDigs->Print();
  }
  GetLoader()->UnloadDigits();
  Printf("totally %i Digits",iTotalDigits);
}
//__________________________________________________________________________________________________
void AliRICH::OccupancyPrint(Int_t iEvtNreq)const
{
//prints occupancy for each chamber in a given event
  Int_t iEvtNmin,iEvtNmax;
  if(iEvtNreq==-1){
    iEvtNmin=0;
    iEvtNmax=gAlice->GetEventsPerRun();
  } else { 
    iEvtNmin=iEvtNreq;iEvtNmax=iEvtNreq+1;
  }
    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;    
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;    
  
//  Info("Occupancy","for event %i",iEvtN);
  if(GetLoader()->LoadHits()) return;
  if(GetLoader()->LoadDigits()) return;

  Int_t totPadsPerChamber = AliRICHParam::NpadsX()*AliRICHParam::NpadsY();  

  
  for(Int_t iEvtN=iEvtNmin;iEvtN<iEvtNmax;iEvtN++){    
    Int_t nDigCh[kNchambers]={0,0,0,0,0,0,0};  
    Int_t iChHits[kNchambers]={0,0,0,0,0,0,0};
    Int_t nPrim[kNchambers]={0,0,0,0,0,0,0};
    Int_t nSec[kNchambers]={0,0,0,0,0,0,0};
    AliInfo(Form("events processed %i",iEvtN));
    if(GetLoader()->GetRunLoader()->GetEvent(iEvtN)) return;    
    AliStack *pStack = GetLoader()->GetRunLoader()->Stack();
    for(Int_t iPrimN=0;iPrimN<GetLoader()->TreeH()->GetEntries();iPrimN++){//prims loop
      GetLoader()->TreeH()->GetEntry(iPrimN);      
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){
        AliRICHHit *pHit = (AliRICHHit *)Hits()->At(iHitN);
        if(pHit->Eloss()>0){
          iChHits[pHit->C()-1]++;
          if(pStack->Particle(pHit->GetTrack())->Rho()<0.01) nPrim[pHit->C()-1]++;else nSec[pHit->C()-1]++;
        }
      }
    }
    GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iCh=0;iCh<fNcham;iCh++) {
      nDigCh[iCh]= ((TClonesArray*)fDig->At(iCh))->GetEntries();
      Double_t occupancy = (Double_t)nDigCh[iCh-1]/(Double_t)totPadsPerChamber;
      Info("Occupancy","for chamber %i = %4.2f %% and charged prim tracks %i and sec. tracks %i with total %i",
        iCh+1,occupancy*100.,nPrim[iCh],nSec[iCh],iChHits[iCh]);
    }
  }
  GetLoader()->UnloadHits();
  GetLoader()->UnloadDigits();
  GetLoader()->GetRunLoader()->UnloadHeader();    
  GetLoader()->GetRunLoader()->UnloadKinematics();    
}
//__________________________________________________________________________________________________
void AliRICH::CluPrint(Int_t iEvtN)const
{
//prints a list of RICH clusters  for a given event
  Printf("List of RICH clusters for event %i",iEvtN);
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->LoadRecPoints()) return;
  
  Int_t iCluCnt=0;
  GetLoader()->TreeR()->GetEntry(0);
  for(Int_t iCham=0;iCham<fNcham;iCham++){
    TClonesArray *pClus=(TClonesArray*)fClu->At(iCham);    iCluCnt+=pClus->GetEntries();    pClus->Print();
  }
  GetLoader()->UnloadRecPoints();
  Printf("totally %i clusters for event %i",iCluCnt,iEvtN);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::DisplayEvent(Int_t iEvtNmin,Int_t iEvtNmax)const
{
// Display digits, reconstructed tracks intersections and RICH rings if available 
  TH2F *pH2[8];

  GetLoader()->LoadDigits();
  
  TLatex t;  t.SetTextSize(0.1);
  TCanvas *pC = new TCanvas("RICHDisplay","RICH Display",0,0,1226,900);  pC->Divide(3,3);  pC->cd(9); t.DrawText(0.2,0.4,"View to IP");  
  gStyle->SetPalette(1);

  
  for(Int_t iCh=1;iCh<=fNcham;iCh++) {
    pH2[iCh] = new TH2F(Form("RichDigH2_%i",iCh),Form("Chamber %i;cm;cm",iCh),165,0,AliRICHParam::PcSizeX(),144,0,AliRICHParam::PcSizeY());
    pH2[iCh]->SetMarkerColor(kGreen); 
    pH2[iCh]->SetMarkerStyle(29); 
    pH2[iCh]->SetMarkerSize(0.4);
    pH2[iCh]->SetStats(kFALSE);
    pH2[iCh]->SetMaximum(300);
  }
  
  if(iEvtNmax>gAlice->GetEventsPerRun()||iEvtNmax==0) iEvtNmax=gAlice->GetEventsPerRun()-1;

  for(Int_t iEvt=iEvtNmin;iEvt<=iEvtNmax;iEvt++) {//events loop
    pC->cd(3);  t.DrawText(0.2,0.4,Form("Event %i",iEvt));        

    GetLoader()->GetRunLoader()->GetEvent(iEvt); //get event
    GetLoader()->TreeD()->GetEntry(0);           //get list of digits 
    for(Int_t iCh=1;iCh<=fNcham;iCh++) {//chambers loop
      pH2[iCh]->Reset();    
      for(Int_t iDig=0;iDig < Digs(iCh)->GetEntries();iDig++) {//digits loop
        AliRICHDigit *pDig = (AliRICHDigit*)Digs(iCh)->At(iDig);
        TVector2 x2=AliRICHParam::Pad2Loc(pDig->Pad());
        pH2[pDig->C()]->Fill(x2.X(),x2.Y(),pDig->Qdc());
      }//digits loop
      if(iCh==1) pC->cd(9);
      if(iCh==2) pC->cd(8);
      if(iCh==3) pC->cd(6);
      if(iCh==4) pC->cd(5);
      if(iCh==5) pC->cd(4);
      if(iCh==6) pC->cd(2);
      if(iCh==7) pC->cd(1);
      pH2[iCh]->Draw("col");
      ReadESD(iEvt,iCh);
      AliRICHParam::DrawSectors();
    }//chambers loop
    pC->Update();
    pC->Modified();

    if(iEvt<iEvtNmax) {gPad->WaitPrimitive();pC->Clear();}
  }//events loop
}//ShowEvent()
//__________________________________________________________________________________________________
void AliRICH::Display()const
{
//Provides fast event display
//For RICH only, full display is .x Display.C    
  Bool_t isHits    =!GetLoader()->LoadHits();
  Bool_t isDigits  =!GetLoader()->LoadDigits();
  Bool_t isClusters=!GetLoader()->LoadRecPoints();
  
  if(!isHits && !isDigits && !isClusters){Error("Exec","No hits digits and clusters. Nothing to display.");return;}
  
  TCanvas *pCanvas = new TCanvas("Display","RICH Display",0,0,600,600);
  
  TH2F *pHitsH2=0,*pDigitsH2=0,*pClustersH2=0;
  
  if(isHits)     pHitsH2     = new TH2F("pHitsH2"  ,  "Event Display;x,cm;y,cm",165,0,AliRICHParam::PcSizeX(),
                                                                                144,0,AliRICHParam::PcSizeY());
  if(pHitsH2)    pHitsH2->SetStats(kFALSE);
  
  if(isDigits)   pDigitsH2   = new TH2F("pDigitsH2"  ,"Event Display",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  if(isClusters) pClustersH2 = new TH2F("pClustersH2","Event Display",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  
  for(Int_t iEvt=0;iEvt<GetLoader()->GetRunLoader()->GetNumberOfEvents();iEvt++){//events Loop
    GetLoader()->GetRunLoader()->GetEvent(iEvt);  
//display all the staff on chamber by chamber basis           
    for(Int_t iCh=1;iCh<=fNcham;iCh++){//chambers loop       
      if(isHits)     pHitsH2    ->Reset();     
      if(isDigits)   pDigitsH2  ->Reset();     
      if(isClusters) pClustersH2->Reset();
//deals with hits
      for(Int_t i=0;i<GetLoader()->TreeH()->GetEntries();i++){//TreeH loop
        GetLoader()->TreeH()->GetEntry(i);
        for(Int_t iHit=0;iHit<Hits()->GetEntries();iHit++){//hits loop
          AliRICHHit *pHit = (AliRICHHit*)Hits()->At(iHit);
          if(pHit->C()==iCh){
            TVector2 hitLocX2 = AliRICHParam::Instance()->Mars2Lors(iCh,pHit->OutX3());
            pHitsH2->Fill(hitLocX2.X(),hitLocX2.Y(),200);
          }//if
        }//hits loop         
      }//TreeH loop
      pHitsH2->SetTitle(Form("event %i chamber %2i",iEvt,iCh));
      pHitsH2->SetMarkerColor(kRed); pHitsH2->SetMarkerStyle(29); pHitsH2->SetMarkerSize(0.4);
      pHitsH2->Draw();
      AliRICHParam::DrawSectors();
      TLatex l; l.SetNDC(); l.SetTextSize(0.02);
      if(!isHits)     {l.SetTextColor(kRed)  ;l.DrawLatex(0.1,0.01,"No Hits"    );}
      if(!isDigits)   {l.SetTextColor(kGreen);l.DrawLatex(0.4,0.01,"No DIGITS"  );}
      if(!isClusters) {l.SetTextColor(kBlue) ;l.DrawLatex(0.8,0.01,"No CLUSTERS");}
      pCanvas->Update();        pCanvas->Modified();       gPad->WaitPrimitive();
//deals with digits      
      if(isDigits){
        GetLoader()->TreeD()->GetEntry(0);
        for(Int_t iDig=0;iDig < Digs(iCh)->GetEntries();iDig++){//digits loop
          AliRICHDigit *pDig = (AliRICHDigit*)Digs(iCh)->At(iDig);
	        TVector2 x2=AliRICHParam::Pad2Loc(pDig->Pad());
	        pDigitsH2->Fill(x2.X(),x2.Y(),100);
        }//digits loop
        pDigitsH2->SetMarkerColor(kGreen); pDigitsH2->SetMarkerStyle(29); pDigitsH2->SetMarkerSize(0.4);
        pDigitsH2->Draw("same");
        pCanvas->Update();        pCanvas->Modified();       gPad->WaitPrimitive();
      }//if(isDigits)      
//deals with clusters      
      if(isClusters){
        GetLoader()->TreeR()->GetEntry(0);
        for(Int_t iClu=0;iClu<Clus(iCh)->GetEntries();iClu++){//clusters loop
          AliRICHCluster *pClu = (AliRICHCluster*)Clus(iCh)->At(iClu);
          pClustersH2->Fill(pClu->X(),pClu->Y(),50);
        }//clusters loop
        pClustersH2->SetMarkerColor(kBlue); pClustersH2->SetMarkerStyle(29);  pClustersH2->SetMarkerSize(0.4);
        pClustersH2->Draw("same");
        pCanvas->Update();        pCanvas->Modified();       gPad->WaitPrimitive();
      }//if(isClusters)
    }//chambers loop
  }//events Loop
  
  delete pCanvas;
  GetLoader()->UnloadHits();
  if(isDigits)   GetLoader()->UnloadDigits();
  if(isClusters) GetLoader()->UnloadRecPoints();
}//Display()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICH::ReadESD(Int_t iEventN, Int_t iChamber)const
{
//
  TFile *pFile=TFile::Open("AliESDs.root","read");
  if(!pFile || !pFile->IsOpen()) {AliInfo("ESD file not open.");return;}      //open AliESDs.root                                                                    
  TTree *pTree = (TTree*) pFile->Get("esdTree");
  if(!pTree){AliInfo("ESD not found.");return;}                               //get ESD tree
  
                                                                 
  AliESD *pESD=new AliESD;  pTree->SetBranchAddress("ESD", &pESD);
  
  pTree->GetEvent(iEventN);
  
  Double_t b = pESD->GetMagneticField()/10.;
  
  Int_t iNtracks=pESD->GetNumberOfTracks();    
  
  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
    Int_t charge = (Int_t)(-TMath::Sign(1.,pTrack->GetSign()*b));
    AliRICHHelix helix(pTrack->X3(),pTrack->P3(),charge,b);
    Int_t iChamberOnRICH=helix.RichIntersect(AliRICHParam::Instance());        
    if(iChamberOnRICH==iChamber) {
      TMarker *trackImpact = new TMarker(helix.PosPc().X(),helix.PosPc().Y(),kStar);
      trackImpact->SetMarkerColor(kRed);
      trackImpact->Draw();
//
      Int_t iChamberRecon = pTrack->GetRICHcluster()/1000000;
      if(iChamberRecon==iChamber) {
        Double_t thetaCer = pTrack->GetRICHsignal();
        if(thetaCer<0) continue;
        TVector3 entrance(helix.PosRad().X(),helix.PosRad().Y(),0);
        Double_t thetaTrack,phiTrack;
        pTrack->GetRICHthetaPhi(thetaTrack,phiTrack);
        TVector3 vectorTrack;
        vectorTrack.SetMagThetaPhi(pTrack->GetP(),thetaTrack,phiTrack);
        AliInfo(Form("Draw ring started for track %i on chamber %i",iTrackN,iChamber));
        AliInfo(Form("ThetaCer %f TrackTheta %f TrackPhi %f Momentum %f",thetaCer,thetaTrack,phiTrack,pTrack->GetP()));
        Double_t dx,dy;
        pTrack->GetRICHdxdy(dx,dy);
        DrawRing(entrance,vectorTrack,thetaCer);
      }
    }
  }
  delete pESD;  pFile->Close();//close AliESDs.root
}
//__________________________________________________________________________________________________
void AliRICH::DrawRing(TVector3 entrance,TVector3 vectorTrack,Double_t thetaCer)const
{
  Double_t xGraph[100],yGraph[100];
  Int_t nPointsToDraw = 0;
  for(Int_t i=0;i<100;i++) {
    Double_t phiCer = 2*TMath::Pi()*i/100;
    TVector3 pos = AliRICHParam::Instance()->ForwardTracing(entrance,vectorTrack,thetaCer,phiCer);
    if(pos.X()==-999) continue;
    xGraph[nPointsToDraw] = pos.X();yGraph[nPointsToDraw] = pos.Y();nPointsToDraw++;
  }
//  AliInfo(Form("Npoints per ring %i",nPointsToDraw));
  TGraph *gra = new TGraph(nPointsToDraw,xGraph,yGraph);
  gra->Draw("C");  
}
//__________________________________________________________________________________________________
void AliRICH::SummaryOfEvent(Int_t iEvtN) const
{
//prints a summary for a given event
  AliInfo(Form("Summary of event %i",iEvtN));
  GetLoader()->GetRunLoader()->GetEvent(iEvtN);    
  if(GetLoader()->GetRunLoader()->LoadHeader()) return;
  if(GetLoader()->GetRunLoader()->LoadKinematics()) return;
  AliStack *pStack=GetLoader()->GetRunLoader()->Stack();
  
  AliGenEventHeader* pGenHeader =  gAlice->GetHeader()->GenEventHeader();
  if(pGenHeader->InheritsFrom("AliGenHijingEventHeader")) {
    AliInfo(Form(" Hijing event with impact parameter b = %.2f (fm)",((AliGenHijingEventHeader*) pGenHeader)->ImpactParameter()));
  }
  Int_t nChargedPrimaries=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) {
    TParticle *pParticle = pStack->Particle(i);
    if(pParticle->IsPrimary()&&pParticle->GetPDG()->Charge()!=0) nChargedPrimaries++;
    }
  AliInfo(Form("Total number of         primaries %i",pStack->GetNprimary()));
  AliInfo(Form("Total number of charged primaries %i",nChargedPrimaries));
  AliInfo(Form("Total n. of tracks in stack(+sec) %i",pStack->GetNtrack()));
  GetLoader()->GetRunLoader()->UnloadHeader();
  GetLoader()->GetRunLoader()->UnloadKinematics();
}
//__________________________________________________________________________________________________
void AliRICH::HitQA(Double_t cut,Double_t cutele,Double_t cutR)
{
// Provides a set of control plots intended primarily for charged particle flux analisys
// Arguments: cut (GeV)    - cut on momentum of any charged particles but electrons, 
//            cetele (GeV) - the same for electrons-positrons
//            cutR (cm)    - cut on production vertex radius (cylindrical system)        
  gBenchmark->Start("HitsAna");
  
  Double_t cutPantiproton    =cut;
  Double_t cutPkaonminus     =cut;
  Double_t cutPpionminus     =cut;
  Double_t cutPmuonminus     =cut;
  Double_t cutPpositron      =cutele;
                    
  Double_t cutPelectron      =cutele;
  Double_t cutPmuonplus      =cut;
  Double_t cutPpionplus      =cut;
  Double_t cutPkaonplus      =cut;
  Double_t cutPproton        =cut;
                       

  TH2F *pEleHitRZ    =new TH2F("EleHitRZ"    ,Form("e^{+} e^{-} hit %s;z[cm];R[cm]" ,GetName())     , 400,-300,300 ,400,-500,500);   //R-z plot 0cm<R<550cm -300cm<z<300cm  
  TH2F *pEleHitRP    =new TH2F("EleHitRP"    ,Form("e^{+} e^{-} hit %s;p[GeV];R[cm]",GetName())     ,1000,-1  ,1   ,400,   0,550);   //R-p plot 0cm<R<550cm -1GeV<p<1GeV 
  TH1F *pEleAllP     =new TH1F("EleAllP"     ,     "e^{+} e^{-} all;p[GeV]"                         ,1000,-1  ,1                );  
  TH1F *pEleHitP     =new TH1F("EleHitP"     ,Form("e^{+} e^{-} hit %s;p[GeV]"      ,GetName())     ,1000,-1  ,1                );   
  TH1F *pMuoHitP     =new TH1F("MuoHitP"     ,Form("#mu^{-} #mu^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pPioHitP     =new TH1F("PioHitP"     ,Form("#pi^{-} #pi^{+} hit %s;p[GeV]"  ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pKaoHitP     =new TH1F("KaoHitP"     ,Form("K^{-} K^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH1F *pProHitP     =new TH1F("ProHitP"     ,Form("p^{-} p^{+} hit %s;p[GeV]"      ,GetName())     ,1000,-4  ,4                ); 
  TH2F *pFlux        =new TH2F("flux"        ,Form("%s flux with Rvertex<%.1fcm"    ,GetName(),cutR),10  ,-5  ,5   , 10,0    ,10); //special text hist
  TH2F *pVertex      =new TH2F("vertex"      ,Form("%s 2D vertex of RICH hit;x;y"   ,GetName())     ,120 ,0   ,600 ,120,0    ,600); //special text hist
  TH1F *pRho         =new TH1F("rho"         ,Form("%s r of RICH hit"               ,GetName())     ,600 ,0   ,600); //special text hist
  pFlux->SetStats(0);
  pFlux->GetXaxis()->SetBinLabel(1 ,Form("p^{-}>%.3fGeV/c"   ,cutPantiproton));        
  pFlux->GetXaxis()->SetBinLabel(2 ,Form("K^{-}>%.3fGeV/c"   ,cutPkaonminus ));        
  pFlux->GetXaxis()->SetBinLabel(3 ,Form("#pi^{-}>%.3fGeV/c" ,cutPpionminus ));      
  pFlux->GetXaxis()->SetBinLabel(4 ,Form("#mu^{-}>%.3fGeV/c" ,cutPmuonminus ));      
  pFlux->GetXaxis()->SetBinLabel(5 ,Form("e^{+}>%.3fGeV/c"   ,cutPpositron  ));        
  
  pFlux->GetXaxis()->SetBinLabel(6 ,Form("e^{-}>%.3fGeV/c"   ,cutPelectron  ));        
  pFlux->GetXaxis()->SetBinLabel(7 ,Form("#mu^{+}>%.3fGeV/c" ,cutPmuonplus  ));      
  pFlux->GetXaxis()->SetBinLabel(8 ,Form("#pi^{+}>%.3fGeV/c" ,cutPpionplus  ));      
  pFlux->GetXaxis()->SetBinLabel(9 ,Form("K^{+}>%.3fGeV/c"   ,cutPkaonplus  ));        
  pFlux->GetXaxis()->SetBinLabel(10,Form("p^{+}>%.3fGeV/c"   ,cutPproton    ));        
  
  pFlux->GetYaxis()->SetBinLabel(1,"sum");  
  pFlux->GetYaxis()->SetBinLabel(2,"ch1");  
  pFlux->GetYaxis()->SetBinLabel(3,"ch2");  
  pFlux->GetYaxis()->SetBinLabel(4,"ch3");  
  pFlux->GetYaxis()->SetBinLabel(5,"ch4");  
  pFlux->GetYaxis()->SetBinLabel(6,"ch5");  
  pFlux->GetYaxis()->SetBinLabel(7,"ch6");  
  pFlux->GetYaxis()->SetBinLabel(8,"ch7");  
  pFlux->GetYaxis()->SetBinLabel(9,"prim"); 
  pFlux->GetYaxis()->SetBinLabel(10,"tot");  
  
//end of hists definition
   
  Int_t iNevents=fLoader->GetRunLoader()->GetAliRun()->GetEventsPerRun(),iCntPrimParts=0,iCntTotParts=0;
//load all needed trees   
  fLoader->LoadHits(); 
  fLoader->GetRunLoader()->LoadHeader(); 
  fLoader->GetRunLoader()->LoadKinematics();  
  
  for(Int_t iEvtN=0;iEvtN < iNevents;iEvtN++){//events loop
    fLoader->GetRunLoader()->GetEvent(iEvtN);
    AliInfo(Form(" %i event processes",fLoader->GetRunLoader()->GetEventNumber()));
    AliStack *pStack= fLoader->GetRunLoader()->Stack(); 
    
    for(Int_t iParticleN=0;iParticleN<pStack->GetNtrack();iParticleN++){//stack loop
      TParticle *pPart=pStack->Particle(iParticleN);

      if(iParticleN%10000==0) AliInfo(Form(" %i particles read",iParticleN));
    
      switch(pPart->GetPdgCode()){
        case kProtonBar: pFlux->Fill(-4.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-4.5,8); break;
        case kKMinus:    pFlux->Fill(-3.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-3.5,8); break;
        case kPiMinus:   pFlux->Fill(-2.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-2.5,8); break;
        case kMuonMinus: pFlux->Fill(-1.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-1.5,8); break;
        case kPositron:  pFlux->Fill(-0.5,9); if(pPart->Rho()<0.01) pFlux->Fill(-0.5,8); pEleAllP->Fill(-pPart->P()); break;
      
        case kElectron:  pFlux->Fill( 0.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 0.5,8); pEleAllP->Fill( pPart->P()); break;      
        case kMuonPlus:  pFlux->Fill( 1.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 1.5,8); break;      
        case kPiPlus:    pFlux->Fill( 2.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 2.5,8); break;      
        case kKPlus:     pFlux->Fill( 3.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 3.5,8); break;      
        case kProton:    pFlux->Fill( 4.5,9); if(pPart->Rho()<0.01) pFlux->Fill( 4.5,8); break;            
      }//switch
    }//stack loop
//now hits analiser        
    for(Int_t iEntryN=0;iEntryN < fLoader->TreeH()->GetEntries();iEntryN++){//TreeH loop
      fLoader->TreeH()->GetEntry(iEntryN);                                  //get current entry (prim)                
      for(Int_t iHitN=0;iHitN < Hits()->GetEntries();iHitN++){//hits loop
        AliRICHHit *pHit = (AliRICHHit*)Hits()->At(iHitN);            //get current hit
        TParticle  *pPart=pStack->Particle(pHit->GetTrack());      //get stack particle which produced the current hit
        
        if(pPart->GetPDG()->Charge()!=0&&pPart->Rho()>0.1) pVertex->Fill(pPart->Vx(),pPart->Vy()); //safe margin for sec.
        if(pPart->GetPDG()->Charge()!=0) pRho->Fill(pPart->Rho()); //safe margin for sec.
        if(pPart->R()>cutR) continue;                                   //cut on production radius (cylindrical system) 
      
        switch(pPart->GetPdgCode()){
          case kProtonBar: if(pPart->P()>cutPantiproton) {pProHitP->Fill(-pPart->P()); pFlux->Fill(-4.5,pHit->C());}break;
          case kKMinus   : if(pPart->P()>cutPkaonminus)  {pKaoHitP->Fill(-pPart->P()); pFlux->Fill(-3.5,pHit->C());}break;
          case kPiMinus  : if(pPart->P()>cutPpionminus)  {pPioHitP->Fill(-pPart->P()); pFlux->Fill(-2.5,pHit->C());}break;
          case kMuonMinus: if(pPart->P()>cutPmuonminus)  {pMuoHitP->Fill(-pPart->P()); pFlux->Fill(-1.5,pHit->C());}break;        
          case kPositron : if(pPart->P()>cutPpositron)   {pEleHitP->Fill(-pPart->P()); pFlux->Fill(-0.5,pHit->C()); 
               pEleHitRP->Fill(-pPart->P(),pPart->R());  pEleHitRZ->Fill(pPart->Vz(),pPart->R()); }break;
          
          case kElectron : if(pPart->P()>cutPelectron)   {pEleHitP->Fill( pPart->P()); pFlux->Fill( 0.5,pHit->C()); 
               pEleHitRP->Fill( pPart->P(),pPart->R());  pEleHitRZ->Fill(pPart->Vz(),pPart->R()); }break;
          case kMuonPlus : if(pPart->P()>cutPmuonplus)   {pMuoHitP->Fill( pPart->P()); pFlux->Fill( 1.5,pHit->C());}break;                     
          case kPiPlus   : if(pPart->P()>cutPpionplus)   {pPioHitP->Fill( pPart->P()); pFlux->Fill( 2.5,pHit->C());}break;           
          case kKPlus    : if(pPart->P()>cutPkaonplus)   {pKaoHitP->Fill( pPart->P()); pFlux->Fill( 3.5,pHit->C());}break;           
          case kProton   : if(pPart->P()>cutPproton)     {pProHitP->Fill( pPart->P()); pFlux->Fill( 4.5,pHit->C());}break;
        }
      }//hits loop      
    }//TreeH loop
    iCntPrimParts +=pStack->GetNprimary();
    iCntTotParts  +=pStack->GetNtrack();
  }//events loop                        
//unload all loaded staff  
  fLoader->UnloadHits();  
  fLoader->GetRunLoader()->UnloadHeader(); 
  fLoader->GetRunLoader()->UnloadKinematics();  
//Calculater some sums
  Stat_t sum=0;
//sum row, sum over rows  
  for(Int_t i=1;i<=pFlux->GetNbinsX();i++){
    sum=0; for(Int_t j=2;j<=8;j++)    sum+=pFlux->GetBinContent(i,j);    
    pFlux->SetBinContent(i,1,sum);
  }
    
//display everything  
  new TCanvas("canvas1",Form("Events %i Nprims=%i Nparticles=%i",iNevents,iCntPrimParts,iCntTotParts),1000,900); pFlux->Draw("text");  gPad->SetGrid();  
//total prims and particles
  TLatex latex; latex.SetTextSize(0.02);
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,10);    latex.DrawLatex(5.1,9.5,Form("%.0f",sum));
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,9);     latex.DrawLatex(5.1,8.5,Form("%.0f",sum));
  for(Int_t iChN=1;iChN<=kNchambers;iChN++) {
    sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,iChN+1);latex.DrawLatex(5.1,iChN+0.5,Form("%.0f",sum));
  }  
  sum=0; for(Int_t i=1;i<=pFlux->GetNbinsX();i++) sum+=pFlux->GetBinContent(i,1);    latex.DrawLatex(5.1,0.5,Form("%.0f",sum));
  
  new TCanvas("cEleAllP"   ,"e" ,200,100); pEleAllP->Draw();
  new TCanvas("cEleHitRP"  ,"e" ,200,100); pEleHitRP->Draw();
  new TCanvas("cEleHitRZ"  ,"e" ,200,100); pEleHitRZ->Draw();
  new TCanvas("cEleHitP"   ,"e" ,200,100); pEleHitP->Draw();
  new TCanvas("cMuoHitP"   ,"mu",200,100); pMuoHitP->Draw();
  new TCanvas("cPioHitP"   ,"pi",200,100); pPioHitP->Draw();
  new TCanvas("cKaoHitP"   ,"K" ,200,100); pKaoHitP->Draw();
  new TCanvas("cProHitP"   ,"p" ,200,100); pProHitP->Draw();
  new TCanvas("cVertex"    ,"2d vertex" ,200,100); pVertex->Draw();
  new TCanvas("cRho"    ,"Rho of sec" ,200,100); pRho->Draw();
  
  gBenchmark->Show("HitsPlots");
}//HitsPlots()
