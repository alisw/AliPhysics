/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliRICHDisplFast.h"
#include "AliRICH.h"
#include "AliRICHChamber.h"
#include "AliRICHParam.h"
#include <AliLoader.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH2F.h>
#include <TParticle.h>

ClassImp(AliRICHDisplFast)
//__________________________________________________________________________________________________
void AliRICHDisplFast::ShowEvent(Int_t iEvtNmin,Int_t iEvtNmax)
{
  TH2F *pDigitsH2[8];
  char titobj[11],titdisp[20];

  AliRICH *pRich = (AliRICH*)gAlice->GetDetector("RICH");
  Bool_t isDigits  =!pRich->GetLoader()->LoadDigits();
  if(!isDigits){Error("ShoEvent","No digits. Nothing to display.");return;}
  
  TCanvas *canvas = new TCanvas("RICHDisplay","RICH Display",0,0,1226,900);
  
  gStyle->SetPalette(1);

  canvas->Divide(3,3);
  
  for(Int_t iChamber=1;iChamber<=7;iChamber++) {
    sprintf(titobj,"pDigitsH2_%i",iChamber);
    sprintf(titdisp,"Chamber  %i",iChamber);
    pDigitsH2[iChamber] = new TH2F(titobj,titdisp,165,0,AliRICHParam::PcSizeX(),144,0,AliRICHParam::PcSizeY());
    pDigitsH2[iChamber]->SetMarkerColor(kGreen); 
    pDigitsH2[iChamber]->SetMarkerStyle(29); 
    pDigitsH2[iChamber]->SetMarkerSize(0.4);
    pDigitsH2[iChamber]->SetStats(kFALSE);
  }
  
  if(iEvtNmax>gAlice->GetEventsPerRun()) iEvtNmax=gAlice->GetEventsPerRun();

  TLatex t;
  t.SetTextSize(0.10);

  TText *tit = new TText(0.1,0.6,"RICH Display");
  tit->SetTextSize(0.10);

  for(Int_t iEventN=iEvtNmin;iEventN<=iEvtNmax;iEventN++) {
        
    canvas->cd(1);
    sprintf(titdisp,"Event Number %i",iEventN);
    t.SetText(0.2,0.4,titdisp);
    t.Draw();
    tit->Draw();
    pRich->GetLoader()->GetRunLoader()->GetEvent(iEventN);
    pRich->GetLoader()->TreeD()->GetEntry(0);
    for(Int_t iChamber=1;iChamber<=7;iChamber++) {
      pDigitsH2[iChamber]->Reset();    
      for(Int_t j=0;j<pRich->Digits(iChamber)->GetEntries();j++) {//digits loop
        AliRICHdigit *pDig = (AliRICHdigit*)pRich->Digits(iChamber)->At(j);
        TVector2 x2=AliRICHParam::Pad2Loc(pDig->Pad());
        pDigitsH2[iChamber]->Fill(x2.X(),x2.Y());
      }//digits loop
      if(iChamber==1) canvas->cd(7);
      if(iChamber==2) canvas->cd(8);
      if(iChamber==3) canvas->cd(4);
      if(iChamber==4) canvas->cd(5);
      if(iChamber==5) canvas->cd(6);
      if(iChamber==6) canvas->cd(2);
      if(iChamber==7) canvas->cd(3);
      pDigitsH2[iChamber]->Draw();
    }
    canvas->Update();
    canvas->Modified();
    
    if(iEvtNmin<iEvtNmax) gPad->WaitPrimitive();
  }
}//ShowEvent()
//__________________________________________________________________________________________________
void AliRICHDisplFast::Exec(Option_t *)
{
  AliRICH *pRich = (AliRICH*)gAlice->GetDetector("RICH");
  Bool_t isHits    =!pRich->GetLoader()->LoadHits();
  Bool_t isDigits  =!pRich->GetLoader()->LoadDigits();
  Bool_t isClusters=!pRich->GetLoader()->LoadRecPoints();
  
  if(!isHits && !isDigits && !isClusters){Error("Exec","No hits digits and clusters. Nothing to display.");return;}
  
  TCanvas *Display = new TCanvas("Display","RICH Display",0,0,600,600);
  gStyle->SetPalette(1);
  
  TH2F *pHitsH2=0,*pDigitsH2=0,*pClustersH2=0;
  
  if(isHits)     pHitsH2     = new TH2F("pHitsH2"  ,  "Event Display;x,cm;y,cm",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  if(pHitsH2)    pHitsH2->SetStats(kFALSE);
  
  if(isDigits)   pDigitsH2   = new TH2F("pDigitsH2"  ,"Event Display",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  if(isClusters) pClustersH2 = new TH2F("pClustersH2","Event Display",165,0,AliRICHParam::PcSizeX(),
                                                                      144,0,AliRICHParam::PcSizeY());
  
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events Loop
    pRich->GetLoader()->GetRunLoader()->GetEvent(iEventN);

  
    Int_t nPrimaries = (Int_t)pRich->GetLoader()->TreeH()->GetEntries();
    TObjArray * Hits = new TObjArray[nPrimaries];
  
    for(Int_t i=0;i<nPrimaries;i++) {
      pRich->GetLoader()->TreeH()->GetEntry(i);
      Int_t nHits = pRich->Hits()->GetEntries();
      for(Int_t k=0;k<nHits;k++)         Hits[i].Add(pRich->Hits()->At(k));
    
    }
//display all the staff on chamber by chamber basis           
    for(Int_t iChamber=1;iChamber<=7;iChamber++){//chambers loop       
      if(isHits)     pHitsH2    ->Reset();     
      if(isDigits)   pDigitsH2  ->Reset();     
      if(isClusters) pClustersH2->Reset();
//deals with hits
      for(Int_t i=0;i<nPrimaries;i++){//prims loop
        pRich->GetLoader()->TreeH()->GetEntry(i);
        Int_t nHits = pRich->Hits()->GetEntries();
        for(Int_t j=0;j<nHits;j++){//hits loop
          AliRICHhit *pHit = (AliRICHhit*)Hits[i].At(j);
          if(pHit->C()==iChamber){
            TVector3 hitGlobX3= pHit->OutX3();
            TVector2 hitLocX2 = pRich->C(iChamber)->Mrs2Pc(hitGlobX3);
            pHitsH2->Fill(hitLocX2.X(),hitLocX2.Y(),200);
          }//if
        }//hits loop         
      }//prims loop
      pHitsH2->SetTitle(Form("event %i chamber %2i",iEventN,iChamber));
      pHitsH2->SetMarkerColor(kRed); pHitsH2->SetMarkerStyle(29); pHitsH2->SetMarkerSize(0.4);
      pHitsH2->Draw();
      AliRICHParam::DrawSectors();
      TLatex l; l.SetNDC(); l.SetTextSize(0.02);
      if(!isHits)     {l.SetTextColor(kRed)  ;l.DrawLatex(0.1,0.01,"No Hits"    );}
      if(!isDigits)   {l.SetTextColor(kGreen);l.DrawLatex(0.4,0.01,"No DIGITS"  );}
      if(!isClusters) {l.SetTextColor(kBlue) ;l.DrawLatex(0.8,0.01,"No CLUSTERS");}
      Display->Update();
      Display->Modified();
      gPad->WaitPrimitive();
//deals with digits      
      if(isDigits){
        pRich->GetLoader()->TreeD()->GetEntry(0);
        for(Int_t j=0;j<pRich->Digits(iChamber)->GetEntries();j++){//digits loop
          AliRICHdigit *pDig = (AliRICHdigit*)pRich->Digits(iChamber)->At(j);
	  TVector2 x2=AliRICHParam::Pad2Loc(pDig->Pad());
	  pDigitsH2->Fill(x2.X(),x2.Y(),100);
        }//digits loop
        pDigitsH2->SetMarkerColor(kGreen); pDigitsH2->SetMarkerStyle(29); pDigitsH2->SetMarkerSize(0.4);
        pDigitsH2->Draw("same");
        Display->Update();
        Display->Modified();       
        gPad->WaitPrimitive();
      }//if(isDigits)      
//deals with clusters      
      if(isClusters){
        pRich->GetLoader()->TreeR()->GetEntry(0);
        for(Int_t j=0;j<pRich->Clusters(iChamber)->GetEntries();j++){//clusters loop
          AliRICHcluster *pClus = (AliRICHcluster*)pRich->Clusters(iChamber)->At(j);
          pClustersH2->Fill(pClus->X(),pClus->Y(),50);
        }//clusters loop
        pClustersH2->SetMarkerColor(kBlue); pClustersH2->SetMarkerStyle(29);  pClustersH2->SetMarkerSize(0.4);
        pClustersH2->Draw("same");
        Display->Update();
        Display->Modified();
        gPad->WaitPrimitive();
      }//if(isClusters)
    }//chambers loop
    delete [] Hits;
  }//events Loop
  pRich->GetLoader()->UnloadHits();
  if(isDigits)   pRich->GetLoader()->UnloadDigits();
  if(isClusters) pRich->GetLoader()->UnloadRecPoints();
}//Exec()
//__________________________________________________________________________________________________
Int_t AliRICHDisplFast::Nparticles(Int_t iPartID,Int_t iEvtN,AliRunLoader *pRL)
{
//counts total number of particles of given type (including secondary) for a given event
  pRL->GetEvent(iEvtN);    
  if(pRL->LoadHeader()) return 0;
  if(pRL->LoadKinematics()) return 0;
  AliStack *pStack=pRL->Stack();
  
  Int_t iCounter=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++){
    if(pStack->Particle(i)->GetPdgCode()==iPartID) iCounter++;
  }
  
  pRL->UnloadHeader();
  pRL->UnloadKinematics();
  return iCounter;
}
