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
#include <TCanvas.h>
#include <TPolyLine.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TH2.h>
#include <TMath.h>
#include <TLatex.h>

ClassImp(AliRICHDisplFast)

//__________________________________________________________________________________________________
void AliRICHDisplFast::Exec()
{
  AliRICH *pRich = (AliRICH*)gAlice->GetDetector("RICH");
  Bool_t isHits    =!pRich->GetLoader()->LoadHits();
  Bool_t isDigits  =!pRich->GetLoader()->LoadDigits();
  Bool_t isClusters=!pRich->GetLoader()->LoadRecPoints();
  
  if(!isHits && !isDigits && !isClusters){Error("Exec","No hits digits and clusters. Nothing to display.");return;}
  
  TCanvas *Display = new TCanvas("Display","RICH Display",0,0,600,600);
  gStyle->SetPalette(1);
  
  TH2F *pHitsH2=0,*pDigitsH2=0,*pClustersH2=0;
  
  if(isHits)     pHitsH2     = new TH2F("pHitsH2"  ,  "Event Display",165,-AliRICHParam::PcSizeX()/2,AliRICHParam::PcSizeX()/2,
                                                             144,-AliRICHParam::PcSizeY()/2,AliRICHParam::PcSizeY()/2);
  if(pHitsH2)    pHitsH2->SetStats(kFALSE);
  
  if(isDigits)   pDigitsH2   = new TH2F("pDigitsH2"  ,"Event Display",165,-AliRICHParam::PcSizeX()/2,AliRICHParam::PcSizeX()/2,
                                                             144,-AliRICHParam::PcSizeY()/2,AliRICHParam::PcSizeY()/2);
  if(isClusters) pClustersH2 = new TH2F("pClustersH2","Event Display",165,-AliRICHParam::PcSizeX()/2,AliRICHParam::PcSizeX()/2,
                                                             144,-AliRICHParam::PcSizeY()/2,AliRICHParam::PcSizeY()/2);

    

  
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
            TVector2 hitLocX2 = pRich->C(iChamber)->Glob2Loc(hitGlobX3);
            pHitsH2->Fill(hitLocX2.X(),hitLocX2.Y(),200);
          }//if
        }//hits loop         
      }//prims loop
      pHitsH2->SetTitle(Form("event %i chamber %2i",iEventN,iChamber));
      pHitsH2->SetMarkerColor(kRed); pHitsH2->SetMarkerStyle(29); pHitsH2->SetMarkerSize(0.4);
      pHitsH2->Draw();
      DrawSectors();
      TLatex l; l.SetNDC(); l.SetTextSize(0.02);
      if(!isHits)     {l.SetTextColor(kRed)  ;l.DrawLatex(0.1,0.01,"No Hits"    );}
      if(!isDigits)   {l.SetTextColor(kGreen);l.DrawLatex(0.4,0.01,"No DIGITS"  );}
      if(!isClusters) {l.SetTextColor(kBlue) ;l.DrawLatex(0.8,0.01,"No CLUSTERS");}
      Display->Update();
      Display->Modified();
      getchar();             
//deals with digits      
      if(isDigits){
        pRich->GetLoader()->TreeD()->GetEntry(0);
        for(Int_t j=0;j<pRich->Digits(iChamber)->GetEntries();j++){//digits loop
          AliRICHdigit *pDig = (AliRICHdigit*)pRich->Digits(iChamber)->At(j);
	  TVector2 x2=AliRICHParam::Pad2Loc(pDig->X(),pDig->Y());
	  pDigitsH2->Fill(x2.X(),x2.Y(),100);
        }//digits loop
        pDigitsH2->SetMarkerColor(kGreen); pDigitsH2->SetMarkerStyle(29); pDigitsH2->SetMarkerSize(0.4);
        pDigitsH2->Draw("same");
        Display->Update();
        Display->Modified();       
        getchar();
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
        getchar();
      }//if(isClusters)
    }//chambers loop
    delete [] Hits;
  }//events Loop
  pRich->GetLoader()->UnloadHits();
  if(isDigits)   pRich->GetLoader()->UnloadDigits();
  if(isClusters) pRich->GetLoader()->UnloadRecPoints();
}//Exec()
//__________________________________________________________________________________________________
void AliRICHDisplFast::DrawSectors() 
{ 
  Double_t x1[5] = {-AliRICHParam::PcSizeX()/2,
                    -AliRICHParam::SectorSizeX()/2-AliRICHParam::DeadZone(),
                    -AliRICHParam::SectorSizeX()/2-AliRICHParam::DeadZone(),
                    -AliRICHParam::PcSizeX()/2,
                    -AliRICHParam::PcSizeX()/2};
  Double_t y1[5] = {AliRICHParam::DeadZone()/2,
                    AliRICHParam::DeadZone()/2,
                    AliRICHParam::PcSizeY()/2,
                    AliRICHParam::PcSizeY()/2,
                    AliRICHParam::DeadZone()/2};  
  Double_t x2[5] = {-AliRICHParam::SectorSizeX()/2,
                     AliRICHParam::SectorSizeX()/2,
                     AliRICHParam::SectorSizeX()/2,
                    -AliRICHParam::SectorSizeX()/2,
                    -AliRICHParam::SectorSizeX()/2};
  Double_t y2[5] = {AliRICHParam::DeadZone()/2,
                    AliRICHParam::DeadZone()/2,
                    AliRICHParam::PcSizeY()/2,
                    AliRICHParam::PcSizeY()/2,
                    AliRICHParam::DeadZone()/2};
  Double_t x3[5] = { AliRICHParam::SectorSizeX()/2+AliRICHParam::DeadZone(),
                     AliRICHParam::PcSizeX()/2,
                     AliRICHParam::PcSizeX()/2,
                     AliRICHParam::SectorSizeX()/2+AliRICHParam::DeadZone(),
                     AliRICHParam::SectorSizeX()/2+AliRICHParam::DeadZone()};
  Double_t y3[5] = {AliRICHParam::DeadZone()/2,
                    AliRICHParam::DeadZone()/2,
                    AliRICHParam::PcSizeY()/2,
                    AliRICHParam::PcSizeY()/2,
                    AliRICHParam::DeadZone()/2};
  Double_t x4[5] = {-AliRICHParam::PcSizeX()/2,
                    -AliRICHParam::SectorSizeX()/2-AliRICHParam::DeadZone(),
                    -AliRICHParam::SectorSizeX()/2-AliRICHParam::DeadZone(),
                    -AliRICHParam::PcSizeX()/2,
                    -AliRICHParam::PcSizeX()/2};
  Double_t y4[5] = {-AliRICHParam::PcSizeY()/2,
                    -AliRICHParam::PcSizeY()/2,
                    -AliRICHParam::DeadZone()/2,
                    -AliRICHParam::DeadZone()/2,
                    -AliRICHParam::PcSizeY()/2};
   Double_t x5[5] = {-AliRICHParam::SectorSizeX()/2,
                     AliRICHParam::SectorSizeX()/2,
                     AliRICHParam::SectorSizeX()/2,
                    -AliRICHParam::SectorSizeX()/2,
                    -AliRICHParam::SectorSizeX()/2};
  Double_t y5[5] = {-AliRICHParam::PcSizeY()/2,
                    -AliRICHParam::PcSizeY()/2,
                    -AliRICHParam::DeadZone()/2,
                    -AliRICHParam::DeadZone()/2,
                    -AliRICHParam::PcSizeY()/2};
  Double_t x6[5] = { AliRICHParam::SectorSizeX()/2+AliRICHParam::DeadZone(),
                     AliRICHParam::PcSizeX()/2,
                     AliRICHParam::PcSizeX()/2,
                     AliRICHParam::SectorSizeX()/2+AliRICHParam::DeadZone(),
                     AliRICHParam::SectorSizeX()/2+AliRICHParam::DeadZone()};
  Double_t y6[5] = {-AliRICHParam::PcSizeY()/2,
                    -AliRICHParam::PcSizeY()/2,
                    -AliRICHParam::DeadZone()/2,
                    -AliRICHParam::DeadZone()/2,
                    -AliRICHParam::PcSizeY()/2};
  TPolyLine *sector1 = new TPolyLine(5,x1,y1);  sector1->SetLineColor(21);  sector1->Draw();
  TPolyLine *sector2 = new TPolyLine(5,x2,y2);  sector2->SetLineColor(21);  sector2->Draw();
  TPolyLine *sector3 = new TPolyLine(5,x3,y3);  sector3->SetLineColor(21);  sector3->Draw();
  TPolyLine *sector4 = new TPolyLine(5,x4,y4);  sector4->SetLineColor(21);  sector4->Draw();
  TPolyLine *sector5 = new TPolyLine(5,x5,y5);  sector5->SetLineColor(21);  sector5->Draw();
  TPolyLine *sector6 = new TPolyLine(5,x6,y6);  sector6->SetLineColor(21);  sector6->Draw();
}//DrawSectors()
