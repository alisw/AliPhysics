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

#include "AliRICH.h"
#include "AliRICHDisplFast.h"
#include "AliRICHParam.h"
#include <AliLoader.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TText.h>

#define NPointsOfRing 201

// Geometry of the RICH at Star...

static const Int_t nPadX      = AliRICHParam::NpadsY();
static const Int_t nPadY      = AliRICHParam::NpadsX();
static const Float_t PadSizeX = AliRICHParam::PadSizeY();
static const Float_t PadSizeY = AliRICHParam::PadSizeX();
static const Float_t spacer   = AliRICHParam::DeadZone();
static const Float_t degree = 180/3.1415926535;

static const Float_t pi = TMath::Pi();

static const Float_t RadiatorWidth = AliRICHParam::FreonThickness();
static const Float_t QuartzWidth   = AliRICHParam::QuartzThickness();
static const Float_t GapWidth      = AliRICHParam::RadiatorToPads();

static const Float_t Xmin = -AliRICHParam::PcSizeY()/2.;
static const Float_t Xmax =  AliRICHParam::PcSizeY()/2.;
static const Float_t Ymin = -AliRICHParam::PcSizeX()/2.;
static const Float_t Ymax =  AliRICHParam::PcSizeX()/2.;


AliRICHDisplFast::AliRICHDisplFast()
{
}

void AliRICHDisplFast::Display()
{

  TText *text;

  TH2F *h2_disp   = new TH2F("h2_disp"  ,"STAR-RICH Event Display",165,Xmin,Xmax,100,Ymin,Ymax);
  TH2F *gHitsH2   = new TH2F("gHitsH2"  ,"STAR-RICH Event Display",165,Xmin,Xmax,100,Ymin,Ymax);
  TH2F *h2_dispad = new TH2F("h2_dispad","STAR-RICH Event Display",161,0.,161.,145,0.,145.);

  TCanvas *Display = new TCanvas("Display","Star Display",0,0,800,800);
    
  // Display event...

  gStyle->SetPalette(1,0);
    
  AliRICH *pRich = (AliRICH*)gAlice->GetDetector("RICH");
  pRich->GetLoader()->LoadDigits();
  pRich->GetLoader()->TreeD()->GetEntry(0);
  pRich->GetLoader()->LoadHits();

  Int_t nPrimaries = (Int_t)pRich->GetLoader()->TreeH()->GetEntries();
  TObjArray Hits[nPrimaries];
  
  for(Int_t i=0;i<nPrimaries;i++) {
    pRich->GetLoader()->TreeH()->GetEntry(i);
    Int_t nHits = pRich->Hits()->GetEntries();
    for(Int_t k=0;k<nHits;k++) {
       Hits[i].Add(pRich->Hits()->At(k));
    }
  }
           
  for(Int_t iChamber=1;iChamber<=7;iChamber++) {
    
     Int_t nDigits = pRich->DigitsOld(iChamber)->GetEntries();
   
     Display->ToggleEventStatus();
     Display->Modified();
     
     text = new TText(0,0,"");
     text->SetTextFont(61);
     text->SetTextSize(0.03);
     text->SetTextAlign(22);                                                       
      
//     Display->Resize();

     h2_disp->Reset();
     gHitsH2->Reset();
//     h2_dispad->Reset();

     Float_t xpad,ypad;
     
// loop for Digits...           
     
     for(Int_t j=0;j<nDigits;j++)
	{
          AliRICHDigit *pDigit = (AliRICHDigit*)pRich->DigitsOld(iChamber)->At(j);
	  AliRICHParam::Pad2Local(pDigit->PadX(),pDigit->PadY(),xpad,ypad);
          Float_t charge = (Float_t)pDigit->Signal();
	  h2_disp->Fill(xpad,ypad,charge);
          h2_dispad->Fill(pDigit->PadX(),pDigit->PadY(),pDigit->Signal());
          pDigit->Print();
	}

      cout << " -----------------" << endl;
      
      for(Int_t i=0;i<nPrimaries;i++) {
        pRich->GetLoader()->TreeH()->GetEntry(i);
        Int_t nHits = pRich->Hits()->GetEntries();

        for(Int_t j=0;j<nHits;j++) {
        
          AliRICHhit *pHit = (AliRICHhit*)Hits[i].At(j);
          cout << " chamber " << iChamber << " hit n. " << j << " ch " << pHit->Chamber() << endl;
          if(pHit->C()==iChamber) {
            TVector3 xyzhit(pHit->X(),pHit->Y(),pHit->Z());
            TVector3 hitlocal = pRich->C(iChamber)->Global2Local(xyzhit);
            gHitsH2->Fill(hitlocal.X(),hitlocal.Y());
//            pHit->Print();
          }
         }         
       }
             
      // end loop
                    

      h2_disp->SetMaximum(200);
      h2_disp->SetStats(0);
      h2_disp->Draw("colz");
      h2_disp->SetTitle(Form("Chamber %2i",iChamber));

/*      
      h2_dispad->SetTitle(Form("Chamber %2i",iChamber));
      h2_dispad->SetMaximum(200);
      h2_dispad->SetStats(0);
      h2_dispad->Draw("colz");
*/
      Display->Update();
      Display->Modified();
     
      getchar();
      gHitsH2->SetTitle(Form("Chamber %2i",iChamber));
      gHitsH2->SetMarkerColor(kBlue);
      gHitsH2->SetMarkerStyle(29);
      gHitsH2->SetMarkerSize(0.4);
      gHitsH2->Draw("same");
//      gHitsH2->Draw();
      Display->Update();
      Display->Modified();
       
      getchar();
    }
}

