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
#include <TParticle.h>
#include <TStyle.h>
#include <TH2.h>
#include <TMath.h>

ClassImp(AliRICHDisplFast)

//__________________________________________________________________________________________________
void AliRICHDisplFast::Exec()
{
  TH2F *pHitsH2     = new TH2F("pHitsH2"  ,  "Event Display",165,-AliRICHParam::PcSizeY()/2,AliRICHParam::PcSizeY()/2,
                                                             144,-AliRICHParam::PcSizeX()/2,AliRICHParam::PcSizeX()/2);
  pHitsH2->SetStats(kFALSE);
  TH2F *pDigitsH2   = new TH2F("pDigitsH2"  ,"Event Display",165,-AliRICHParam::PcSizeY()/2,AliRICHParam::PcSizeY()/2,
                                                             144,-AliRICHParam::PcSizeX()/2,AliRICHParam::PcSizeX()/2);
  TH2F *pClustersH2 = new TH2F("pClustersH2","Event Display",165,-AliRICHParam::PcSizeY()/2,AliRICHParam::PcSizeY()/2,
                                                             144,-AliRICHParam::PcSizeX()/2,AliRICHParam::PcSizeX()/2);

  TCanvas *Display = new TCanvas("Display","RICH Display",0,0,600,600);
    
  gStyle->SetPalette(1);

  AliRICH *pRich = (AliRICH*)gAlice->GetDetector("RICH");
  pRich->GetLoader()->LoadRecPoints();
  pRich->GetLoader()->LoadDigits();
  pRich->GetLoader()->LoadHits();
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events Loop
    pRich->GetLoader()->GetRunLoader()->GetEvent(iEventN);
    pRich->GetLoader()->TreeD()->GetEntry(0);
    pRich->GetLoader()->TreeR()->GetEntry(0);

  
    Int_t nPrimaries = (Int_t)pRich->GetLoader()->TreeH()->GetEntries();
    TObjArray * Hits = new TObjArray[nPrimaries];
  
    for(Int_t i=0;i<nPrimaries;i++) {
      pRich->GetLoader()->TreeH()->GetEntry(i);
      Int_t nHits = pRich->Hits()->GetEntries();
      for(Int_t k=0;k<nHits;k++)         Hits[i].Add(pRich->Hits()->At(k));
    
    }
           
    for(Int_t iChamber=1;iChamber<=7;iChamber++){//modules loop
    
     Int_t nDigits   = pRich->Digits(iChamber)->GetEntries();
     Int_t nClusters = pRich->Clusters(iChamber)->GetEntries();
   
     pHitsH2->Reset();     pDigitsH2->Reset();     pClustersH2->Reset();

     Double_t xpad,ypad;

      for(Int_t i=0;i<nPrimaries;i++){//prims loop
        pRich->GetLoader()->TreeH()->GetEntry(i);
        Int_t nHits = pRich->Hits()->GetEntries();
        for(Int_t j=0;j<nHits;j++){//hits loop
          AliRICHhit *pHit = (AliRICHhit*)Hits[i].At(j);
          if(pHit->C()==iChamber){
            TVector3 xyzhit(pHit->X(),pHit->Y(),pHit->Z());
            TVector3 hitlocal = pRich->C(iChamber)->Glob2Loc(xyzhit);
            pHitsH2->Fill(hitlocal.X(),hitlocal.Y(),200);
          }//if
        }//hits loop         
      }//prims loop
     
      for(Int_t j=0;j<nDigits;j++){//digits loop
        AliRICHdigit *pDigit = (AliRICHdigit*)pRich->Digits(iChamber)->At(j);
	AliRICHParam::Pad2Loc(pDigit->X(),pDigit->Y(),xpad,ypad);
	pDigitsH2->Fill(xpad,ypad,100);
      }//digits loop
        
      for(Int_t j=0;j<nClusters;j++){//clusters loop
        AliRICHcluster *pCluster = (AliRICHcluster*)pRich->Clusters(iChamber)->At(j);
        pClustersH2->Fill(pCluster->X(),pCluster->Y(),50);
      }//clusters loop

      pHitsH2->SetTitle(Form("event %i module %2i",iEventN,iChamber));
      pHitsH2->SetMarkerColor(kRed); pHitsH2->SetMarkerStyle(29); pHitsH2->SetMarkerSize(0.4);
      pHitsH2->Draw();
      Display->Update();
      Display->Modified();
      getchar();
             
      pDigitsH2->SetMarkerColor(kGreen); pDigitsH2->SetMarkerStyle(29); pDigitsH2->SetMarkerSize(0.4);
      pDigitsH2->Draw("same");
      Display->Update();
      Display->Modified();       
      getchar();
      
      pClustersH2->SetMarkerColor(kBlue); pClustersH2->SetMarkerStyle(29);  pClustersH2->SetMarkerSize(0.4);
      pClustersH2->Draw("same");
      Display->Update();
      Display->Modified();
      getchar();
     }//modules loop
    delete [] Hits;
  }////events Loop
  pRich->GetLoader()->UnloadRecPoints();
  pRich->GetLoader()->UnloadDigits();
  pRich->GetLoader()->UnloadHits();
}//Exec()
//__________________________________________________________________________________________________
