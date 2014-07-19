/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliAnalysisTaskV2AllChAOD class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskV2AllChAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliPIDCombined.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include <TMCProcess.h>

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskV2AllChAOD)

//________________________________________________________________________
AliAnalysisTaskV2AllChAOD::AliAnalysisTaskV2AllChAOD(const char *name) : AliAnalysisTaskSE(name),  
  fAOD(0x0),
  fTrackCuts(0x0),
  fEventCuts(0x0),
  fIsMC(0),
  fCharge(0),
  fVZEROside(0),
  fOutput(0x0),
  fOutput_lq(0x0),
  fOutput_sq(0x0),
  fnCentBins(20),
  fnQvecBins(100),
  fIsQvecCalibMode(0),
  fQvecUpperLim(100),
  fCutLargeQperc(90.),
  fCutSmallQperc(10.),
  fEtaGapMin(-0.5),
  fEtaGapMax(0.5),
  fTrkBit(272),
  fEtaCut(0.8),
  fMinPt(0),
  fMaxPt(20.0),
  fMinTPCNcls(70),
  fResSP(0),
  fEta_vs_Phi_bef(0),
  fEta_vs_Phi(0),
  fResSP_lq(0),
  fResSP_sq(0)
{
  
  for (Int_t i = 0; i< 9; i++){
    fv2SPGap1A[i] = 0;
    fv2SPGap1B[i] = 0;

    fSinGap1Aq[i] = 0;
    fCosGap1Aq[i] = 0;
    fSinGap1Bq[i] = 0;
    fCosGap1Bq[i] = 0;

    fSinGap1A[i] = 0;
    fCosGap1A[i] = 0;
    fSinGap1B[i] = 0;
    fCosGap1B[i] = 0;
    
    //large q
    fv2SPGap1A_lq[i] = 0;
    fv2SPGap1B_lq[i] = 0;
    
    fSinGap1Aq_lq[i] = 0;
    fCosGap1Aq_lq[i] = 0;
    fSinGap1Bq_lq[i] = 0;
    fCosGap1Bq_lq[i] = 0;
    
    fSinGap1A_lq[i] = 0;
    fCosGap1A_lq[i] = 0;
    fSinGap1B_lq[i] = 0;
    fCosGap1B_lq[i] = 0;
    
    //small q
    fv2SPGap1A_sq[i] = 0;
    fv2SPGap1B_sq[i] = 0;
    
    fSinGap1Aq_sq[i] = 0;
    fCosGap1Aq_sq[i] = 0;
    fSinGap1Bq_sq[i] = 0;
    fCosGap1Bq_sq[i] = 0;
    
    fSinGap1A_sq[i] = 0;
    fCosGap1A_sq[i] = 0;
    fSinGap1B_sq[i] = 0;
    fCosGap1B_sq[i] = 0;
    
  }
    
  // Default constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskV2AllChAOD::UserCreateOutputObjects()
{
  // create output objects
  fOutput=new TList();
  fOutput->SetOwner();
  fOutput->SetName("fOutput");
  
  fOutput_lq=new TList();
  fOutput_lq->SetOwner();
  fOutput_lq->SetName("fOutput_lq");
  
  fOutput_sq=new TList();
  fOutput_sq->SetOwner();
  fOutput_sq->SetName("fOutput_sq");
  
  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  
  // binning common to all the THn
  //change it according to your needs + move it to global variables -> setter/getter
//   Double_t ptBins[] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
//   const Int_t nptBins = 31;
  Double_t ptBins[] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  const Int_t nptBins = 33;
  
  fResSP = new TProfile("fResSP", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput->Add(fResSP);
  
  fEta_vs_Phi_bef = new TH2F("fEta_vs_Phi_bef","eta vs phi distribution before eta gap;#eta;#phi",200.,-1.,1.,350.,0.,7.);
  fOutput->Add(fEta_vs_Phi_bef);
  
  fEta_vs_Phi = new TH2F("fEta_vs_Phi","eta vs phi distribution;#eta;#phi",200.,-1.,1.,350.,0.,7.);
  fOutput->Add(fEta_vs_Phi);
  
  //large q resolution
  fResSP_lq = new TProfile("fResSP_lq", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput_lq->Add(fResSP_lq);
  
  //small q resolution
  fResSP_sq = new TProfile("fResSP_sq", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput_sq->Add(fResSP_sq);
  
  for (Int_t iC = 0; iC < 9; iC++){

    fv2SPGap1A[iC] = new TProfile(Form("fv2SPGap1A_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1A[iC]);

    fv2SPGap1B[iC] = new TProfile(Form("fv2SPGap1B_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1B[iC]);

    fSinGap1Aq[iC] = new TProfile(Form("fSinGap1Aq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1Aq[iC]);
      
    fCosGap1Aq[iC] = new TProfile(Form("fCosGap1Aq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1Aq[iC]);
    
    fSinGap1Bq[iC] = new TProfile(Form("fSinGap1Bq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1Bq[iC]);
    
    fCosGap1Bq[iC] = new TProfile(Form("fCosGap1Bq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1Bq[iC]);

    fSinGap1A[iC] = new TProfile(Form("fSinGap1A_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1A[iC]);
      
    fCosGap1A[iC] = new TProfile(Form("fCosGap1A_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1A[iC]);
    
    fSinGap1B[iC] = new TProfile(Form("fSinGap1B_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1B[iC]);
    
    fCosGap1B[iC] = new TProfile(Form("fCosGap1B_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1B[iC]);
    
    //large q
    fv2SPGap1A_lq[iC] = new TProfile(Form("fv2SPGap1A_lq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_lq->Add(fv2SPGap1A_lq[iC]);

    fv2SPGap1B_lq[iC] = new TProfile(Form("fv2SPGap1B_lq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_lq->Add(fv2SPGap1B_lq[iC]);

    fSinGap1Aq_lq[iC] = new TProfile(Form("fSinGap1Aq_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1Aq_lq[iC]);
      
    fCosGap1Aq_lq[iC] = new TProfile(Form("fCosGap1Aq_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1Aq_lq[iC]);
    
    fSinGap1Bq_lq[iC] = new TProfile(Form("fSinGap1Bq_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1Bq_lq[iC]);
    
    fCosGap1Bq_lq[iC] = new TProfile(Form("fCosGap1Bq_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1Bq_lq[iC]);

    fSinGap1A_lq[iC] = new TProfile(Form("fSinGap1A_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1A_lq[iC]);
      
    fCosGap1A_lq[iC] = new TProfile(Form("fCosGap1A_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1A_lq[iC]);
    
    fSinGap1B_lq[iC] = new TProfile(Form("fSinGap1B_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1B_lq[iC]);
    
    fCosGap1B_lq[iC] = new TProfile(Form("fCosGap1B_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1B_lq[iC]);
    
    //small q
    fv2SPGap1A_sq[iC] = new TProfile(Form("fv2SPGap1A_sq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_sq->Add(fv2SPGap1A_sq[iC]);

    fv2SPGap1B_sq[iC] = new TProfile(Form("fv2SPGap1B_sq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_sq->Add(fv2SPGap1B_sq[iC]);

    fSinGap1Aq_sq[iC] = new TProfile(Form("fSinGap1Aq_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1Aq_sq[iC]);
      
    fCosGap1Aq_sq[iC] = new TProfile(Form("fCosGap1Aq_sq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1Aq_sq[iC]);
    
    fSinGap1Bq_sq[iC] = new TProfile(Form("fSinGap1Bq_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1Bq_sq[iC]);
    
    fCosGap1Bq_sq[iC] = new TProfile(Form("fCosGap1Bq_sq_%d", iC), "p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1Bq_sq[iC]);

    fSinGap1A_sq[iC] = new TProfile(Form("fSinGap1A_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1A_sq[iC]);
      
    fCosGap1A_sq[iC] = new TProfile(Form("fCosGap1A_sq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1A_sq[iC]);
    
    fSinGap1B_sq[iC] = new TProfile(Form("fSinGap1B_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1B_sq[iC]);
    
    fCosGap1B_sq[iC] = new TProfile(Form("fCosGap1B_sq_%d", iC), "p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1B_sq[iC]);
  };
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fOutput_lq  );
  PostData(5, fOutput_sq  );
}

//________________________________________________________________________

void AliAnalysisTaskV2AllChAOD::UserExec(Option_t *)
{
  //Printf("An event");
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
  
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }
  
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection

  Double_t Qvec=0.;//in case of MC we save space in the memory
  if(!fIsMC){
    if(fIsQvecCalibMode){
      if(fVZEROside==0)Qvec=fEventCuts->GetqV0A();
      else if (fVZEROside==1)Qvec=fEventCuts->GetqV0C();
    }
    else Qvec=fEventCuts->GetQvecPercentile(fVZEROside);
  }
  
  Double_t Cent=fEventCuts->GetCent();
  
  Short_t centV0 = -1;
  if ((Cent > 0) && (Cent <= 5.0))
      centV0 = 0; 
    else if ((Cent > 5.0) && (Cent <= 10.0))
      centV0 = 1;
    else if ((Cent > 10.0) && (Cent <= 20.0))
      centV0 = 2;
    else if ((Cent > 20.0) && (Cent <= 30.0))
      centV0 = 3;   
    else if ((Cent > 30.0) && (Cent <= 40.0))
      centV0 = 4; 
    else if ((Cent > 40.0) && (Cent <= 50.0))
      centV0 = 5;  
    else if ((Cent > 50.0) && (Cent <= 60.0))
      centV0 = 6;
    else if ((Cent > 60.0) && (Cent <= 70.0))
      centV0 = 7;         			       
    else if ((Cent > 70.0) && (Cent <= 80.0))
      centV0 = 8; 
    

  Double_t QxGap1A = 0., QyGap1A = 0.;
  Double_t QxGap1B = 0., QyGap1B = 0.;
  Int_t multGap1A = 0, multGap1B = 0;
  
  for (Int_t loop = 0; loop < 2; loop++){

    //main loop on tracks
    for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
      AliAODTrack* track = fAOD->GetTrack(iTracks);
      if(fCharge != 0 && track->Charge() != fCharge) continue;//if fCharge != 0 only select fCharge 
      if (!fTrackCuts->IsSelected(track,kTRUE)) continue; //track selection (rapidity selection NOT in the standard cuts)
    
      fEta_vs_Phi_bef->Fill( track->Eta(), track->Phi() );
  
      if (loop == 0) {
	
	if (track->Eta() > fEtaGapMax){
	  QxGap1A += TMath::Cos(2.*track->Phi());
          QyGap1A += TMath::Sin(2.*track->Phi());
          multGap1A++;

          fSinGap1Aq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
          fCosGap1Aq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  
	  fEta_vs_Phi->Fill( track->Eta(), track->Phi() );
                    
          if (Qvec > fCutLargeQperc){
	    fSinGap1Aq_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1Aq_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
      
          if (Qvec < fCutSmallQperc){
	    fSinGap1Aq_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1Aq_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
	}
    
      if (track->Eta() < fEtaGapMin){
	QxGap1B += TMath::Cos(2.*track->Phi());
        QyGap1B += TMath::Sin(2.*track->Phi());
        multGap1B++;
                    
        fCosGap1Bq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
        fSinGap1Bq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	
	fEta_vs_Phi->Fill( track->Eta(), track->Phi() );
                    
        if (Qvec > fCutLargeQperc){
	  fSinGap1Bq_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	  fCosGap1Bq_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
        }
      
        if (Qvec < fCutSmallQperc){
	  fSinGap1Bq_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	  fCosGap1Bq_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
        }
      }
  
    } else {
      
        //eval v2 scalar product
        if (track->Eta() < fEtaGapMin && multGap1A > 0){
          Double_t v2SPGap1A = (TMath::Cos(2.*track->Phi())*QxGap1A + TMath::Sin(2.*track->Phi())*QyGap1A)/(Double_t)multGap1A;
          fv2SPGap1A[centV0]->Fill(track->Pt(), v2SPGap1A);

	  fSinGap1A[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
          fCosGap1A[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
      
          if (Qvec > fCutLargeQperc){
	    fv2SPGap1A_lq[centV0]->Fill(track->Pt(), v2SPGap1A);
	    fSinGap1A_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1A_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
      
          if (Qvec < fCutSmallQperc){
	    fv2SPGap1A_sq[centV0]->Fill(track->Pt(), v2SPGap1A);
	    fSinGap1A_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1A_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
        }
      
        if (track->Eta() > fEtaGapMax && multGap1B > 0){
          Double_t v2SPGap1B = (TMath::Cos(2.*track->Phi())*QxGap1B + TMath::Sin(2.*track->Phi())*QyGap1B)/(Double_t)multGap1B;
          fv2SPGap1B[centV0]->Fill(track->Pt(), v2SPGap1B);
	  
          fCosGap1B[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          fSinGap1B[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
      
          if (Qvec > fCutLargeQperc){
	    fv2SPGap1B_lq[centV0]->Fill(track->Pt(), v2SPGap1B);
	    fSinGap1B_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1B_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
      
          if (Qvec < fCutSmallQperc){
	    fv2SPGap1B_sq[centV0]->Fill(track->Pt(), v2SPGap1B);
	    fSinGap1B_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1B_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
        }
      }// end else 
    } // end loop on tracks
  } // end loop
  
  
  if (multGap1A > 0 && multGap1B > 0){
    Double_t res = (QxGap1A*QxGap1B + QyGap1A*QyGap1B)/(Double_t)multGap1A/(Double_t)multGap1B;
    fResSP->Fill((Double_t)centV0, res);
        
    if (Qvec > fCutLargeQperc)
      fResSP_lq->Fill((Double_t)centV0, res);

    if (Qvec < fCutSmallQperc)
      fResSP_sq->Fill((Double_t)centV0, res);
  }
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fOutput_lq  );
  PostData(5, fOutput_sq  );
}

//_________________________________________________________________
Bool_t  AliAnalysisTaskV2AllChAOD::GetDCA(const AliAODTrack* trk, Double_t * p){
  
  //AliAODTrack::DCA(): for newest AOD fTrack->DCA() always gives -999. This should fix.
  //FIXME should update EventCuts?
  //FIXME add track->GetXYZ(p) method
  
  double xyz[3],cov[3];
  
  if (!trk->GetXYZ(xyz)) { // dca is not stored
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(trk);
    AliVEvent* ev = (AliVEvent*)trk->GetEvent();
    if (!ev) {/*printf("Event is not connected to the track\n");*/ return kFALSE;}
    if (!etp.PropagateToDCA(ev->GetPrimaryVertex(), ev->GetMagneticField(),999,xyz,cov)) return kFALSE; // failed, track is too far from vertex
  }
  p[0] = xyz[0];
  p[1] = xyz[1];
  return kTRUE;

}

//_________________________________________________________________
void   AliAnalysisTaskV2AllChAOD::Terminate(Option_t *)
{
  // Terminate
}
