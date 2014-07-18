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
  fQxGap1A(0),
  fQyGap1A(0),
  fmultGap1A(0),
  fQxGap1B(0),
  fQyGap1B(0),
  fmultGap1B(0),
  fResSP_lq(0),
  fQxGap1A_lq(0),
  fQyGap1A_lq(0),
  fmultGap1A_lq(0),
  fQxGap1B_lq(0),
  fQyGap1B_lq(0),
  fmultGap1B_lq(0),
  fResSP_sq(0),
  fQxGap1A_sq(0),
  fQyGap1A_sq(0),
  fmultGap1A_sq(0),
  fQxGap1B_sq(0),
  fQyGap1B_sq(0),
  fmultGap1B_sq(0)
{
  
  for (Int_t i = 0; i< 9; i++){
    fv2SPGap1A[i] = 0;
    fh2v2SPGap1A[i] = 0;
    fv2SPGap1B[i] = 0;
    fh2v2SPGap1B[i] = 0;

    fSinGap1A[i] = 0;
    fCosGap1A[i] = 0;
    fSinGap1B[i] = 0;
    fCosGap1B[i] = 0;
    
    //large q
    fv2SPGap1A_lq[i] = 0;
    fv2SPGap1B_lq[i] = 0;
    fSinGap1A_lq[i] = 0;
    fCosGap1A_lq[i] = 0;
    fSinGap1B_lq[i] = 0;
    fCosGap1B_lq[i] = 0;
    
    //small q
    fv2SPGap1A_sq[i] = 0;
    fv2SPGap1B_sq[i] = 0;
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
  
    fQxGap1A = new TProfile("fQxGap1A", "Qx mean; centrality; <Qx>", 9, -0.5, 8.5);
    fOutput->Add(fQxGap1A);
    fQyGap1A = new TProfile("fQyGap1A", "Qy mean; centrality; <Qy>", 9, -0.5, 8.5);
    fOutput->Add(fQyGap1A);
    fmultGap1A = new TProfile("fmultGap1A", " Multiplicity B; centrality; M_{A}", 9, -0.5, 8.5);
    fOutput->Add(fmultGap1A);
  
    fQxGap1B= new TProfile("fQxGap1B", "QxB mean; centrality; <Qx>", 9, -0.5, 8.5);
    fOutput->Add(fQxGap1B);
    fQyGap1B = new TProfile("fQyGap1B", "QyB mean; centrality; <Qy>", 9, -0.5, 8.5);
    fOutput->Add(fQyGap1B);
    fmultGap1B = new TProfile("fmultGap1B", " Multiplicity B; centrality; M_{B}", 9, -0.5, 8.5);
    fOutput->Add(fmultGap1B);
  
  //large q resolution
  fResSP_lq = new TProfile("fResSP_lq", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput_lq->Add(fResSP_lq);
  
    fQxGap1A_lq = new TProfile("fQxGap1A_lq", "Qx mean; centrality; <Qx>", 9, -0.5, 8.5);
    fOutput_lq->Add(fQxGap1A_lq);
    fQyGap1A_lq = new TProfile("fQyGap1A_lq", "Qy mean; centrality; <Qy>", 9, -0.5, 8.5);
    fOutput_lq->Add(fQyGap1A_lq);
    fmultGap1A_lq = new TProfile("fmultGap1A_lq", " Multiplicity B; centrality; M_{A}", 9, -0.5, 8.5);
    fOutput_lq->Add(fmultGap1A_lq);
  
    fQxGap1B_lq= new TProfile("fQxGap1B_lq", "QxB mean; centrality; <Qx>", 9, -0.5, 8.5);
    fOutput_lq->Add(fQxGap1B_lq);
    fQyGap1B_lq = new TProfile("fQyGap1B_lq", "QyB mean; centrality; <Qy>", 9, -0.5, 8.5);
    fOutput_lq->Add(fQyGap1B_lq);
    fmultGap1B_lq = new TProfile("fmultGap1B_lq", " Multiplicity B; centrality; M_{B}", 9, -0.5, 8.5);
    fOutput_lq->Add(fmultGap1B_lq);
  
  //small q resolution
  fResSP_sq = new TProfile("fResSP_sq", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput_sq->Add(fResSP_sq);
  
    fQxGap1A_sq = new TProfile("fQxGap1A_sq", "Qx mean; centrality; <Qx>", 9, -0.5, 8.5);
    fOutput_sq->Add(fQxGap1A_sq);
    fQyGap1A_sq = new TProfile("fQyGap1A_sq", "Qy mean; centrality; <Qy>", 9, -0.5, 8.5);
    fOutput_sq->Add(fQyGap1A_sq);
    fmultGap1A_sq = new TProfile("fmultGap1A_sq", " Multiplicity B; centrality; M_{A}", 9, -0.5, 8.5); 
    fOutput_sq->Add(fmultGap1A_sq);
  
    fQxGap1B_sq= new TProfile("fQxGap1B_sq", "QxB mean; centrality; <Qx>", 9, -0.5, 8.5);
    fOutput_sq->Add(fQxGap1B_sq);
    fQyGap1B_sq = new TProfile("fQyGap1B_sq", "QyB mean; centrality; <Qy>", 9, -0.5, 8.5);
    fOutput_sq->Add(fQyGap1B_sq);
    fmultGap1B_sq = new TProfile("fmultGap1B_sq", " Multiplicity B; centrality; M_{B}", 9, -0.5, 8.5);
    fOutput_sq->Add(fmultGap1B_sq);

  for (Int_t iC = 0; iC < 9; iC++){

    fv2SPGap1A[iC] = new TProfile(Form("fv2SPGap1A_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1A[iC]);

    fh2v2SPGap1A[iC] = new TH2F(Form("fh2v2SPGap1A_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins,200,-10,10);
    fOutput->Add(fh2v2SPGap1A[iC]);
 
    fv2SPGap1B[iC] = new TProfile(Form("fv2SPGap1B_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1B[iC]);

    fh2v2SPGap1B[iC] = new TH2F(Form("fh2v2SPGap1B_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins,200,-10,10);
    fOutput->Add(fh2v2SPGap1B[iC]);
      
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
    
  
      if (loop == 0) {
	
	if (track->Eta() > fEtaGapMax){
	  QxGap1A += TMath::Cos(2.*track->Phi());
          QyGap1A += TMath::Sin(2.*track->Phi());
          multGap1A++;

          fSinGap1A[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
          fCosGap1A[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
                    
          if (Qvec > fCutLargeQperc){
	    fSinGap1A_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1A_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
      
          if (Qvec < fCutSmallQperc){
	    fSinGap1A_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	    fCosGap1A_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
	  }
	}
    
      if (track->Eta() < fEtaGapMin){
	QxGap1B += TMath::Cos(2.*track->Phi());
        QyGap1B += TMath::Sin(2.*track->Phi());
        multGap1B++;
                    
        fCosGap1B[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
        fSinGap1B[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
                    
        if (Qvec > fCutLargeQperc){
	  fSinGap1B_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	  fCosGap1B_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
        }
      
        if (Qvec < fCutSmallQperc){
	  fSinGap1B_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
	  fCosGap1B_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
        }
      }
  
    } else {
      
        //eval v2 scalar product
        if (track->Eta() < fEtaGapMin && multGap1A > 0){
          Double_t v2SPGap1A = (TMath::Cos(2.*track->Phi())*QxGap1A + TMath::Sin(2.*track->Phi())*QyGap1A)/(Double_t)multGap1A;
          fv2SPGap1A[centV0]->Fill(track->Pt(), v2SPGap1A);
          fh2v2SPGap1A[centV0]->Fill(track->Pt(), v2SPGap1A);
      
          if (Qvec > fCutLargeQperc)
	    fv2SPGap1A_lq[centV0]->Fill(track->Pt(), v2SPGap1A);
      
         if (Qvec < fCutSmallQperc)
	    fv2SPGap1A_sq[centV0]->Fill(track->Pt(), v2SPGap1A);
        }
      
        if (track->Eta() > fEtaGapMax && multGap1B > 0){
          Double_t v2SPGap1B = (TMath::Cos(2.*track->Phi())*QxGap1B + TMath::Sin(2.*track->Phi())*QyGap1B)/(Double_t)multGap1B;
          fv2SPGap1B[centV0]->Fill(track->Pt(), v2SPGap1B);
          fh2v2SPGap1B[centV0]->Fill(track->Pt(), v2SPGap1B);
      
          if (Qvec > fCutLargeQperc)
	    fv2SPGap1B_lq[centV0]->Fill(track->Pt(), v2SPGap1B);
      
          if (Qvec < fCutSmallQperc)
	    fv2SPGap1B_sq[centV0]->Fill(track->Pt(), v2SPGap1B);
        }
      }// end else 
    } // end loop on tracks
  } // end loop
  
  
  if (multGap1A > 0 && multGap1B > 0){
    Double_t res = (QxGap1A*QxGap1B + QyGap1A*QyGap1B)/(Double_t)multGap1A/(Double_t)multGap1B;
    fResSP->Fill((Double_t)centV0, res);
      fQxGap1A->Fill((Double_t)centV0, QxGap1A);
      fQyGap1A->Fill((Double_t)centV0, QyGap1A);
      fmultGap1A->Fill((Double_t)centV0, multGap1A);
      fQxGap1B->Fill((Double_t)centV0, QxGap1B);
      fQyGap1B->Fill((Double_t)centV0, QyGap1B);
      fmultGap1B->Fill((Double_t)centV0, multGap1B);
        
    if (Qvec > fCutLargeQperc)
      fResSP_lq->Fill((Double_t)centV0, res);
        fQxGap1A_lq->Fill((Double_t)centV0, QxGap1A);
        fQyGap1A_lq->Fill((Double_t)centV0, QyGap1A);
        fmultGap1A_lq->Fill((Double_t)centV0, multGap1A);
        fQxGap1B_lq->Fill((Double_t)centV0, QxGap1B);
        fQyGap1B_lq->Fill((Double_t)centV0, QyGap1B);
        fmultGap1B_lq->Fill((Double_t)centV0, multGap1B);

    if (Qvec < fCutSmallQperc)
      fResSP_sq->Fill((Double_t)centV0, res);
        fQxGap1A_sq->Fill((Double_t)centV0, QxGap1A);
        fQyGap1A_sq->Fill((Double_t)centV0, QyGap1A);
        fmultGap1A_sq->Fill((Double_t)centV0, multGap1A);
        fQxGap1B_sq->Fill((Double_t)centV0, QxGap1B);
        fQyGap1B_sq->Fill((Double_t)centV0, QyGap1B);
        fmultGap1B_sq->Fill((Double_t)centV0, multGap1B);
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
