/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//=========================================================================//

#include "TMath.h"
#include "TAxis.h"
#include "TSystem.h" 
#include "TProfile.h" 
#include "TH2F.h" 
#include "TH3F.h" 
#include "TFile.h" 
#include "TPRegexp.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliTracker.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliEbyEPidRatioQA.h"

using namespace std;

ClassImp(AliEbyEPidRatioQA)
//________________________________________________________________________
AliEbyEPidRatioQA::AliEbyEPidRatioQA() :
  AliEbyEPidRatioBase("QA", "QA"),

  fHnQA(NULL) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliEbyEPidRatioQA",10);
}

//________________________________________________________________________
AliEbyEPidRatioQA::~AliEbyEPidRatioQA() {
  // Destructor

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioQA::CreateHistograms() {
 
  Int_t    binHnQA[17] = {AliEbyEPidRatioHelper::fgkfHistNBinsCent, AliEbyEPidRatioHelper::fgkfHistNBinsEta,  //       cent |       eta
			  AliEbyEPidRatioHelper::fgkfHistNBinsRap,  AliEbyEPidRatioHelper::fgkfHistNBinsPhi,  //          y |       phi
			  AliEbyEPidRatioHelper::fgkfHistNBinsPt,   AliEbyEPidRatioHelper::fgkfHistNBinsPt,   //         pt |    pInner
			  AliEbyEPidRatioHelper::fgkfHistNBinsSign, 500,                                      //       sign | TPCsignal |  
			  50, 50, 50,                                                                         //  nSigmaITS | nSigmaTPC |  nSigmaTOF
			  50, 50, 50, 50,                                                                     //  DCAr      |      DCAz | nSigmaDCAr | nSigmaDCAz 
			  3,4};                                                                                 //  MCisProbe  
  
  Double_t minHnQA[17] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[0], AliEbyEPidRatioHelper::fgkfHistRangeEta[0], 
			  AliEbyEPidRatioHelper::fgkfHistRangeRap[0],  AliEbyEPidRatioHelper::fgkfHistRangePhi[0], 
			  AliEbyEPidRatioHelper::fgkfHistRangePt[0],   AliEbyEPidRatioHelper::fgkfHistRangePt[0],   
			  AliEbyEPidRatioHelper::fgkfHistRangeSign[0], 30, 
			  -10., -10., -10.,
			  -10., -10., -10., -10., 
			  -0.5,-0.5};
  
  Double_t maxHnQA[17] = {AliEbyEPidRatioHelper::fgkfHistRangeCent[1], AliEbyEPidRatioHelper::fgkfHistRangeEta[1], 
			  AliEbyEPidRatioHelper::fgkfHistRangeRap[1],  AliEbyEPidRatioHelper::fgkfHistRangePhi[1], 
			  AliEbyEPidRatioHelper::fgkfHistRangePt[1],   AliEbyEPidRatioHelper::fgkfHistRangePt[1],
			  AliEbyEPidRatioHelper::fgkfHistRangeSign[1], 500,
			  10., 10., 10.,
			  10.,  10., 10., 10., 
			  2.5,3.5};
  

  fHnQA = new THnSparseF("hnQA", "cent:eta:y:phi:pt:pInner:sign:TPCsignal:nSigmaITS:nSigmaTPC:nSigmaTOF:DCAr:DCAz:nSigmaDCAr:nSigmaDCAz:MCisProbe:PID",
				 17, binHnQA, minHnQA, maxHnQA);
  
  fHnQA->Sumw2();
  
  fHnQA->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnQA->GetAxis(1)->SetTitle("#eta");                         //  eta  [-0.9, 0.9]
  fHnQA->GetAxis(2)->SetTitle("#it{y}");                       //  rapidity [-0.5, 0.5]
  fHnQA->GetAxis(3)->SetTitle("#varphi");                      //  phi  [ 0. , 2Pi]
  fHnQA->GetAxis(4)->SetTitle("#it{p}_{T} (GeV/#it{c})");      //  pt   [ 0.46, 2.22]
  fHnQA->GetAxis(5)->SetTitle("#it{p}_{Inner} (GeV/#it{c})");  //  pInner [ 0.1, 3.0]
  fHnQA->GetAxis(6)->SetTitle("sign");                         //  -1 | 0 | +1 
  
  fHnQA->GetAxis(7)->SetTitle("TPC signal");                   //  TPCsignal [30, 500]
  fHnQA->GetAxis(8)->SetTitle("n #sigma ITS");                 //  nSigma ITS [-10, 10]
  fHnQA->GetAxis(9)->SetTitle("n #sigma TPC");                 //  nSigma TPC [-10, 10]
  fHnQA->GetAxis(10)->SetTitle("n #sigma TOF");                //  nSigma TOF [-10, 10]
  
  fHnQA->GetAxis(11)->SetTitle("DCAr");                        //  DCAr [-10, 10]
  fHnQA->GetAxis(12)->SetTitle("DCAz");                        //  DCAz [-10, 10]
  fHnQA->GetAxis(13)->SetTitle("n #sigma #sqrt(Cdd)/DCAr");    //  nSigma DCAr [-10, 10]
  fHnQA->GetAxis(14)->SetTitle("n #sigma #sqrt(Czz)/DCAz");    //  nSigma DCAz [-10, 10]
  fHnQA->GetAxis(15)->SetTitle("MCisProbe");                   //  0 | 1 (isProbeParticle) | 2 (isProbeParticle wrong sign) 
  fHnQA->GetAxis(16)->SetTitle("PID");                         //  0 | 1 | 2 | 3
  
  fHelper->BinLogAxis(fHnQA, 4);
  fHelper->BinLogAxis(fHnQA, 5);

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioQA::Process() {
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);
  
  // -- Track Loop
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted - ESD
    if (fESD && !fESDTrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
      continue;
    
    // -- Check if accepted - AOD
    if (fAOD){
      AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
      
      if (!trackAOD) {
	AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
	continue;
      }
      if (!trackAOD->TestFilterBit(fAODtrackCutBit))
	continue;
      
      // -- Check if in pT and eta range (is done in ESDTrackCuts for ESDs)
      if(!(track->Pt() > ptRange[0] && track->Pt() <= ptRange[1] && TMath::Abs(track->Eta()) <= etaRange[1]))
	continue;
    }

    Int_t gPdgCode = 0;

    Int_t iPid = 0;
    Double_t pid[3];
    if      (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kPion)))  {  iPid = 1; gPdgCode = 211;}
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kKaon)))  {  iPid = 2; gPdgCode = 321;}
    else if (fHelper->IsTrackAcceptedPID(track, pid, (AliPID::kProton))){  iPid = 3; gPdgCode = 2212;}
    else iPid = 0;

    Double_t yP;
    if ((iPid != 0) && !fHelper->IsTrackAcceptedRapidity(track, yP, iPid))
      continue;

    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;
    
    // -- Check if probe particle
    Int_t isProbeParticle = 0; 
    if (fIsMC) {
      Int_t label  = TMath::Abs(track->GetLabel());
      
      AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
      if (particle) {
	if (TMath::Abs(particle->PdgCode()) == gPdgCode) {
	  ++isProbeParticle;
	  if (particle->PdgCode() != (track->Charge()*gPdgCode))
	    ++isProbeParticle;
	}
      }
    }

    // -- Get dca r/z
    Float_t dca[] = {0.,0.};     // dca_xy, dca_z,
    Float_t cov[] = {0.,0.,0.}; // sigma_xy, sigma_xy_z, sigma_z
    if (fESD)
      (static_cast<AliESDtrack*>(track))->GetImpactParameters(dca,cov);

    Float_t dcaRoverCdd = ( TMath::Sqrt(cov[0]) != 0. )  ? dca[0]/TMath::Sqrt(cov[0]) : -9.99;
    Float_t dcaZoverCzz = ( TMath::Sqrt(cov[2]) != 0. )  ? dca[1]/TMath::Sqrt(cov[2]) : -9.99;

    if (iPid == 0)  
      yP = track->Eta();
   

    // cout << pid[0] << "  " << pid[1] << " " << pid[2] << yP << "  " << iPid << endl;

    if (iPid != 0) { 
      Double_t aTrack[17] = {
	Double_t(fCentralityBin),               //  0 centrality 
	track->Eta(),                           //  1 eta
	yP,                                     //  2 rapidity
	track->Phi(),                           //  3 phi
	track->Pt(),                            //  4 pt
	track->GetTPCmomentum(),                //  5 pInner
	track->Charge(),                        //  6 sign
	track->GetTPCsignal(),                  //  7 TPC dE/dx
	pid[0],                                 //  8 n Sigma ITS
	pid[1],                                 //  9 n Sigma TPC
	pid[2],                                 // 10 n Sigma TOF
	dca[0],                                 // 11 dca r
	dca[1],                                 // 12 dca z
	dcaRoverCdd,                            // 13 sqrt(cov[dd])
	dcaZoverCzz,                            // 14 sqrt(cov[zz])
	isProbeParticle,                        // 15 isProbe
	0                                    // 16 PID
      };
      
      fHnQA->Fill(aTrack);
    }
    
    Double_t aTrack[17] = {
      Double_t(fCentralityBin),               //  0 centrality 
      track->Eta(),                           //  1 eta
      yP,                                     //  2 rapidity
      track->Phi(),                           //  3 phi
      track->Pt(),                            //  4 pt
      track->GetTPCmomentum(),                //  5 pInner
      track->Charge(),                        //  6 sign
      track->GetTPCsignal(),                  //  7 TPC dE/dx
      pid[0],                                 //  8 n Sigma ITS
      pid[1],                                 //  9 n Sigma TPC
      pid[2],                                 // 10 n Sigma TOF
      dca[0],                                 // 11 dca r
      dca[1],                                 // 12 dca z
      dcaRoverCdd,                            // 13 sqrt(cov[dd])
      dcaZoverCzz,                            // 14 sqrt(cov[zz])
      isProbeParticle,                        // 15 isProbe
      iPid                                    // 16 PID
    };

   fHnQA->Fill(aTrack);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}
