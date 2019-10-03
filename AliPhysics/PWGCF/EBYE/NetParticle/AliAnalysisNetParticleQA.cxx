//-*- Mode: C++ -*-

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

#include "AliAnalysisNetParticleQA.h"

using namespace std;

/**
 * Class for for NetParticle QA
 * -- Create input for QA
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

ClassImp(AliAnalysisNetParticleQA)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleQA::AliAnalysisNetParticleQA() :
  AliAnalysisNetParticleBase("QA", "QA"),

  fHnQA(NULL) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliAnalysisNetParticleQA",10);
}

//________________________________________________________________________
AliAnalysisNetParticleQA::~AliAnalysisNetParticleQA() {
  // Destructor

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */


//________________________________________________________________________
void AliAnalysisNetParticleQA::CreateHistograms() {
  // -- Add histograms to outlist

  // ------------------------------------------------------------------
  // -- Create THnSparseF - QA
  // ------------------------------------------------------------------

  Int_t    binHnQA[16] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent, AliAnalysisNetParticleHelper::fgkfHistNBinsEta,  //       cent |       eta
			  AliAnalysisNetParticleHelper::fgkfHistNBinsRap,  AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,  //          y |       phi
			  AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   //         pt |    pInner
			  AliAnalysisNetParticleHelper::fgkfHistNBinsSign, 500,                                             //       sign | TPCsignal |  
			  50, 50, 50,                                                                                       //  nSigmaITS | nSigmaTPC |  nSigmaTOF
			  50, 50, 50, 50,                                                                                   //  DCAr      |      DCAz | nSigmaDCAr | nSigmaDCAz 
			  3};                                                                                               //  MCisProbe  
  
  Double_t minHnQA[16] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			  AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   
			  AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 30, 
			  -10., -10., -10.,
			  -10., -10., -10., -10., 
			  -0.5};
  
  Double_t maxHnQA[16] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			  AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			  AliAnalysisNetParticleHelper::fgkfHistRangeSign[1], 500,
			  10., 10., 10.,
			  10.,  10., 10., 10., 
			  2.5};
  

  fHnQA = new THnSparseF("hnQA", "cent:eta:y:phi:pt:pInner:sign:TPCsignal:nSigmaITS:nSigmaTPC:nSigmaTOF:DCAr:DCAz:nSigmaDCAr:nSigmaDCAz:MCisProbe",
				 16, binHnQA, minHnQA, maxHnQA);
  
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
  
  fHelper->BinLogAxis(fHnQA, 4);
  fHelper->BinLogAxis(fHnQA, 5);

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleQA::Process() {
  // -- Process ESD/AOD tracks and fill QA histogram

  // -- Get ranges for AOD particles
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);
  
  // -- Track Loop
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Check track cuts
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

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

    // -- Check if accepted in rapidity window -- for identified particles
    //    for charged particles -- eta check is done in the trackcuts
    Double_t yP;
    if (fHelper->GetUsePID() && !fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted with thighter DCA cuts
    // -- returns kTRUE for AODs for now : MW
    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;

    // -- Check if accepted by PID from TPC or TPC+TOF
    Double_t pid[3];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Fill Track
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    // -- Check if probe particle
    Int_t isProbeParticle = 0; 
    if (fIsMC) {
      Int_t label  = TMath::Abs(track->GetLabel());
      
      AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
      if (particle) {
	if (TMath::Abs(particle->PdgCode()) == fPdgCode) {
	  ++isProbeParticle;
	  if (particle->PdgCode() != (track->Charge()*fPdgCode))
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

    if (!fHelper->GetUsePID())  
      yP = track->Eta();
    
    Double_t aTrack[16] = {
      Double_t(fCentralityBin),               //  0 centrality 
      track->Eta(),                           //  1 eta
      yP,                                     //  2 rapidity
      track->Phi(),                           //  3 phi
      track->Pt(),                            //  4 pt
      track->GetTPCmomentum(),                //  5 pInner
      static_cast<Double_t>(track->Charge()),                        //  6 sign
      track->GetTPCsignal(),                  //  7 TPC dE/dx
      pid[0],                                 //  8 n Sigma ITS
      pid[1],                                 //  9 n Sigma TPC
      pid[2],                                 // 10 n Sigma TOF
      dca[0],                                 // 11 dca r
      dca[1],                                 // 12 dca z
      dcaRoverCdd,                            // 13 sqrt(cov[dd])
      dcaZoverCzz,                            // 14 sqrt(cov[zz])
      static_cast<Double_t>(isProbeParticle)                         // 15 isProbe
    };

   fHnQA->Fill(aTrack);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}
