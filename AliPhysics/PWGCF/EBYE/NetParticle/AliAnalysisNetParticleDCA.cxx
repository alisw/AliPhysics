//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisNetParticleDCA.h"

using namespace std;

/**
 * Class for NetParticle Distributions
 * -- DCA distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

ClassImp(AliAnalysisNetParticleDCA)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleDCA::AliAnalysisNetParticleDCA() :
  AliAnalysisNetParticleBase("DCA", "DCA"),
  fESDTrackCutsBkg(NULL),
  fHnDCA(NULL) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliAnalysisNetParticleDCA",10);
}

//________________________________________________________________________
AliAnalysisNetParticleDCA::~AliAnalysisNetParticleDCA() {
  // Destructor
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleDCA::Process() {
  // -- Process event
  //    Process ESD/AOD tracks and fill DCA THnSparse

  // -- Get ranges for AOD particles
  Float_t etaRange[2];
  fESDTrackCutsBkg->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCutsBkg->GetPtRange(ptRange[0],ptRange[1]); 

  // -- Track Loop
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted - ESD
    if (fESD && !fESDTrackCutsBkg->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
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

    // -- Check if accepted by PID from TPC or TPC+TOF
    Double_t pid[3];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- Check if accepted with thighter DCA cuts
    // -- returns kTRUE for AODs for now : MW
    Bool_t isDCArAccepted = fHelper->IsTrackAcceptedDCA(track);

    // -- Check if accepted with thighter DCA cuts
    // ?!?!? How to mimic this in AODs?
    if (fESD && !fESDTrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
      isDCArAccepted = kFALSE;
    
    // -- Check for contamination 
    Int_t contIdx = (fIsMC) ? GetContIdxTrack(TMath::Abs(track->GetLabel()), track->Charge()) : 1;
    
    // -- Get DCAs (dca_r, dca_z, sigma_xy, sigma_xy_z, sigma_z)
    Float_t dca[2], cov[3]; // 
    if (fESD)
      (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
    else 
      dca[0] = 1.; 

    // -- Fill THnSparse
    Double_t hnDCA[9] = {fCentralityBin, track->Eta(), yP, track->Phi(), track->Pt(), static_cast<Double_t>(track->Charge()), static_cast<Double_t>(contIdx), static_cast<Double_t>(isDCArAccepted), dca[0]};
    fHnDCA->Fill(hnDCA);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}      

/*
 * ---------------------------------------------------------------------------------
 *                                 Private Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleDCA::CreateHistograms() {
  // -- Create histograms

  // ------------------------------------------------------------------
  // -- Create THnSparseD - DCA
  // ------------------------------------------------------------------

  Int_t    binHnDCA[9] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent, AliAnalysisNetParticleHelper::fgkfHistNBinsEta,        //     cent |       etaRec
			  AliAnalysisNetParticleHelper::fgkfHistNBinsRap,  AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,        //     yRec |       phiRec
			  AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   AliAnalysisNetParticleHelper::fgkfHistNBinsSign,       //    ptRec |         sign
			  4,  2,  77};                                                                                            //  contPart| DCArAccepted | DCAr

  Double_t minHnDCA[9] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			  AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 
			  0.5, -0.5, -3.};
  
  Double_t maxHnDCA[9] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			  AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1], 
			  4.5, 1.5, 3.};


  fHnDCA = new THnSparseD("hnDCA", "cent:etaRec:yRec:phiRec:ptRec:sign:contPart:contSign:DCArAccepted:DCAr", 9, binHnDCA, minHnDCA, maxHnDCA);
  
  fHnDCA->Sumw2();
  fHnDCA->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnDCA->GetAxis(1)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnDCA->GetAxis(2)->SetTitle("#it{y}_{Rec}");                 //  rapidity  [-0.5, 0.5]
  fHnDCA->GetAxis(3)->SetTitle("#varphi_{Rec} (rad)");          //  phi  [ 0. , 2Pi]
  fHnDCA->GetAxis(4)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})");  //  pT   [ 0.2, 2.6]
  fHnDCA->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 

  fHnDCA->GetAxis(6)->SetTitle("contPart");                     //  1  primary | 2 missId | 3 from WeakDecay | 4 p from Material
  fHnDCA->GetAxis(7)->SetTitle("DCArAccepted");                 //  0 not accepted | 1 accepted 
  fHnDCA->GetAxis(8)->SetTitle("DCAr");                         //  DCAr [-3, 3]

  fHelper->BinLogAxis(fHnDCA,  4, fESDTrackCuts);
  fHelper->BinLogAxis(fHnDCA,  4, fESDTrackCutsBkg);

  // -- Set binning for DCAr
  Double_t binsDCAr[77] = {-3.,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.285,-0.27,-0.255,-0.24,-0.225,-0.21,-0.195,-0.18,-0.165,-0.15,-0.135,-0.12,-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0.,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3.};
  fHnDCA->GetAxis(8)->Set(76, binsDCAr);

  // ------------------------------------------------------------------
  
  return;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDCA::GetContIdxTrack(Int_t label, Int_t sign) {
  // Check if particle is contamination or correctly identified for ESDs and AODs
  // Check for missidentified primaries and secondaries
  // return  1  primary | 2 missId | 3 from WeakDecay | 4 from Material | -1 unknown

  Int_t contIdx = -1;

  AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
  if (!particle)
    return contIdx;

  Bool_t isPhysicalPrimary = (fESD) ? fStack->IsPhysicalPrimary(label): (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
  Bool_t isSecondaryFromWeakDecay = (fESD) ? fStack->IsSecondaryFromWeakDecay(label) : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
  Bool_t isSecondaryFromMaterial  = (fESD) ? fStack->IsSecondaryFromMaterial(label)  : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
  
  // -- Check if primary 
  //    > if PID required check -> for the correct (signed pdgcode) particle
  //    > no PID just check for primary 
  if (isPhysicalPrimary) {
    if (fHelper->GetUsePID()) {
      // -- Check if correctly identified 
      if (particle->PdgCode() == (sign*fPdgCode))
	contIdx = 1;   
      // -- MissIdentification
      else 
	contIdx = 2;
    }
    else
      contIdx = 1;   
  }

  // -- Check if secondaries from material or weak decay
  else if(isSecondaryFromWeakDecay)
    contIdx = 3;
  else if (isSecondaryFromMaterial)
    contIdx = 4;
  else
    contIdx = -1;

  return contIdx;
}



