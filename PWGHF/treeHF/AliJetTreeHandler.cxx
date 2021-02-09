/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/

/**
 * \class AliJetTreeHandler
 * \brief Helper class to handle a tree for cut optimisation and MVA analyses, based heavily on AliHFTreeHandler and AliAnalysisTaskDmesonJets
 *
 * \author James Mulligan <james.mulligan@berkeley.edu>
 * \date Feb 15 2019
 */

#include <TVector2.h>

#include "AliJetTreeHandler.h"

//________________________________________________________________
/// \cond CLASSIMP
ClassImp(AliJetTreeHandler);
/// \endcond

//________________________________________________________________
// Default constructor
AliJetTreeHandler::AliJetTreeHandler():
  TObject(),
  fTreeJet(nullptr),
  fTreeJetConstituent(nullptr),
  fFillJetConstituentTree(false),
  fJetContainer(nullptr),
  fMinJetPtCorr(0.),
  fFillJetEtaPhi(false),
  fFillPtCorr(false),
  fFillPtUncorr(false),
  fFillArea(false),
  fFillNConstituents(false),
  fFillZLeading(false),
  fFillRadialMoment(false),
  fFillpTD(false),
  fFillMass(false),
  fFillMatchingJetID(false),
  fTrackPt(-999.),
  fTrackEta(-999.),
  fTrackPhi(-999.),
  fRunNumber(0),
  fEventID(0),
  fJetID(999),
  fPtCorr(-999.),
  fEta(-999.),
  fPhi(-999.),
  fPtUncorr(-999.),
  fArea(-999.),
  fN(0),
  fZLeading(-999.),
  fRadialMoment(-999.),
  fpTD(-999.),
  fMass(-999.),
  fMatchedJetID(999)
{
}

//________________________________________________________________
// Destructor
AliJetTreeHandler::~AliJetTreeHandler()
{
  if(fTreeJet) delete fTreeJet;
  if(fTreeJetConstituent) delete fTreeJetConstituent;
}

/**
 * Create jet TTree, with completely flat structure.
 */
//________________________________________________________________
TTree* AliJetTreeHandler::BuildJetTree(TString name, TString title)
{
  if(fTreeJet) {
    delete fTreeJet;
    fTreeJet=0x0;
  }
  fTreeJet = new TTree(name.Data(),title.Data());
  
  // Create branches for each jet variable
  
  fTreeJet->Branch("run_number", &fRunNumber);
  fTreeJet->Branch("ev_id",&fEventID);
  fTreeJet->Branch("jet_id",&fJetID);

  if (fFillJetEtaPhi) {
    fTreeJet->Branch("Eta",&fEta);
    fTreeJet->Branch("Phi",&fPhi);
  }
  
  if (fFillPtCorr) {
    fTreeJet->Branch("PtCorr",&fPtCorr);
  }

  if (fFillPtUncorr) {
    fTreeJet->Branch("PtUncorr",&fPtUncorr);
  }
  
  if (fFillArea) {
    fTreeJet->Branch("Area",&fArea);
  }
  
  if (fFillNConstituents) {
    fTreeJet->Branch("N",&fN);
  }
  
  if (fFillZLeading) {
    fTreeJet->Branch("ZLeading", &fZLeading);
  }
  
  if (fFillRadialMoment) {
    fTreeJet->Branch("RadialMoment", &fRadialMoment);
  }
  
  if (fFillpTD) {
    fTreeJet->Branch("pTD", &fpTD);
  }
  
  if (fFillMass) {
    fTreeJet->Branch("Mass", &fMass);
  }
  
  if (fFillMatchingJetID) {
    fTreeJet->Branch("MatchedJetID", &fMatchedJetID);
  }
  
  return fTreeJet;
}

/**
 * Create jet TTree, with completely flat structure.
 */
//________________________________________________________________
TTree* AliJetTreeHandler::BuildJetConstituentTree(TString name, TString title)
{
  if(fTreeJetConstituent) {
    delete fTreeJetConstituent;
    fTreeJetConstituent=0x0;
  }
  fTreeJetConstituent = new TTree(name.Data(),title.Data());
  
  // Create branches for each jet variable
  fTreeJetConstituent->Branch("run_number", &fRunNumber);
  fTreeJetConstituent->Branch("ev_id",&fEventID);
  fTreeJetConstituent->Branch("jet_id",&fJetID);
  fTreeJetConstituent->Branch("TrackPt",&fTrackPt);
  fTreeJetConstituent->Branch("TrackEta",&fTrackEta);
  fTreeJetConstituent->Branch("TrackPhi",&fTrackPhi);
 
  return fTreeJetConstituent;
}

/**
 * Set tree variables and fill them
 */
//________________________________________________________________
void AliJetTreeHandler::FillTree(int runNumber, int eventID, int eventID_Ext, Long64_t eventID_Long)
{
  
  fRunNumber = runNumber;
  fEventID = eventID;
  fEventIDExt = eventID_Ext;
  fEventIDLong = eventID_Long;
  
  for (const auto jet : fJetContainer->accepted()) {
    
    // Check if jet passes min pt threshold
    if (GetJetPt(jet) < fMinJetPtCorr) {
      continue;
    }
    
    /////////////////////////////////////////////
    // Fill jet tree
    
    // Set jet variables
    fJetID = jet->GetLabel();
    SetJetVariables(jet);
    
    // Fill jet tree
    fTreeJet->Fill();
    
    /////////////////////////////////////////////
    // Fill jet constituent tree (if enabled)
    
    if (fFillJetConstituentTree) {
      
      // Loop through tracks
      const int nTracks = jet->GetNumberOfTracks();
      for (int i = 0; i < nTracks; ++i) {
        
        // Set jet constituent variables
        const AliVParticle* track = jet->Track(i);
        SetJetConstituentVariables(track);
        
        // Fill jet constituent tree
        fTreeJetConstituent->Fill();
        
      }
    }
    
  }
  
}

/**
 * Set jet tree variables
 */
//________________________________________________________________
void AliJetTreeHandler::SetJetVariables(const AliEmcalJet* jet)
{
  
  // Set jet variables
  if (fFillJetEtaPhi) {
    fEta = jet->Eta();
    fPhi = jet->Phi_0_2pi();
  }
  if (fFillPtCorr) {
    fPtCorr = GetJetPt(jet);
  }
  if (fFillPtUncorr) {
    fPtUncorr = jet->Pt();
  }
  if (fFillArea) {
    fArea = jet->Area();
  }
  if (fFillNConstituents) {
    fN = static_cast<unsigned short int>(jet->GetNumberOfConstituents());
  }
  if (fFillZLeading) {
    fZLeading = fJetContainer->GetZLeadingCharged(jet);
  }
  if (fFillRadialMoment) {
    fRadialMoment = RadialMoment(jet);
  }
  if (fFillpTD) {
    fpTD = PTD(jet);
  }
  if (fFillMass) {
    fMass = jet->M();
  }
  
  // Get matched jet (assumes the matches have been filled by a previous task)
  if (fFillMatchingJetID) {
    
    unsigned short int matchedJetLabel = 999;
    const AliEmcalJet* matchedJet = jet->ClosestJet();
    if (matchedJet) {
      matchedJetLabel = static_cast<unsigned short int>(matchedJet->GetLabel());
    }
    fMatchedJetID = matchedJetLabel;
    
  }

}

/**
* Set jet constituent tree variables
*/
//________________________________________________________________
void AliJetTreeHandler::SetJetConstituentVariables(const AliVParticle* track)
{
  
  fTrackPt = track->Pt();
  fTrackEta = track->Eta();
  fTrackPhi = TVector2::Phi_0_2pi(track->Phi());
  
}

/**
 * Set Jet ID for accepted jets in each event: 0, 1, 2, ..., N
 * If jet is not accepted, it will have ID 999.
 */
//________________________________________________________________
void AliJetTreeHandler::SetJetLabels()
{
  // Reset all labels
  for (auto jet : fJetContainer->all()) {
    jet->SetLabel(999);
  }

  // Assign jet labels
  unsigned short int jetID = 0;
  for (auto jet : fJetContainer->accepted()) {
    
    // Check if jet passes min pt threshold
    if (GetJetPt(jet) < fMinJetPtCorr) {
      continue;
    }
    
    jet->SetLabel(jetID);
    jetID++;
  }
}

/**
 * Get pT of jet -- background subtracted
 */
//________________________________________________________________
Double_t AliJetTreeHandler::GetJetPt(const AliEmcalJet* jet)
{
  Float_t rhoVal = 0;
  if (fJetContainer->GetRhoParameter()) {
    rhoVal = fJetContainer->GetRhoVal();
  }
  
  Float_t pTCorr = jet->Pt() - rhoVal * jet->Area();
  return pTCorr;
}

//________________________________________________________________
Double_t AliJetTreeHandler::PTD(const AliEmcalJet *jet){
  
  Double_t numeratorSquared=0;
  Double_t denominator=0;
  
  for (Int_t i=0; i< jet->GetNumberOfTracks(); i++){
    
    AliVParticle *particle = static_cast<AliVParticle*>(jet->Track(i));
    if (!particle) continue;
    
    numeratorSquared += particle->Pt() * particle->Pt();
    denominator += particle->Pt();
    
  }
  Double_t numerator = TMath::Sqrt(numeratorSquared);
  
  if(TMath::Abs(denominator) > 1e-3) {
    return numerator/denominator;
  }
  else {
    return -1;
  }
}

//________________________________________________________________
Double_t AliJetTreeHandler::RadialMoment(const AliEmcalJet* jet){
  
  Double_t numerator = 0;
  Double_t denominator=0;
  
  for (Int_t i=0; i< jet->GetNumberOfTracks(); i++){
    
    const AliVParticle* particle = static_cast<AliVParticle*>(jet->Track(i));
    if(!particle) continue;
    
    numerator += particle->Pt() * DeltaR(jet, particle);
    denominator += particle->Pt();
    
  }
  
  if(TMath::Abs(denominator) > 1e-3) {
    return numerator/denominator;
  }
  else {
    return -1;
  }
  
}

//________________________________________________________________
Double_t AliJetTreeHandler::DeltaR(const AliEmcalJet* jet, const AliVParticle* part) {
  
  Double_t jetEta = jet->Eta();
  Double_t jetPhi = jet->Phi_0_2pi();

  Double_t partEta = part->Eta();
  Double_t partPhi = TVector2::Phi_0_2pi(part->Phi());
  
  return TMath::Sqrt( (jetEta-partEta)*(jetEta-partEta) + (jetPhi-partPhi)*(jetPhi-partPhi) );
  
}
