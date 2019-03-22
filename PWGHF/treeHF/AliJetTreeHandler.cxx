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
  fTreeVar(nullptr),
  fJetContainer(nullptr),
  fFillPtUncorr(false),
  fFillArea(true),
  fFillNConstituents(true),
  fFillZLeading(true),
  fFillRadialMoment(true),
  fFillpTD(true),
  fFillMass(true),
  fFillMatchingJetID(false),
  fPtCorr(),
  fEta(),
  fPhi(),
  fPtUncorr(),
  fArea(),
  fN(),
  fZLeading(),
  fRadialMoment(),
  fpTD(),
  fMass(),
  fMatchedJetID()
{
}

//________________________________________________________________
// Destructor
AliJetTreeHandler::~AliJetTreeHandler()
{
  if(fTreeVar) delete fTreeVar;
}

/**
 * Create jet TTree, with a branch of vectors for each jet variable.
 * There will be one entry in each vector for each jet that is found.
 */
//________________________________________________________________
TTree* AliJetTreeHandler::BuildTree(TString name, TString title)
{
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=0x0;
  }
  fTreeVar = new TTree(name.Data(),title.Data());
  
  // Create branches for each jet variable
  fTreeVar->Branch("PtCorr",&fPtCorr);
  fTreeVar->Branch("Eta",&fEta);
  fTreeVar->Branch("Phi",&fPhi);
  
  if (fFillPtUncorr) {
    fTreeVar->Branch("PtUncorr",&fPtUncorr);
  }
  
  if (fFillArea) {
    fTreeVar->Branch("Area",&fArea);
  }
  
  if (fFillNConstituents) {
    fTreeVar->Branch("N",&fN);
  }
  
  if (fFillZLeading) {
    fTreeVar->Branch("ZLeading", &fZLeading);
  }
  
  if (fFillRadialMoment) {
    fTreeVar->Branch("RadialMoment", &fRadialMoment);
  }
  
  if (fFillpTD) {
    fTreeVar->Branch("pTD", &fpTD);
  }
  
  if (fFillMass) {
    fTreeVar->Branch("Mass", &fMass);
  }
  
  if (fFillMatchingJetID) {
    fTreeVar->Branch("MatchedJetID", &fMatchedJetID);
  }
  
  return fTreeVar;
}

/**
 * Set jet tree variables
 */
//________________________________________________________________
bool AliJetTreeHandler::SetJetVariables()
{

  for (const auto jet : fJetContainer->accepted()) {
    
    fPtCorr.push_back(GetJetPt(jet));
    fEta.push_back(jet->Eta());
    fPhi.push_back(jet->Phi_0_2pi());
    
    if (fFillPtUncorr) {
      fPtUncorr.push_back(jet->Pt());
    }
    
    if (fFillArea) {
      fArea.push_back(jet->Area());
    }
    
    if (fFillNConstituents) {
      fN.push_back(jet->GetNumberOfConstituents());
    }
    
    if (fFillZLeading) {
      fZLeading.push_back(fJetContainer->GetZLeadingCharged(jet));
    }
    
    if (fFillRadialMoment) {
      fRadialMoment.push_back(RadialMoment(jet));
    }
    
    if (fFillpTD) {
      fpTD.push_back(PTD(jet));
    }
    
    if (fFillMass) {
      fMass.push_back(jet->M());
    }
    
    // Get matched jet (assumes the matches have been filled by a previous task)
    if (fFillMatchingJetID) {
      
      int matchedJetLabel = -1;
      const AliEmcalJet* matchedJet = jet->ClosestJet();
      if (matchedJet) {
        matchedJetLabel = matchedJet->GetLabel();
      }
      fMatchedJetID.push_back(matchedJetLabel);
      
    }

  }
  
  return true;
}

/**
 * Fill jet tree, and reset all vectors
 */
//________________________________________________________________
void AliJetTreeHandler::FillTree() {

  fTreeVar->Fill();
  
  // Reset all vectors
  fPtCorr.clear();
  fEta.clear();
  fPhi.clear();
  
  if (fFillPtUncorr) {
    fPtUncorr.clear();
  }
  
  if (fFillArea) {
    fArea.clear();
  }
  
  if (fFillNConstituents) {
    fN.clear();
  }
  
  if (fFillZLeading) {
    fZLeading.clear();
  }
  
  if (fFillRadialMoment) {
    fRadialMoment.clear();
  }
  
  if (fFillpTD) {
    fpTD.clear();
  }
  
  if (fFillMass) {
    fMass.clear();
  }
  
  if (fFillMatchingJetID) {
    fMatchedJetID.clear();
  }
  
}

/**
 *  If filling jet matching info (for MC), loop through jets and set fLabel,
 *  which specifies the index of the jet in the tree variable std::vectors.

 */
//________________________________________________________________
void AliJetTreeHandler::SetJetLabels()
{
  int i = 0;
  for (auto jet : fJetContainer->accepted()) {
    jet->SetLabel(i);
    i++;
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
