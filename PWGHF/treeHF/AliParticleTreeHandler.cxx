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
 * \class AliParticleTreeHandler
 * \brief Helper class to handle a tree for cut optimisation and MVA analyses, based heavily on AliHFTreeHandler and AliAnalysisTaskDmesonJets
 *
 * \author James Mulligan <james.mulligan@berkeley.edu>
 * \date June 9 2019
 */

#include <AliTLorentzVector.h>

#include "AliParticleTreeHandler.h"

//________________________________________________________________
/// \cond CLASSIMP
ClassImp(AliParticleTreeHandler);
/// \endcond

//________________________________________________________________
// Default constructor
AliParticleTreeHandler::AliParticleTreeHandler():
  TObject(),
  fTreeParticle(nullptr),
  fParticleContainer(nullptr),
  fParticlePt(-999.),
  fParticleEta(-999.),
  fParticlePhi(-999.),
  fRunNumber(0),
  fEventID(0), 
  fEventIDExt(0), 
  fEventIDLong(0) 
{
}

//________________________________________________________________
// Destructor
AliParticleTreeHandler::~AliParticleTreeHandler()
{
  if(fTreeParticle) delete fTreeParticle;
}

/**
 * Create particle TTree, with completely flat structure.
 */
//________________________________________________________________
TTree* AliParticleTreeHandler::BuildTree(TString name, TString title)
{
  if(fTreeParticle) {
    delete fTreeParticle;
    fTreeParticle=0x0;
  }
  fTreeParticle = new TTree(name.Data(),title.Data());
  
  // Create branches for each particle variable
  fTreeParticle->Branch("run_number", &fRunNumber);
  fTreeParticle->Branch("ev_id",&fEventID);
  fTreeParticle->Branch("ev_id_ext",&fEventIDExt);
  fTreeParticle->Branch("ev_id_long",&fEventIDLong);
  fTreeParticle->Branch("ParticlePt",&fParticlePt);
  fTreeParticle->Branch("ParticleEta",&fParticleEta);
  fTreeParticle->Branch("ParticlePhi",&fParticlePhi);
  
  return fTreeParticle;
}

/**
 * Set tree variables and fill them
 */
//________________________________________________________________
void AliParticleTreeHandler::FillTree(int runNumber, int eventID, int eventID_Ext, Long64_t eventID_Long)
{
  
  fRunNumber = runNumber;
  fEventID = eventID;
  fEventIDExt = eventID_Ext;
  fEventIDLong = eventID_Long;
  
  AliTLorentzVector partVec;
  for (const auto particleIterator : fParticleContainer->accepted_momentum()) {

    // Get particle four-vector
    partVec.Clear();
    partVec = particleIterator.first;
    
    // Set particle variables
    fParticleEta = partVec.Eta();
    fParticlePhi = partVec.Phi_0_2pi();
    fParticlePt = partVec.Pt();
    
    // Fill jet tree
    fTreeParticle->Fill();

  }
  
}
