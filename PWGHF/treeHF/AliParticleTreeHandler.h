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

#ifndef ALIPARTICLETREEHANDLER_H
#define ALIPARTICLETREEHANDLER_H

/**
 * \class AliParticleTreeHandler
 * \brief Helper class to handle a tree for cut optimisation and MVA analyses, based heavily on AliHFTreeHandler and AliAnalysisTaskDmesonJets
 *
 * \author James Mulligan <james.mulligan@berkeley.edu>
 * \date June 9 2019
 */

#include "TTree.h"

#include "AliParticleContainer.h"

//________________________________________________________________
//****************************************************************
class AliParticleTreeHandler : public TObject
{
  public:

    AliParticleTreeHandler();
    virtual ~AliParticleTreeHandler();

    // Core methods
    TTree* BuildTree(TString name, TString title);
    void FillTree(int runNumber, int eventID, int eventID_Ext, Long64_t eventID_Long);
  
    // Setters
    void SetParticleContainer(AliParticleContainer* partCont) { fParticleContainer = partCont; }

  protected:
  
    // The TreeHandler stores:
    //     - A particle tree, filled once per particle (i.e. track or MCParticle).
    //
    // The trees contain: event id, run number, particle pT, particle eta, particle phi.
    //
    // Each tree structure is completely flat -- the branches are all primitive types.
  
    TTree*                       fTreeParticle;            ///< Tree with particle variables
  
    AliParticleContainer*        fParticleContainer;       //!<! Particle container for this tree
  
    // Track quantities.
    float                        fParticlePt;              //!<! Pt of particle
    float                        fParticleEta;             //!<! Eta of particle
    float                        fParticlePhi;             //!<! Phi of particle (0 < phi < 2pi)
  
    // Event quantities
    int                          fRunNumber;               //!<! run number
    int                          fEventID;                 //!<! event ID (unique identifier when run number is fixed), first 32 bits of fEventIDLong
    int                          fEventIDExt;                 //!<! event ID (unique identifier when run number is fixed), second 32 bits of fEventIDLong
    Long64_t                     fEventIDLong;                 //!<! event ID (unique identifier when run number is fixed), full 64 bits of fEventIDLong

  /// \cond CLASSIMP
  ClassDef(AliParticleTreeHandler,1); ///
  /// \endcond
};

#endif
