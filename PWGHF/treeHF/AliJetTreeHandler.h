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

#ifndef ALIJETTREEHANDLER_H
#define ALIJETTREEHANDLER_H

/**
 * \class AliJetTreeHandler
 * \brief Helper class to handle a tree for cut optimisation and MVA analyses, based heavily on AliHFTreeHandler and AliAnalysisTaskDmesonJets
 *
 * \author James Mulligan <james.mulligan@berkeley.edu>
 * \date Feb 15 2019
 */

#include "TTree.h"

#include "AliJetContainer.h"

//________________________________________________________________
//****************************************************************
class AliJetTreeHandler : public TObject
{
  public:

    AliJetTreeHandler();
    virtual ~AliJetTreeHandler();

    // Core methods
    TTree* BuildTree(TString name, TString title);
    bool SetJetVariables();
    void FillTree(); //to be called for each event
  
    // Setters
    void SetJetContainer(AliJetContainer* jetCont) { fJetContainer = jetCont; }
    void SetFillPtUncorr(bool b) { fFillPtUncorr = b; }
    void SetFillArea(bool b) { fFillArea = b; }
    void SetFillNConstituents(bool b) { fFillNConstituents = b; }
    void SetFillZLeading(bool b) { fFillZLeading = b; }
    void SetFillRadialMoment(bool b) { fFillRadialMoment = b; }
    void SetFillpTD(bool b) { fFillpTD = b; }
    void SetFillMass(bool b) { fFillMass = b; }
    void SetFillMatchingJetID(bool b) { fFillMatchingJetID = b; }
  
    // Utility functions
    void SetJetLabels();
    Double_t GetJetPt(const AliEmcalJet* jet);
    Double_t PTD(const AliEmcalJet *jet);
    Double_t RadialMoment(const AliEmcalJet* jet);
    Double_t DeltaR(const AliEmcalJet* jet, const AliVParticle* part);

  protected:
  
    TTree*                       fTreeVar;                 ///< Tree with compressed jet objects
    AliJetContainer*             fJetContainer;            //!<! Jet container for this tree
  
    // Flags specifying what info to fill to the tree
    bool                         fFillPtUncorr;            ///< Pt of the jet (GeV/c) (not background subtracted)
    bool                         fFillArea;                ///< Area
    bool                         fFillNConstituents;       ///< N constituents
    bool                         fFillZLeading;            ///< ZLeading
    bool                         fFillRadialMoment;        ///< Radial moment
    bool                         fFillpTD;                 ///< pT,D
    bool                         fFillMass;                ///< Mass
    bool                         fFillMatchingJetID;       ///< jet matching

    // Jet parameters to be stored in the tree.
    // Each branch in the tree consists of a vector of a given variable,
    // where the i^th entry in the vector corresponds to the i^th jet.
  
    // Basic jet quantities (always filled)
    std::vector<float>           fPtCorr;                  //!<! Pt of the jet after subtracting average background
    std::vector<float>           fEta;                     //!<! Eta of the jet
    std::vector<float>           fPhi;                     //!<! Phi of the jet (0 < phi < 2pi)

    // Other jet quantities
    std::vector<float>           fPtUncorr;                //!<! Pt of the jet (GeV/c) (not background subtracted)
    std::vector<float>           fArea;                    //!<! Area of the jet
  
    // Jet substructure observables
    std::vector<float>           fN;                       //!<! Number of jet constituents
    std::vector<float>           fZLeading;                //!<! z of leading track
    std::vector<float>           fRadialMoment;            //!<! Radial moment (not background subtracted)
    std::vector<float>           fpTD;                     //!<! Momentum dispersion (not background subtracted)
    std::vector<float>           fMass;                    //!<! Jet mass (not background subtracted)

    // Jet matching
    std::vector<float>           fMatchedJetID;            //!<! Index of matched jet in the matching container's std::vectors
  
  /// \cond CLASSIMP
  ClassDef(AliJetTreeHandler,1); ///
  /// \endcond
};

#endif
