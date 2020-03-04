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
    TTree* BuildJetTree(TString name, TString title);
    TTree* BuildJetConstituentTree(TString name, TString title);
    void FillTree(int runNumber, int eventID, int eventID_Ext, Long64_t eventID_Long);
    void SetJetVariables(const AliEmcalJet* jet);
    void SetJetConstituentVariables(const AliVParticle* track);
    void SetJetLabels();
  
    // Setters
    void SetJetContainer(AliJetContainer* jetCont) { fJetContainer = jetCont; }
    void SetMinJetPtCorr(double pt) { fMinJetPtCorr = pt; }
    void SetFillJetConstituentTree(bool b) { fFillJetConstituentTree = b; }

    void SetFillJetEtaPhi(bool b) { fFillJetEtaPhi = b; }
    void SetFillPtCorr(bool b) { fFillPtCorr = b; }
    void SetFillPtUncorr(bool b) { fFillPtUncorr = b; }
    void SetFillArea(bool b) { fFillArea = b; }
    void SetFillNConstituents(bool b) { fFillNConstituents = b; }
    void SetFillZLeading(bool b) { fFillZLeading = b; }
    void SetFillRadialMoment(bool b) { fFillRadialMoment = b; }
    void SetFillpTD(bool b) { fFillpTD = b; }
    void SetFillMass(bool b) { fFillMass = b; }
    void SetFillMatchingJetID(bool b) { fFillMatchingJetID = b; }
  
    // Utility functions
    Double_t GetJetPt(const AliEmcalJet* jet);
    Double_t PTD(const AliEmcalJet *jet);
    Double_t RadialMoment(const AliEmcalJet* jet);
    Double_t DeltaR(const AliEmcalJet* jet, const AliVParticle* part);

  protected:
  
    // The TreeHandler stores:
    //     - A jet tree, filled once per jet. By default, it stores only: event id, jet id.
    //       Additional variables can be activated by the setters in the TreeCreator
    //     - (optional) A jet constituent tree, filled once per jet constituent.
    //
    // Each tree structure is completely flat -- the branches are all primitive types
  
    TTree*                       fTreeJet;                 ///< Tree with jet variables
    TTree*                       fTreeJetConstituent;      ///< Tree with jet constituent variables (optional)
  
    bool                         fFillJetConstituentTree;  ///< Store pT,eta,phi of all tracks inside the jet
  
    AliJetContainer*             fJetContainer;            //!<! Jet container for this tree
    double                       fMinJetPtCorr;            ///< Min jet Pt (background subtracted) to fill jet into tree
  
    // Fill jet tree according to the below flags. By default, it only contains: event id, jet id
    bool                         fFillJetEtaPhi;           ///< Jet eta, phi
    bool                         fFillPtCorr;              ///< Pt of the jet (GeV/c) (background subtracted)
    bool                         fFillPtUncorr;            ///< Pt of the jet (GeV/c) (not background subtracted)
    bool                         fFillArea;                ///< Area
    bool                         fFillNConstituents;       ///< N constituents
    bool                         fFillZLeading;            ///< ZLeading
    bool                         fFillRadialMoment;        ///< Radial moment
    bool                         fFillpTD;                 ///< pT,D
    bool                         fFillMass;                ///< Mass
    bool                         fFillMatchingJetID;       ///< jet matching
  
    // Jet constituent quantities.
    float                        fTrackPt;                 //!<! Pt of track
    float                        fTrackEta;                //!<! Eta of track
    float                        fTrackPhi;                //!<! Phi of track (0 < phi < 2pi)
  
    // Jet quantities
    int                          fRunNumber;               //!<! run number
    int                          fEventID;                 //!<! event ID (unique identifier when run number is fixed), first 32 bits of fEventIDLong
    int                          fEventIDExt;              //!<! event ID (unique identifier when run number is fixed), second 32 bits of fEventIDLong
    Long64_t                     fEventIDLong;             //!<! event ID (unique identifier when run number is fixed), full 4 bits of fEventIDLong
    unsigned short int           fJetID;                   //!<! jet ID (i^th jet in each event)

    float                        fPtCorr;                  //!<! Pt of the jet after subtracting average background
    float                        fEta;                     //!<! Eta of the jet
    float                        fPhi;                     //!<! Phi of the jet (0 < phi < 2pi)

    // Other jet quantities
    float                        fPtUncorr;                //!<! Pt of the jet (GeV/c) (not background subtracted)
    float                        fArea;                    //!<! Area of the jet
  
    // Jet substructure observables
    unsigned short int           fN;                       //!<! Number of jet constituents
    float                        fZLeading;                //!<! z of leading track
    float                        fRadialMoment;            //!<! Radial moment (not background subtracted)
    float                        fpTD;                     //!<! Momentum dispersion (not background subtracted)
    float                        fMass;                    //!<! Jet mass (not background subtracted)

    // Jet matching
    unsigned short int           fMatchedJetID;            //!<! JetID of matched jet in the matching jet tree
  
  /// \cond CLASSIMP
  ClassDef(AliJetTreeHandler,4); ///
  /// \endcond
};

#endif
