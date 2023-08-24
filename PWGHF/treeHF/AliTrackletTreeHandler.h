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

#ifndef ALITRACKLETTREEHANDLER_H
#define ALITRACKLETTREEHANDLER_H

/**
 * \class AliTrackletTreeHandler
 * \brief Helper class to handle a tree for cut optimisation and MVA analyses, based heavily on AliHFTreeHandler and AliAnalysisTaskDmesonJets
 *
 * \author Nima Zardoshti nima.zardoshti@cern.ch
 *         Benedikt Volkel benedikt.volkel@cern.ch
 *         Cristina Terrevoli cristina.terrevoli@cern.ch
 *         Gian Michele Innocenti ginnocen@cern.ch
 * \date June 27 2019
 */

#include "TTree.h"
#include "AliAODTracklets.h"
//________________________________________________________________
//****************************************************************
class AliTrackletTreeHandler : public TObject
{
  public:

    AliTrackletTreeHandler();
    virtual ~AliTrackletTreeHandler();

    // Core methods
    TTree* BuildTree(TString name, TString title);
    void FillTree(int runNumber, int eventID, int eventID_Ext, Long64_t eventID_Long);

    // Setters
    void SetTrackletContainer(AliAODTracklets* TrackletContainer) { fTrackletContainer = TrackletContainer; }

  protected:

    TTree*                       fTreeTracklet;            ///< Tree with tracklet variables

    AliAODTracklets*        fTrackletContainer;       //!<! Tracklet container for this tree

    // Track quantities.
    float                        fTrackletEta;             //!<! Eta of tracklet
    float                        fTrackletPhi;             //!<! Phi of tracklet (0 < phi < 2pi)

    // Event quantities
    int                          fRunNumber;               //!<! run number
    int                          fEventID;                 //!<! event ID (unique identifier when run number is fixed)
    int                          fEventIDExt;                 //!<! event ID (unique identifier when run number is fixed)
    Long64_t                     fEventIDLong;                 //!<! event ID (unique identifier when run number is fixed)

  /// \cond CLASSIMP
  ClassDef(AliTrackletTreeHandler,1); ///
  /// \endcond
};

#endif
