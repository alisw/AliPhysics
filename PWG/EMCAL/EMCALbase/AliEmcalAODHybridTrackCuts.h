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
#ifndef ALIEMCALAODHYBRIDTRACKCUTS_H
#define ALIEMCALAODHYBRIDTRACKCUTS_H

#include "AliEmcalCutBase.h"

namespace PWG {

namespace EMCAL {
    
/**
 * @class AliEmcalAODHybridTrackCuts
 * @brief Cut class selecting hybrid tracks using the IsHybrid function
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 6, 2017
 * 
 * Selection of hybrid tracks is expressed as a cut class inheriting from
 * AliVCuts. This functionality is needed in virtual track selections in 
 * order to overcome special treatment of AOD tracks. Internally the class
 * fully relies on the function IsHybridTrackGlobalConstrainedGlobal from
 * AliAODTrack.
 */
 class AliEmcalAODHybridTrackCuts : public AliEmcalCutBase {
 public:

  /**
   * @brief Dummy constructor
   */
  AliEmcalAODHybridTrackCuts();

  /**
   * @brief Main constructor
   * 
   * @param name Name of the hybrid track cuts
   */
  AliEmcalAODHybridTrackCuts(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalAODHybridTrackCuts() {}

  /**
   * @brief Run track selection of hybrid tracks
   * 
   * @param o Object (AliAODTrack) to be tested
   * @return Track selection result with the selection status and the hybrid track type
   */
  virtual AliEmcalTrackSelResultPtr IsSelected(TObject *o);

  /**
   * @brief Switch on/off selection of hybrid tracks without ITSrefit
   * 
   * Only valid for productions which use hybrid track definitions including 
   * tracks without ITSrefit
   * 
   * @param doReject If true hybrid tracks without ITSrefit are rejected
   */
  void SetSelectNonITSrefitTracks(bool doReject) { fSelectNonITSrefitTracks = doReject; }

  /**
   * @brief Set the filterbits used to distinguish the different hybrid track types
   * 
   * @param globalfilterbit Filterbit for global hybrid tracks
   * @param constrainedfilterbit Filterbit for constrained hybrid tracks (+ non-refit hybrid tracks if available)
   */
  void SetHybridFilterBits(Int_t globalfilterbit, Int_t constrainedfilterbit) {
    fHybridFilterBits[0] = globalfilterbit;
    fHybridFilterBits[1] = constrainedfilterbit;
  }

private:
  Bool_t                            fSelectNonITSrefitTracks;  ///< Select non-refit tracks
  Int_t                             fHybridFilterBits[2];      ///< Bit numbers for various hybrid filter bits

  /// \cond CLASSIMP
  ClassDef(AliEmcalAODHybridTrackCuts, 1);
  /// \endcond
  
};

}

}
#endif
