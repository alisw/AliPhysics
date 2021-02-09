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

/**
 * @class TestAliEmcalAODFilterBitCuts
 * @brief Unit test for AOD hybrid track cuts
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 18, 2017
 * 
 * Unit test class for AOD hybrid track cuts. Covering
 * - 2010 definition with non-refit tracks
 * - 2010 definition without non-refit tracks
 * - 2011 definition
 * In each test for each supported category one track is prepared as hybrid track of
 * a given category, and one track is prepared as non-hybrid track. For passing the test
 * in all cases the selection status must match (i.e. hybrid tracks true, non-hybrid tracks 
 * false). In addition for hybrid tracks a user object of type AliEmcalTrackSelResultHybrid
 * must be provided storing the correct track type. Non-hybrid track selection results must 
 * not carry a user object.
 */
class TestAliEmcalAODHybridTrackCuts : public TObject {
public:
  TestAliEmcalAODHybridTrackCuts();
  
  /**
   * @brief Destructor
   * 
   * Deleting track selection objects
   */
  virtual ~TestAliEmcalAODHybridTrackCuts();

  /**
   * @brief Initializing track selection objects
   */
  void Init();

  /**
   * @brief Run all unit tests for the class AliEmcalAODHybridTrackCuts
   * 
   * @return true All tests passed
   * @return false At least one failure observed
   */
  bool RunAllTests() const;

  /**
   * @brief Test for hybrid tracks according to the 2010 definition including non-refit tracks
   * 
   * Preparation of 4 AOD tracks
   * 1) Hybrid flag set, refit true, cat 1 track
   * 2) Hybrid flag set, refit true, cat 2 track
   * 3) Hybrid flag set, refit false, cat 2 track
   * 4) Hybrid flag not set
   * Track selection must determine correcty the selection status and for true hybrid tracks also the category must match
   * 
   * @return true  All tests passed
   * @return false At least one failure observed
   */
  bool TestDef2010wRefit() const;

  /**
   * @brief Test for hybrid tracks according to the 2010 definition excluding non-refit tracks
   * 
   * Preparation of 4 AOD tracks
   * 1) Hybrid flag set, refit true, cat 1 track
   * 2) Hybrid flag set, refit true, cat 2 track
   * 3) Hybrid flag set, refit false, cat 2 track
   * 4) Hybrid flag not set
   * Track selection must determine correcty the selection status (excluding non-hybrid track and hybrid track without refit) 
   * and for true hybrid tracks also the category must match
   * 
   * @return true  All tests passed
   * @return false At least one failure observed
   */
  bool TestDef2010woRefit() const;

  /**
   * @brief Test for hybrid tracks according to the 2011 definition
   * 
   * Preparation of 3 AOD tracks
   * 1) Hybrid flag set, cat 1 track
   * 2) Hybrid flag set, cat 2 track
   * 3) Hybrid flag not set
   * Track selection must determine correctly the selection status and for true hybrid tracks also the category must match
   * 
   * @return true  All tests passed
   * @return false At least one failure observed 
   */
  bool TestDef2011() const;

private:
  AliEmcalAODHybridTrackCuts    *fDef2010wRefit;        ///< Hybrid track definition from 2010 including non-refit tracks
  AliEmcalAODHybridTrackCuts    *fDef2010woRefit;       ///< Hybrid track definition from 2010 excluding non-refit tracks
  AliEmcalAODHybridTrackCuts    *fDef2011;              ///< Hybrid track definition from 2011 e

  TestAliEmcalAODHybridTrackCuts(const TestAliEmcalAODHybridTrackCuts &);
  TestAliEmcalAODHybridTrackCuts &operator=(const TestAliEmcalAODHybridTrackCuts &);

  /// \cond CLASSIMP
  ClassDef(TestAliEmcalAODHybridTrackCuts, 1);
  /// \endcond
};

}

}
#endif
