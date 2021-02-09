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
#ifndef ALIEMCALTRACKSELRESULTHYBRID_H
#define ALIEMCALTRACKSELRESULTHYBRID_H

#include <TObject.h>

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalTrackSelResultHybrid
 * @brief Selection result of the hybrid track selection
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 21, 2017
 * 
 * The object stores the type of the hybrid track determined by the AliEmcalESD/AODHybridTrackCuts.
 * For the various hybrid track types see @ref HybridType_t. The object is produced by the
 * hybrid track cuts object for each track and added as optional payload to the resulting
 * AliEmcalTrackSelResultPtr. Users can obtain the payload object from the AliEmcalTrackSelResultPtr
 * and obtain the hybrid track type by one of the various getters.
 */
class AliEmcalTrackSelResultHybrid : public TObject {
public:
  /**
   * @enum HybridType_t
   * @brief Various types of hybrid tracks
   * 
   * | CAT  | Hybrid track type         | Properties                                                                    |
   * |------|---------------------------|-------------------------------------------------------------------------------|
   * |  0   | kUndefined                | Not a hybrid track                                                            |
   * |  I   | kHybridGlobal             | Global hybrid track (with SPD any condition). No vertex constraint            |
   * |  IIa | kHybridConstrainedTrue    | No SPD cluster, no alive SPD module, probably primary. With vertex constraint |
   * |  IIb | kHybridConstrainedFake    | No SPD cluster but alive SPD module, probably secondary. No vertex constraint |
   * |  III | kHybridConstrainedNoRefit | Constrained track failing in refit, treated as primary particle               |
   */
  enum HybridType_t {
    kUndefined,                               ///< Undefined - no hybrid track
    kHybridGlobal,                            ///< Global hybrid track                                      (type I)
    kHybridConstrainedTrue,                   ///< Complementary hybrid track without alive SPD module      (type IIa)
    kHybridConstrainedFake,                   ///< Complementary hybrid track with alive SPD module         (type IIb)
    kHybridConstrainedNoITSrefit              ///< Complementary hybrid track without ITS refit             (type III)
  };

  /**
   * @brief Dumy constructor
   */

  AliEmcalTrackSelResultHybrid();

  /**
   * @brief Constructor with track type
   * 
   * For track types see @ref HybridType_t
   */
  AliEmcalTrackSelResultHybrid(HybridType_t tracktype);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTrackSelResultHybrid() {}

  /**
   * @brief Set the type of the hybrid track
   * 
   * For track types see @ref HybridType_t
   * @param tracktype Type of the hybrid track
   */
  void SetHybridTrackType(HybridType_t tracktype) { fHybridTrackType = tracktype; }

  /**
   * @brief Get the type of the hybrid track
   * 
   * For track types see @ref HybridType_t
   * @return Type of the hybrid track
   */
  HybridType_t GetHybridTrackType() const { return fHybridTrackType; }

  /**
   * @brief Check whether track is selected as hybrid track
   * @return True if the track is selected as hybrid track, false otherwise
   */
  Bool_t IsHybridTrack() const { return fHybridTrackType != kUndefined; }

  /**
   * @brief Check whether track is selected as global hybrid track
   * @return True if the track is selected as global hybrid track, false otherwise
   */
  Bool_t IsHybridTrackGlobal() const { return fHybridTrackType == kHybridGlobal; }

  /**
   * @brief Check whether track is selected as any complementary (true or fake) hybrid track
   * @return True if the track is selected as complementary hybrid track, false otherwise
   */
  Bool_t IsHybridTrackConstrained() const { return fHybridTrackType == kHybridConstrainedTrue || fHybridTrackType == kHybridConstrainedFake; }

  /**
   * @brief Check whether track is selected as true complementary hybrid track
   * 
   * For true complementary hybrid tracks both SPD modules were off
   * @return True if the track is selected as true complementary hybrid track, false otherwise
   */
  Bool_t IsHybridTrackConstrainedTrue() const { return fHybridTrackType == kHybridConstrainedTrue; } 

  /**
   * @brief Check whether track is selected as fake complementary hybrid track
   * 
   * For fake complementary hybrid tracks at least one SPD module was active
   * @return True if the track is selected as fake complementary hybrid track, false otherwise
   */
  Bool_t IsHybridTrackConstrainedFake() const { return fHybridTrackType == kHybridConstrainedFake; }

  /**
   * @brief Check whether track is selected as non-refit complementary hybrid track
   * @return True if the track is selected as non-refit complementary hybrid track, false otherwise
   */
  Bool_t IsHybridTrackNonRefit() const { return fHybridTrackType == kHybridConstrainedNoITSrefit; }

private:
  HybridType_t                  fHybridTrackType;           ///< Hybrid track type

  /// \cond CLASSIMP
  ClassDef(AliEmcalTrackSelResultHybrid, 1);
  /// \endcond
};

}
}

#endif