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
#ifndef ALIESDTRACKCUTSWRAPPER_H
#define ALIESDTRACKCUTSWRAPPER_H

#include "AliVCuts.h"

class AliESDtrackCuts;

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalESDTrackCutsWrapper
 * @brief Wrapper class around AliESDtrackCuts object with focus on AOD analysis
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 6, 2017
 * 
 * The function IsSelected is not properly implemented for AliAODTracks in AliESDtrackCuts,
 * though a minimal functionality is provided. Instead AliESDtrackCuts offer the functions
 * AcceptTrack for AliESDtracks, automatically used in IsSelected of AliESDtrackCuts, and 
 * AcceptVTrack for other track objects inheriting from AliVTrack. The class delegates the 
 * IsSelect function to the proper Select function for the input data type. The focus is 
 * on virtual track selection frameworks which abstractize the track selection in order to
 * be blind to the input data type.
 */
class AliEmcalESDtrackCutsWrapper : public AliVCuts {
public:

  /**
   * @brief Dummy constructor
   */
  AliEmcalESDtrackCutsWrapper();

  /**
   * @brief Main constructor
   * 
   * @param name Name of the track cuts
   * @param trackcuts Underlying track cuts object
   */
  AliEmcalESDtrackCutsWrapper(const char *name, AliESDtrackCuts *trackcuts = NULL);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalESDtrackCutsWrapper();

  void SetESDTrackCuts(AliESDtrackCuts *trackcuts) { fTrackCuts = trackcuts; }
  AliESDtrackCuts *GetTrackCuts() const { return fTrackCuts; }

  /**
   * @brief Checking whether track is accepted by the underlying AliESDtrackCuts
   * 
   * @param o Object to be checked
   * @return true Track is accepted by the track cuts
   * @return false Track is rejected by the track cuts or input type is unsupported
   */
  virtual bool IsSelected(TObject *o);

private:

  AliESDtrackCuts         *fTrackCuts;      ///< Underlying track cuts object

  /// \cond CLASSIMP
  ClassDef(AliEmcalESDtrackCutsWrapper, 1);
  /// \endcond
};

}

}
#endif
