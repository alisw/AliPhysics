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
#ifndef ALIEMCALTASKTRACKSELECTIONESD_H_
#define ALIEMCALTASKTRACKSELECTIONESD_H_

#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelResultPtr.h"

class TList;
class AliVCuts;
class AliVTrack;

/**
 * @class AliEmcalTrackSelectionESD
 * @brief Implementation of virtual track selection for ESDs
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @date Jul 24, 2015
 *
 * Implementation of the track selection for the analysis on ESDs using
 * AliESDtrackCuts as underlying structure
 */
class AliEmcalTrackSelectionESD: public AliEmcalTrackSelection {
public:

  /**
   * @brief Dummy constructor
   */
	AliEmcalTrackSelectionESD();

	/**
	 * @brief Constructor with cuts
	 */
	AliEmcalTrackSelectionESD(AliVCuts *cuts);

	/**
	 * @brief Constructor, initalising track cuts depending on the requested type of filtering
	 *
	 * @param[in] type Track filtering type
	 * @param[in] period  Period string (e.g. LHC11h)
	 */
	AliEmcalTrackSelectionESD(ETrackFilterType_t type, const char* period = "");

	/**
	 * @brief Destructor
	 *
	 * Cleaning up memory. AliESDtrackCuts objects which are
	 * stored in the QA output are not handled in the destructor
	 * as ownership changed.
	 */
	virtual ~AliEmcalTrackSelectionESD() {}

	/**
	 * @brief Automatically generates track cuts depending on the requested type of filtering
	 *
	 * @param[in] type    Track filtering type
	 * @param[in] period  Period string (e.g. LHC11h)
	 */
	virtual void GenerateTrackCuts(ETrackFilterType_t type, const char* period = "");

	/**
	 * @brief Check whether track is accepted.
	 *
	 * Iterates over all cuts assigned to the track selection.
	 * @param[in] trk Track to check
	 * @return true if selected, false otherwise
	 */
	virtual PWG::EMCAL::AliEmcalTrackSelResultPtr IsTrackAccepted(AliVTrack * const trk);

  virtual void SaveQAObjects(TList *outputList);

	/// \cond CLASSIMP
	ClassDef(AliEmcalTrackSelectionESD,1);
	/// \endcond
};

#endif /* ALIEMCALPTTASKTRACKSELECTIONESD_H_ */
