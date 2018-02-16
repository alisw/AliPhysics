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
#ifndef ALIEMCALTRACKSELECTIONAOD_H_
#define ALIEMCALTRACKSELECTIONAOD_H_

#include <AliEmcalTrackSelection.h>
#include "AliEmcalTrackSelResultPtr.h"

class AliVCuts;
class AliVTrack;

class AliEmcalTrackSelResultHybrid;

/**
 * @class AliEmcalTrackSelectionAOD
 * @brief Implement virtual track selection for AOD analysis
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @date Jul 24, 2015
 *
 * Implementation of track selection in case the analysis runs on AODs
 * For the moment it uses the AliESDtrackCuts and converts AOD tracks to
 * ESD tracks, which might change in the future when an AOD track selection
 * framework becomes available.
 */
class AliEmcalTrackSelectionAOD: public AliEmcalTrackSelection {
public:

  /**
   * @brief Main constructor

   * Initializes fields with 0 (or NULL). For ROOT I/O, not intended
   * to be used by the users.
   */
	AliEmcalTrackSelectionAOD();

	/**
	 * @brief Constructor for periods

	 * Initializing track cuts depending on the requested type of filtering
	 * @param[in] type Track filtering type
	 * @param[in] period  Period string (e.g. LHC11h)
	 */
	AliEmcalTrackSelectionAOD(ETrackFilterType_t type, const char* period = "");

	/**
   * @brief Main Constructor

   * Initalizing also track cuts and filter bits. In case the initial cuts
   * is a nullpointer, only filter bits are used for the track selection. This constructor is
   * intended to be used by the users.
   *
   * @param[in] cuts Inital track cut object (of type AliESDtrackCuts, can be a nullpointer)
   * @param[in] filterbits Filter bits required
   */
	AliEmcalTrackSelectionAOD(AliVCuts *cuts, UInt_t filterbits);

	/**
	 * @brief Destructor
	 */
	virtual ~AliEmcalTrackSelectionAOD() {}

	/**
	 * @brief Automatically generates track cuts depending on the requested type of filtering
	 * @param[in] type Track filtering type
	 */
	virtual void GenerateTrackCuts(ETrackFilterType_t type, const char* /*period*/ = "");

	/**
	 * @brief Performing track selection
	 *
	 * Function checks whether track is accepted under the given track selection cuts.
	 * The function can handle AliAODTrack and AliPicoTrack, while for AliPico track an
	 * AliAODTrack is expected to be the underlying structure. If it is not possible to
	 * access an AOD track from the input track, the object will not be selected. Otherwise
	 * first the status bits are checked (if requested), and if further track cuts (of type
	 * AliESDtrackCuts) are provided, the track is converted to an ESD track for further checks.
	 *
	 * @param[in] trk Track to check
	 * @return true if selected, false otherwise
	 */
	virtual PWG::EMCAL::AliEmcalTrackSelResultPtr IsTrackAccepted(AliVTrack * const trk);

	/**
	 * @brief Add a new filter bit to the track selection.
	 *
	 * Multiple filter bits can be set at the same time
	 * (via the bitwise or operator |).
	 *
	 * @param filterbits Filter bits used to select tracks
	 */
	void AddFilterBit(UInt_t filterbits);

	/**
	 * @brief Returns the hybrid filter bits according to a hard-coded look-up table
	 * @param[in] bits c-array of 2 elements where the bits are returned
	 * @param[in] period data taking period
	 * @return true if successful, false otherwise
	 */
	static Bool_t GetHybridFilterBits(Char_t bits[], TString period);

private:

	/// \cond CLASSIMP
	ClassDef(AliEmcalTrackSelectionAOD, 2);
	/// \endcond
};

namespace PWG {

namespace EMCAL{

/**
 * @class TestAliEmcalTrackSelectionAOD
 * @brief Unit test for the class AliEmcalTrackSelectionAOD
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 19, 2018 
 */
class TestAliEmcalTrackSelectionAOD : public TObject {
public:

	/**
	 * @brief Constructor
	 */
	TestAliEmcalTrackSelectionAOD();

	/**
	 * @brief Destructor
	 */
	virtual ~TestAliEmcalTrackSelectionAOD();

	/**
	 * @brief Init test suite
	 * 
	 * Create track selection objects for the various selections supported in the test suite
	 */
	void Init();

	/**
	 * @brief Run all tests
	 * 
	 * @return true  All tests passed
	 * @return false At least one test failed
	 */
	bool RunAllTests() const;
	bool TestHybridDef2010wRefit() const;
	bool TestHybridDef2010woRefit() const;
	bool TestHybridDef2011() const;
	bool TestTPConly() const;

private:
	/**
	 * @brief Extract hybrid track user object from a track selection result ptr
	 * 
	 * Tool used to extract recursively a user object of type AliEmcalTrackSelResultHybrid. 
	 * In case the user object is of type AliEmcalTrackSelResultCombined it tries to find
	 * the hybrid user information within the results of the combined track selection result.
	 * 
	 * @param data Track selection result pointer to be inspected
	 * @return Pointer to the hybrid track selection user object (if existing), nullptr otherwise
	 */
	const AliEmcalTrackSelResultHybrid 			*FindHybridSelectionResult(const AliEmcalTrackSelResultPtr &data) const;

	AliEmcalTrackSelectionAOD 				*fTrackSelHybrid2010wRefit;				///< Hybrid tracks from 2010 including non-refit tracks
	AliEmcalTrackSelectionAOD					*fTrackSelHybrid2010woRefit;			///< Hybrid tracks from 2010 excluding non-refit tracks
	AliEmcalTrackSelectionAOD					*fTrackSelHybrid2011;							///< Hybrid tracks from 2011
	AliEmcalTrackSelectionAOD					*fTrackSelTPConly;								///< TPConly tracks

	TestAliEmcalTrackSelectionAOD(const TestAliEmcalTrackSelectionAOD &);
	TestAliEmcalTrackSelectionAOD &operator=(const TestAliEmcalTrackSelectionAOD &);

	/// \cond CLASSIMP
	ClassDef(TestAliEmcalTrackSelectionAOD, 1);
	/// \endcond
};

}

}

#endif /* ALIEMCALTRACKSELECTIONAOD_H_ */
