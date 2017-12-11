#ifndef ALIEMCALTRACKSELECTIONAOD_H_
#define ALIEMCALTRACKSELECTIONAOD_H_
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliEmcalTrackSelection.h>
#include "AliESDtrackCuts.h"
#include "AliEmcalTrackSelResultPtr.h"

class AliVCuts;
class AliVTrack;

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
	void AddFilterBit(UInt_t filterbits) { fFilterBits |= filterbits; }

	/**
	 * @brief Returns the hybrid filter bits according to a hard-coded look-up table
	 * @param[in] bits c-array of 2 elements where the bits are returned
	 * @param[in] period data taking period
	 * @return true if successful, false otherwise
	 */
	static Bool_t GetHybridFilterBits(Char_t bits[], TString period);

private:
	UInt_t			fFilterBits;				    ///< Track filter bits
	Bool_t      fFilterHybridTracks;    ///< Filter hybrid tracks using AliAODTrack::IsHybridGlobalConstrainedGlobal
	Bool_t      fFilterTPCTracks;       ///< Filter TPC-only tracks using AliAODTrack::IsHybridGlobalConstrainedGlobal
	Char_t      fHybridFilterBits[2];   ///< Filter bits of hybrid tracks

	/// \cond CLASSIMP
	ClassDef(AliEmcalTrackSelectionAOD, 2);
	/// \endcond
};

#endif /* ALIEMCALTRACKSELECTIONAOD_H_ */
