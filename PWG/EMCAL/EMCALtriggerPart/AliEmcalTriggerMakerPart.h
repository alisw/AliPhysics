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
#ifndef ALIEMCALTRIGGERMAKERPART_H
#define ALIEMCALTRIGGERMAKERPART_H

#include <TObject.h>

#include "AliEmcalTriggerPartGammaAlgorithm.h"
#include "AliEmcalTriggerPartJetAlgorithm.h"
#include "AliEmcalTriggerPartChannelMap.h"
#include "AliEmcalTriggerPartBadChannelContainer.h"
#include "AliEmcalTriggerPartMapping.h"
#include "AliEmcalTriggerPartSetup.h"

namespace PWG {

namespace EMCAL {

namespace TriggerPart {

class AliEmcalTriggerMakerPart : public TObject {
public:
	enum {
		kMinRowPHOS = 0,
		kMaxRowPHOS = 35,
		kMinEtaPHOS = 16,
		kMaxEtaPHOS = 31
	};
	
	AliEmcalTriggerMakerPart();
	virtual ~AliEmcalTriggerMakerPart();

	void 																					Reset();
	void 																					FindPatches();
	std::vector<AliEmcalTriggerPartRawPatch>			GetPatches(const int what = AliEmcalTriggerPartRawPatch::kAny);
	AliEmcalTriggerPartRawPatch 									GetMaxGammaEMCAL();
	AliEmcalTriggerPartRawPatch 									GetMaxGammaDCALPHOS();
	AliEmcalTriggerPartRawPatch 									GetMaxJetEMCAL();
	AliEmcalTriggerPartRawPatch 									GetMaxJetDCALPHOS();
	AliEmcalTriggerPartRawPatch 									GetMaxJetEMCAL8x8();
	AliEmcalTriggerPartRawPatch 									GetMaxJetDCALPHOS8x8();

	double 																				GetMedian(std::vector<AliEmcalTriggerPartRawPatch> v);

	double 																				GetMedianGammaEMCAL();
	double 																				GetMedianGammaDCALPHOS();
	double 																				GetMedianJetEMCAL();
	double 																				GetMedianJetDCALPHOS();
	double 																				GetMedianJetEMCAL8x8();
	double 																				GetMedianJetDCALPHOS8x8();

	void FillChannelMap(double eta, double phi, double energy);

	/**
	 * Get the map of the EMCAL trigger channels
	 * @return Map of the EMCAL trigger channels
	 */
	const AliEmcalTriggerPartChannelMap &GetEMCALChannels() const { return fTriggerChannelsEMCAL; }

	/**
	 * Get the map of the EMCAL trigger channels
	 * @return Map of the EMCAL trigger channels
	 */
	const AliEmcalTriggerPartChannelMap &GetDCALPHOSChannels() const { return fTriggerChannelsDCALPHOS; }

	/**
	 * Get the mapping between eta and phi on the one side and row and col in the EMCAL / DCAL on the other side
	 * @return Mapping for EMCAL and DCAL/PHOS trigger channels
	 */
	const AliEmcalTriggerPartMapping &GetTriggerChannelMapping() const { return fTriggerMapping; }

											/**
	 					* Setup trigger patch finders
	 * @param setup Configuration of 					the trigger patch finders									 */	
	void SetTriggerSetup(AliEmcalTriggerPartSetup	 &setup) { fTriggerSetup = setup; }

	/*		*
	 * Accept pa														tches 100% in PHOS
	 * @param d	o													Accept switch whether we accept or not
	 */
	void SetAcceptPHOSPatches(bool doAccept) { fAcceptPHOSPatches = doAccept; }

	/**
	 * Add bad channel position in EMCAL in row and col to the list of bad channels.
	 * Bad channels will be ignored when setting the energy
	 * @param col Column of the bad channel
	 * @param row Row of the bad channel
	 */
	void AddBadChannelEMCAL(int col, int row) { fBadChannelsEMCAL.AddChannel(col, row); }

	/**
	 * Add bad channel position in EMCAL in row and col to the list of bad channels.
	 * Bad channels will be ignored when setting the energy
	 * @param col Column of the bad channel
	 * @param row Row of the bad channel
	 */
	void AddBadChannelDCALPHOS(int col, int row) { fBadChannelsDCALPHOS.AddChannel(col, row); }

	/**
	 * Get bad channels container
	 */
	AliEmcalTriggerPartBadChannelContainer GetBadChannelContainerEMCAL() const
	{
		return fBadChannelsEMCAL;
	}

	/**
	 * Get bad channels container
	 */
	AliEmcalTriggerPartBadChannelContainer GetBadChannelContainerDCALPHOS() const
	{
		return fBadChannelsDCALPHOS;
	}

	/**
	 * Check whether patch is in PHOS
	 * @param col
	 * @param row
	 * @param size
	 * @return
	 */
	bool IsPHOSPatch(int col, int row, int size);

private:
	AliEmcalTriggerPartJetAlgorithm								fJetTrigger;								///< Algorithm finding jet patches on a trigger channel map
	AliEmcalTriggerPartGammaAlgorithm							fGammaTrigger;							///< Algorithm finding gamma patches on a trigger channel map
	AliEmcalTriggerPartChannelMap									fTriggerChannelsEMCAL;			///< Trigger channels for the EMCAL
	AliEmcalTriggerPartChannelMap									fTriggerChannelsDCALPHOS;		///< Trigger channels for the combination DCAL-PHOS
	AliEmcalTriggerPartMapping										fTriggerMapping;						///< Mapping between trigger channels and eta and phi
	AliEmcalTriggerPartSetup											fTriggerSetup;							///< Setup of the EMCAL / DCAL-PHOS trigger algorithms
	AliEmcalTriggerPartBadChannelContainer				fBadChannelsEMCAL;					///< Map with bad EMCAL channels
	AliEmcalTriggerPartBadChannelContainer				fBadChannelsDCALPHOS;				///< Map with bad DCAL-PHOS channels
	bool																					fHasRun;
	bool 																					fAcceptPHOSPatches;					///< Accept patches 100% in PHOS
	std::vector<AliEmcalTriggerPartRawPatch>			fGammaEMCAL;
	std::vector<AliEmcalTriggerPartRawPatch>			fGammaDCALPHOS;
	std::vector<AliEmcalTriggerPartRawPatch>			fJetEMCAL;
	std::vector<AliEmcalTriggerPartRawPatch>			fJetDCALPHOS;
	std::vector<AliEmcalTriggerPartRawPatch>			fJetEMCAL8x8;
	std::vector<AliEmcalTriggerPartRawPatch>			fJetDCALPHOS8x8;

	ClassDef(AliEmcalTriggerMakerPart, 1);
};

}
}
}


#endif /* SlLImMCALRIGGERMAKERPaARTH */
