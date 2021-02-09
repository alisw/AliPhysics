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
#include <algorithm>
#include "AliEmcalTriggerPartChannelMap.h"
#include "AliEmcalTriggerPartJetAlgorithm.h"
#include "AliEmcalTriggerPartSetup.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartJetAlgorithm);

using namespace PWG::EMCAL::TriggerPart;

/**
 * Constructor
 */
AliEmcalTriggerPartJetAlgorithm::AliEmcalTriggerPartJetAlgorithm():
  AliEmcalTriggerPartAlgorithm()
{
}

/**
 * Destructor
 */
AliEmcalTriggerPartJetAlgorithm::~AliEmcalTriggerPartJetAlgorithm() {
}

/**
 * Gamma trigger algorithm
 * 1. Loop over all rows (- patchsize) to get the starting position of the patch
 * 2. Loop over ADC values in the 16x16 window
 * 3. Sorting of the trigger patches so that the highest energetic patch (main patch is the first)
 * 4. Fill the output trigger object
 * @param channes Input channel map
 * @return vector with trigger patches
 */
std::vector<AliEmcalTriggerPartRawPatch> AliEmcalTriggerPartJetAlgorithm::FindPatches(const AliEmcalTriggerPartChannelMap *channels) const {
	std::vector<AliEmcalTriggerPartRawPatch> rawpatches;

	double adcsum(0);
	for(unsigned char irow = 0; irow < channels->GetNumberOfRows() - 15; irow+=4){
		for(unsigned char icol = 0; icol < channels->GetNumberOfCols() - 15; icol+=4){
			// 16x16 window
			adcsum = 0;
			for(unsigned char jrow = 0; jrow < 16; jrow++)
				for(unsigned char jcol = 0; jcol < 16; jcol++)
					adcsum += channels->GetADC(icol + jcol, irow + jrow);

			// make decision, low and high threshold
			int triggerBits(0);
			if(adcsum > fTriggerSetup->GetThresholdJetHigh()) triggerBits |= 1 << fTriggerSetup->GetTriggerBitConfiguration().GetJetHighBit();
			if(adcsum > fTriggerSetup->GetThresholdJetLow()) triggerBits |= 1 << fTriggerSetup->GetTriggerBitConfiguration().GetJetLowBit();

			// Set special bit
			if(triggerBits){
				AliEmcalTriggerPartRawPatch patch(icol, irow, adcsum, triggerBits);
				patch.SetPatchSize(16);
				rawpatches.push_back(patch);
			}
		}
	}

	// sort patches so that the main patch appears first
	std::sort(rawpatches.begin(), rawpatches.end());
	return rawpatches;
}

std::vector<AliEmcalTriggerPartRawPatch> AliEmcalTriggerPartJetAlgorithm::FindPatches8x8(const AliEmcalTriggerPartChannelMap *channels) const {
	std::vector<AliEmcalTriggerPartRawPatch> rawpatches;

	double adcsum(0);
	for(unsigned char irow = 0; irow < channels->GetNumberOfRows() - 8-1; irow+=4){
		for(unsigned char icol = 0; icol < channels->GetNumberOfCols() - 8-1; icol+=4){
			// 8x8 window
			adcsum = 0;
			for(unsigned char jrow = 0; jrow < 8; jrow++)
				for(unsigned char jcol = 0; jcol < 8; jcol++)
					adcsum += channels->GetADC(icol + jcol, irow + jrow);

			// make decision, low and high threshold
			int triggerBits(0);
			if(adcsum > fTriggerSetup->GetThresholdJetHigh()) triggerBits |= 1 << fTriggerSetup->GetTriggerBitConfiguration().GetJetHighBit();
			if(adcsum > fTriggerSetup->GetThresholdJetLow()) triggerBits |= 1 << fTriggerSetup->GetTriggerBitConfiguration().GetJetLowBit();

			// Set special bit
			if(triggerBits){
				AliEmcalTriggerPartRawPatch patch(icol, irow, adcsum, triggerBits);
				patch.SetPatchSize(8);
				rawpatches.push_back(patch);
			}
		}
	}

	// sort patches so that the main patch appears first
	std::sort(rawpatches.begin(), rawpatches.end());
	return rawpatches;
}
