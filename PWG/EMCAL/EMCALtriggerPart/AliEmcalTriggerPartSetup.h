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
#ifndef ALIEMCALTRIGGERPARTSETUP_H
#define ALIEMCALTRIGGERPARTSETUP_H

#include <TObject.h>
#include "AliEmcalTriggerPartBitConfig.h"

namespace PWG {

namespace EMCAL {

namespace TriggerPart {

/**
 * @class AliEmcalTriggerPartSetup
 */
class AliEmcalTriggerPartSetup : public TObject {
public:
	AliEmcalTriggerPartSetup();
	AliEmcalTriggerPartSetup(const AliEmcalTriggerPartSetup &p);
	AliEmcalTriggerPartSetup &operator=(const AliEmcalTriggerPartSetup &p);
	/**
	 * Destructor
	 */
	virtual ~AliEmcalTriggerPartSetup() {}

	double GetThresholdJetLow() const { return fThresholds[2]; }
	double GetThresholdJetHigh() const { return fThresholds[0]; }

	double GetThresholdGammaLow() const { return fThresholds[3]; }
	double GetThresholdGammaHigh() const { return fThresholds[1]; }

	const AliEmcalTriggerPartBitConfig &GetTriggerBitConfiguration() const { return fTriggerBitConfig; }

	void SetThresholds( double i0, double i1, double i2, double i3 ) {
		fThresholds[0] = i0; fThresholds[1] = i1; fThresholds[2] = i2; fThresholds[3] = i3;}
	void SetTriggerBitConfig(const AliEmcalTriggerPartBitConfig &ref){
		fTriggerBitConfig.Initialise(ref);
	}

	void Clean();

protected:
	double                                  fThresholds[4];                 ///< per event L1 online thresholds in ADC counts
	AliEmcalTriggerPartBitConfig         		fTriggerBitConfig;              ///< Trigger bit configuration

	ClassDef(AliEmcalTriggerPartSetup, 1);
};

}
}
}

#endif
