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
#include "AliEmcalTriggerPartBitConfig.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartBitConfig);
ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartBitConfigOld);
ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartBitConfigNew);

using namespace PWG::EMCAL::TriggerPart;

/**
 * Dummy constructor for the configuraiton base classes, not to be callled
 */
AliEmcalTriggerPartBitConfig::AliEmcalTriggerPartBitConfig():
TObject(),
fL0Bit(-1),
fJHighBit(-1),
fJLowBit(-1),
fGHighBit(-1),
fGLowBit(-1),
fTriggerTypesEnd(-1)
{
}

/**
 * Constructor initialising the configurations. Used by the inheriting classes
 * @param l0bit	Bit for Level0 (not used here)
 * @param jhighbit Bit for Jet High (EJ1)
 * @param jlowbit Bit for Jet Low (EJ2)
 * @param ghighbit Bit for Gamma High (EG1)
 * @param glowbit Bit for Gamma Low (EG2)
 * @param mcoffset Monte Carlo offset bit (not used here)
 */
AliEmcalTriggerPartBitConfig::AliEmcalTriggerPartBitConfig(
		int l0bit,
		int jhighbit,
		int jlowbit,
		int ghighbit,
		int glowbit,
		int mcoffset):
			TObject(),
			fL0Bit(l0bit),
			fJHighBit(jhighbit),
			fJLowBit(jlowbit),
			fGHighBit(ghighbit),
			fGLowBit(glowbit),
			fTriggerTypesEnd(mcoffset)
{
}

/**
 * Initialise from other object
 * @param ref Reference used to initialize this object
 */
void AliEmcalTriggerPartBitConfig::Initialise(const AliEmcalTriggerPartBitConfig& ref) {
	fL0Bit = ref.GetLevel0Bit();
	fJHighBit = ref.GetJetHighBit();
	fJLowBit = ref.GetJetLowBit();
	fGHighBit = ref.GetGammaHighBit();
	fGLowBit = ref.GetGammaLowBit();
	fTriggerTypesEnd = ref.GetTriggerTypesEnd();
}

/**
 * Settings for the 2-bit configuration
 */
AliEmcalTriggerPartBitConfigOld::AliEmcalTriggerPartBitConfigOld():
    		AliEmcalTriggerPartBitConfig(0,2,2,1,1,3)
{
}

/**
 * Settings for the 4-bit configuration
 */
AliEmcalTriggerPartBitConfigNew::AliEmcalTriggerPartBitConfigNew():
    		AliEmcalTriggerPartBitConfig(0,3,4,1,2,5)
{
}

