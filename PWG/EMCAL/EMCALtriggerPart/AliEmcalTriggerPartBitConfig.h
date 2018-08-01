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
#ifndef ALIEMCALTRIGGERPARTBITCONFIG_H
#define ALIEMCALTRIGGERPARTBITCONFIG_H

#include <exception>
#include <sstream>
#include <string>
#include <TObject.h>

namespace PWG {

namespace EMCAL {

namespace TriggerPart {

class AliEmcalTriggerPartBitConfig : public TObject {
public:
	class InvalidConfigurationException : public std::exception{
	public:
		InvalidConfigurationException(std::string bitname):
			fMessage("")
		{
			std::stringstream msgbuilder;
			msgbuilder << "Invalid trigger configuration: " << bitname << " bit < 0";
			fMessage = msgbuilder.str();
		}
		virtual ~InvalidConfigurationException() throw() {}

		const char *what() const throw() {
			return fMessage.c_str();
		}

	private:
		std::string 				fMessage;
	};

	AliEmcalTriggerPartBitConfig();
	AliEmcalTriggerPartBitConfig(int l0bit, int j1bit, int j2bit, int g1bit, int g2bit, int mcoffset);
	virtual ~AliEmcalTriggerPartBitConfig() {}

	void Initialise(const AliEmcalTriggerPartBitConfig &ref);

	int GetLevel0Bit() const { if(fL0Bit < 0) throw InvalidConfigurationException("Level0"); return fL0Bit; }
	int GetJetHighBit() const { if(fJHighBit < 0) throw InvalidConfigurationException("JetHigh"); return fJHighBit; }
	int GetJetLowBit() const { if(fJLowBit < 0) throw InvalidConfigurationException("JetLow"); return fJLowBit; }
	int GetGammaHighBit() const { if(fGHighBit < 0) throw InvalidConfigurationException("GammaHigh"); return fGHighBit; }
	int GetGammaLowBit() const { if(fGLowBit < 0) throw InvalidConfigurationException("GammaLow"); return fGLowBit; }
	int GetTriggerTypesEnd() const {if(fTriggerTypesEnd < 0) throw InvalidConfigurationException("MCOffset"); return fTriggerTypesEnd; }

protected:
	int fL0Bit;      			///< Level0 bit
	int fJHighBit;   			///< Jet High bit
	int fJLowBit;    			///< Jet Low bit
	int fGHighBit;   			///< Gamma High bit
	int fGLowBit;    			///< Gamma Low bit
	int fTriggerTypesEnd;   	///< Monte-Carlo offset

	ClassDef(AliEmcalTriggerPartBitConfig, 1);
};

class AliEmcalTriggerPartBitConfigOld : public AliEmcalTriggerPartBitConfig{
public:
	AliEmcalTriggerPartBitConfigOld();
	virtual ~AliEmcalTriggerPartBitConfigOld() {}

	ClassDef(AliEmcalTriggerPartBitConfigOld, 1);
};

class AliEmcalTriggerPartBitConfigNew : public AliEmcalTriggerPartBitConfig{
public:
	AliEmcalTriggerPartBitConfigNew();
	virtual ~AliEmcalTriggerPartBitConfigNew() {}

	ClassDef(AliEmcalTriggerPartBitConfigNew, 1);
};

}
}
}


#endif /* ALIEMCALTRIGGERPARTBITCONFIG_H */
