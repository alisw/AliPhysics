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
#ifndef ALIEMCALTRIGGERPARTALGORITHM_H
#define ALIEMCALTRIGGERPARTALGORITHM_H

#include <vector>
#include <TObject.h>

namespace PWG {
	
namespace EMCAL {

namespace TriggerPart {

class AliEmcalTriggerPartPatchContainer;
class AliEmcalTriggerPartChannelMap;
class AliEmcalTriggerPartSetup;

/**
 * @class AliEmcalTriggerPartRawPatch
 * @brief Helper structure for raw patches
 *
 * Stores minimal information of reconstructed patches, such as row, column, and
 * reconstructed amplitude
 */
class AliEmcalTriggerPartRawPatch : public TObject{
	/**
	 * Default constructor, initializing all values with -1
	 */
public:
	enum EPatchtype_t {
		kAny,
		kEMCALpatch,
		kDCALPHOSpatch,
		kEMCALpatchGA,
		kEMCALpatchJE,		
		kEMCALpatchJE8x8,		
		kDCALpatchGA,		
		kDCALpatchJE,				
		kDCALpatchJE8x8,				
		kUndefPatch
	};

	AliEmcalTriggerPartRawPatch():
		TObject(),
		fCol(0),
		fRow(0),
		fADC(-1.),
		fTriggerBits(0),
		fPatchSize(0),
		fPatchType(kUndefPatch)
	{}
	/**
	 * Constructor, initializing position and amplitude
	 * @param col Starting col
	 * @param row Starting row
	 * @param adc Patch amplitude
	 */
	AliEmcalTriggerPartRawPatch(unsigned char col, unsigned char row, double adc, unsigned int triggerBits):
		TObject(),
		fCol(col),
		fRow(row),
		fADC(adc),
		fTriggerBits(triggerBits),
		fPatchSize(0),
		fPatchType(kUndefPatch)
	{}
	/**
	 * Destructor
	 */
	virtual ~AliEmcalTriggerPartRawPatch() {}

	/**
	 * Comparison operator, comparing to other in terms of ADC value
	 * @param other Object to compare with
	 * @return True if this adc is smaller, false otherwise
	 */
	bool operator<(const AliEmcalTriggerPartRawPatch & other) const {
		return fADC < other.fADC;
	}

	/**
	 * Set the type of the patch (EMCAL or DCAL-PHOS)
	 * @param ptype Type of the patch
	 */
	void SetPatchType(EPatchtype_t ptype) { fPatchType = ptype; }

	/**
	 * Set the size of the patch
	 * @param patchsize
	 */
	void SetPatchSize(unsigned char patchsize) { fPatchSize = patchsize; }

	/**
	 * Get starting row of the patch
	 * @return starting row
	 */
	unsigned char GetRowStart() const { return fRow; }
	/**
	 * Get Starting column of the patch
	 * @return starting column
	 */
	unsigned char GetColStart() const { return fCol; }
	/**
	 * Get the patch amplitude
	 * @return the patch amplitude
	 */
	double GetADC() const { return fADC; }

	/**
	 * Get the size of the patch
	 * @return size of the patch
	 */
	unsigned char GetPatchSize() const { return fPatchSize; }

	/**
	 * Get the trigger bits
	 */
	int GetTriggerBits() const { return fTriggerBits; }

	/**
	 * Check whether patch is of type EMCAL
	 * @return True if patch is of type EMCAL, false otherwise
	 */
	bool IsEMCAL() const { return fPatchType == kEMCALpatch; }

	/**
	 * Check whether patch is of type DCAL-PHOS
	 * @return True if patch is of type DCAL-PHOS, false otherwise
	 */
	bool IsDCALPHOS() const { return fPatchType == kDCALPHOSpatch; }

	/**
	 * Get unique ID of the patch, calculated, from col, row and subregion size
	 * @return Unique ID of the patch;
	 */
	int GetID() const;

private:
	unsigned char        	fCol;           ///< Lower left column in col-row coordinate space
	unsigned char      		fRow;           ///< Lower left row in col-row coordinate space
	double      					fADC;           ///< ADC value of the raw patch
	unsigned int         	fTriggerBits;   ///< Tigger bit settings
	unsigned char					fPatchSize;	    ///< Patch size
	EPatchtype_t	  			fPatchType;		 	///< Type of the trigger patch

	ClassDef(AliEmcalTriggerPartRawPatch, 1);
};

/**
 * @class AliEmcalTriggerPartAlgorithm
 * @brief Base class for EMCAL trigger algorithms
 *
 * Base class for trigger algorithm implementations for the EMCAL trigger.
 */
class AliEmcalTriggerPartAlgorithm : public TObject {
public:
	AliEmcalTriggerPartAlgorithm();
	virtual ~AliEmcalTriggerPartAlgorithm() {}

	virtual std::vector<PWG::EMCAL::TriggerPart::AliEmcalTriggerPartRawPatch> FindPatches(const AliEmcalTriggerPartChannelMap * channels) const = 0;
	/**
	 * Set the trigger channel ADC map used to create the trigger patches
	 * @param inputdata input
	 */
	void SetTriggerSetup(AliEmcalTriggerPartSetup *triggersetup) { fTriggerSetup = triggersetup; }

protected:
	AliEmcalTriggerPartSetup  		              *fTriggerSetup;       ///< Trigger setup data

	ClassDef(AliEmcalTriggerPartAlgorithm, 1);
};

}

}

}

#endif /* ALIEMCALTRIGGERPARTALGORITHM_H */
