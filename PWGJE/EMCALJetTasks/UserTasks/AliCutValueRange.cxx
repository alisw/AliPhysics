/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
 * Class defining a range in which a value to be checked is valid. Can be used
 * as a cut. In case a negative comparison (value valid only outside this range)
 * is desired, this is handled when setting the object to negate (function Negate()).
 * The class is a template, expecting the comparison operators to be overloaded.
 *
 *   Author: Markus Fasel
 */

#include "AliCutValueRange.h"

templateClassImp(EMCalTriggerPtAnalysis::AliCutValueRange)

namespace EMCalTriggerPtAnalysis {

	//______________________________________________________________________________
	template<typename t>
	AliCutValueRange<t>::AliCutValueRange():
	fNegate(false)
	{
		/*
		 * Dummy constructor, producing a range open to both sides
		 */
		fHasLimit[0] = fHasLimit[1] = false;
	}

	//______________________________________________________________________________
	template<typename t>
	AliCutValueRange<t>::AliCutValueRange(t min, t max):
	fNegate(false)
	{
		/*
		 * Constructor, producing a range closed to both sides
		 *
		 * @param min: lower limit
		 * @param max: upper limit
		 */
		fLimits[0] = min;
		fLimits[1] = max;
		fHasLimit[0] = fHasLimit[1] = true;
	}

	//______________________________________________________________________________
	template<typename t>
	AliCutValueRange<t>::AliCutValueRange(t limit, bool isUpper):
	fNegate(false)
	{
		/*
		 * Constructor, producing a range closed to both sides
		 *
		 * @param limit: the limit to be set
		 * @param isUpper: defining whether the limit is the upper (case true) or lower limit
		 */
		if(isUpper){
			fLimits[1] = limit;
			fHasLimit[0] = false;
			fHasLimit[1] = true;
		} else {
			fLimits[0] = limit;
			fHasLimit[0] = true;
			fHasLimit[1] = false;
		}
	}

	//______________________________________________________________________________
	template<typename t>
	bool AliCutValueRange<t>::IsInRange(t value) const {
		/*
		 * Check whether value is within a given range
		 *
		 * @param value: value to be checked
		 * @return: comparison result
		 */
		bool result = true;
		if(fHasLimit[0] && fHasLimit[1]){
			// Double-sided limited
			result = fNegate ? (value < fLimits[0] || value > fLimits[1]) : (value > fLimits[0] && value < fLimits[1]);
		} else if(fHasLimit[1]) {
			// only upper bound
			result = fNegate ? (value > fLimits[1]) : (value < fLimits[1]);
		} else if(fHasLimit[0]){
			// only lower bound
			result = fNegate ? (value < fLimits[0]) : (value > fLimits[0]);
		}
		return result;
	}

	template class AliCutValueRange<int>;
	template class AliCutValueRange<double>;
	template class AliCutValueRange<float>;

}

