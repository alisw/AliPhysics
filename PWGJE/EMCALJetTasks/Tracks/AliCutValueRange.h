#ifndef ALICUTVALUERANGE_H
#define ALICUTVALUERANGE_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iosfwd>
#include <TObject.h>

namespace PWGJE {

namespace EMCALJetTasks{

template<typename t>
class AliCutValueRange;

template<typename t>
std::ostream &operator<<(std::ostream &stream, const AliCutValueRange<t> &val);

/**
 * @class AliCutValueRange
 * @brief Class containing a range for a value to cut on
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since Dec 12, 2014
 *
 * Class defining a range in which a value to be checked is valid. Can be used
 * as a cut. In case a negative comparison (value valid only outside this range)
 * is desired, this is handled when setting the object to negate (function Negate()).
 * The class is a template, expecting the comparison operators to be overloaded.
 *
 */
template<typename t>
class AliCutValueRange : public TObject {
public:
	/**
	 * Dummy constructor, producing a range open to both sides
	 */
	AliCutValueRange();
	/**
	 * Constructor, producing a range closed to both sides
	 *
	 * @param[in] min lower limit
	 * @param[in] max upper limit
	 */
	AliCutValueRange(t min, t max);
	/**
	 * Constructor, producing a range closed to both sides
	 *
	 * @param[in] limit the limit to be set
	 * @param[in] isUpper defining whether the limit is the upper (case true) or lower limit
	 */
	AliCutValueRange(t limit, bool isUpper);
	/**
	 * Destructor, nothing to do
	 */
	virtual ~AliCutValueRange() {}

	friend std::ostream &operator<< <>(std::ostream &stream, const AliCutValueRange<t> &val);

	void SetLimits(t min, t max){
		fLimits[0] = min;
		fLimits[1] = max;
		fHasLimit[0] = fHasLimit[1] = true;
	}
	void UnsetLimits(){ fHasLimit[0] = fHasLimit[1] = false; }
	/**
	 * Set a limit
	 * @param[in] value Value of the limit
	 * @param[in] isUpper Indicator for upper (true) or lower (false) limit
	 */
	void SetLimit(t value, bool isUpper){
		int bin = isUpper ? 1 : 0;
		fLimits[bin] = value;
		fHasLimit[bin] = true;
	}
	/**
	 * Make cut range open in one direction
	 * @param[in] isUpper Indictator for upper (true) or lower (false)
	 */
	void UnsetLimit(bool isUpper){
		int bin = isUpper ? 1 : 0;
		fHasLimit[bin] = false;
	}
	/**
	 * Define that the decision should be negated
	 */
	void Negate() { fNegate = true; }
	/**
	 * Define that decision should not be negated
	 */
	void SetPositive() { fNegate = false; }
	/**
	 * Define whether we use larger equal instead of larger
	 * @param[in] doUse if true we use larger equal, otherwise larger
	 */
	void SetUseLargerEqual(Bool_t doUse = true) { fUseLargerEqual = doUse; }
	/**
	 * Define whether we use smaller equal instead of smaller
	 * @param[in] doUse if true we use smaller equal, otherwise smaller
	 */
	void SetUseSmallerEqual(Bool_t doUse = true) { fUseSmallerEqual = doUse; }
	/**
	 * Check whether value is within a given range
	 *
	 * @param[in] value value to be checked
	 * @return comparison result
	 */
	bool IsInRange(t value) const;
  /**
   * Get the minimum value of the cut range
   * @return minimum value of the cut range
   */
	t GetMinimum() const { return fLimits[0]; }
	/**
	 * Get the maximum value of the cut range
	 * @return maximum value of the cut range
	 */
	t GetMaximum() const { return fLimits[1]; }

	/**
	 * Print cut values to the ostream
	 * @param[in] stream stream used for printout
   */
	void PrintStream(std::ostream &stream) const;

private:
	t       fLimits[2];                 ///< Specifies the limit in either of the direction (not used unless fHasLimit of that direction is true)
	bool    fHasLimit[2];               ///< Specifies whether limit in any of the two directions is set
	bool    fNegate;                    ///< Defines whether result should be inverted
	bool    fUseSmallerEqual;           ///< Use smaller equal for upper bound (true by default)
	bool    fUseLargerEqual;            ///< Use larger equal for lower bound (true by default)

	ClassDef(AliCutValueRange, 1);     // Value range for cuts
};

}

}

#endif
