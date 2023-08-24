/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include "AliCutValueRange.h"

templateClassImp(PWGJE::EMCALJetTasks::AliCutValueRange)

namespace PWGJE {
  namespace EMCALJetTasks {

  template<typename t>
  AliCutValueRange<t>::AliCutValueRange():
  fNegate(false),
  fUseSmallerEqual(true),
  fUseLargerEqual(true)
  {
    fHasLimit[0] = fHasLimit[1] = false;
  }

  template<typename t>
  AliCutValueRange<t>::AliCutValueRange(t min, t max):
  fNegate(false),
  fUseSmallerEqual(true),
  fUseLargerEqual(true)
  {
    fLimits[0] = min;
    fLimits[1] = max;
    fHasLimit[0] = fHasLimit[1] = true;
  }

  template<typename t>
  AliCutValueRange<t>::AliCutValueRange(t limit, bool isUpper):
  fNegate(false),
  fUseSmallerEqual(true),
  fUseLargerEqual(true)
  {
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

  template<typename t>
  bool AliCutValueRange<t>::IsInRange(t value) const {
    bool result = true;
    if(fHasLimit[0] && fHasLimit[1]){
      // Double-sided limited, correct choise of comparison operator
      Bool_t withinUpper = fUseSmallerEqual ? (value <= fLimits[1]) : (value < fLimits[1]),
          withinLower = fUseLargerEqual ? (value >= fLimits[0]) : (value > fLimits[0]);
      result = withinLower && withinUpper;
      if(fNegate)
        result = !result;
    } else if(fHasLimit[1]) {
      // only upper bound, correct choise of comparison operator
      result = fUseSmallerEqual ? (value <= fLimits[1]) : (value < fLimits[1]);
      if(fNegate)
        result = !result;
    } else if(fHasLimit[0]){
      // only lower bound, correct choise of comparison operator
      result = fUseLargerEqual ? (value >= fLimits[0]) : (value > fLimits[0]);
      if(fNegate)
        result = !result;
    }
    return result;
  }

  template<typename t>
  void AliCutValueRange<t>::PrintStream(std::ostream &stream) const {
    stream << "[" << fLimits[0] << "|" << fLimits[1] << "]";
  }

  template<typename t>
  std::ostream &operator<<(std::ostream &stream, const AliCutValueRange<t> &val){
    val.PrintStream(stream);
    return stream;
  }

  /*
   * Template specifications
   */
	template class AliCutValueRange<int>;
	template class AliCutValueRange<double>;
	template class AliCutValueRange<float>;
	template std::ostream &operator<<(std::ostream &stream, const AliCutValueRange<int> &val);
	template std::ostream &operator<<(std::ostream &stream, const AliCutValueRange<double> &val);
	template std::ostream &operator<<(std::ostream &stream, const AliCutValueRange<float> &val);
}

}

