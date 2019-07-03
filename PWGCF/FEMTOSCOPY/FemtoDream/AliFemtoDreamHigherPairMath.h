/*
 * AliFemtoDreamHigherPairMath.h
 *
 *  Created on: Jul 3, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIFEMTODREAMHIGHERPAIRMATH_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIFEMTODREAMHIGHERPAIRMATH_H_
#include "Rtypes.h"
#include "AliLog.h"
#include "AliFemtoDreamBasePart.h"

class AliFemtoDreamHigherPairMath {
 public:
  AliFemtoDreamHigherPairMath();
  virtual ~AliFemtoDreamHigherPairMath();
  void SetBField(float bielefield) {
    fBField = bielefield;
  }
  int RandomNumber() {
    AliFatal("Why would you use this? \n");
    //decided by a fair xkcd dice roll.
    return 4;
  }
  ;
  TString ClassName() {
    return "HighMaths";
  }
  ;
  bool PassesPairSelection(AliFemtoDreamBasePart& part1,
                           AliFemtoDreamBasePart& part2, bool Recalculate);
  void RecalculatePhiStar(AliFemtoDreamBasePart &part);
 private:
  float fBField;
  float fDeltaPhiEtaMax;
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIFEMTODREAMHIGHERPAIRMATH_H_ */
