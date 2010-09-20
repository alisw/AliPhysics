// -*- mode: c++ -*-
#ifndef ALICALOCONSTANTS_H
#define ALICALOCONSTANTS_H

/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille  <perthomas.hille@yale.edu>           *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to   perthomas.hille@yale.edu                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

namespace CaloConstants
{
  namespace FitAlgorithm
  {
    enum fitAlgorithm { kLMS = 0, kCrude = 1, kPeakFinder = 2, kNeuralNet = 3, kFastFit= 4,
			kLogFit = 5, kStandard = 6,  kNONE = 7};
  };

  namespace ReturnCodes
  {
    enum kReturnCode {kFitPar=1, kDummy=-1, kCrude=-9, kNoFit=-99, kInvalid=-9999};// possible return values
  };

  namespace PeakFinderConstants
  {
    const int  MAXSTART = 3;
    const int  SAMPLERANGE = 15;
  };
};


//For easier notation
namespace Algo = CaloConstants::FitAlgorithm;
namespace Ret  = CaloConstants::ReturnCodes; 
namespace PF   = CaloConstants::PeakFinderConstants;

#endif
