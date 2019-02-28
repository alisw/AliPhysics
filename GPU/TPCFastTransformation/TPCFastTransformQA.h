//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file  TPCFastTransformManager.h
/// \brief Definition of TPCFastTransformManager class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_TPCFASTTRANSFORMQA_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_TPCFASTTRANSFORMQA_H

#include "GPUCommonDef.h"
#include "TPCFastTransformManager.h"
#include <cmath>
#include <iostream>

#include "Rtypes.h"
#include "TString.h"
#include "AliTPCTransform.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The TPCFastTransformQA class does performance check for TPCFastTransformation object
///

class TPCFastTransformQA
{
 public:
  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  TPCFastTransformQA();

  /// Copy constructor: disabled
  TPCFastTransformQA(const TPCFastTransformQA&) CON_DELETE;

  /// Assignment operator: disabled
  TPCFastTransformQA& operator=(const TPCFastTransformQA&) CON_DELETE;

  /// Destructor
  ~TPCFastTransformQA() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// create fast transformation and perform a quality check
  int doQA(Long_t TimeStamp);

  /// create perform quality check
  int doQA(const TPCFastTransform& fastTransform);

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);
  TString mError; ///< error string
};

inline int TPCFastTransformQA::storeError(int code, const char* msg)
{
  mError = msg;
  std::cout << msg << std::endl;
  return code;
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
