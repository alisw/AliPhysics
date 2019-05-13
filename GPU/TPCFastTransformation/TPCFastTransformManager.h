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

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_TPCFASTTRANSFORMMANAGER_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_TPCFASTTRANSFORMMANAGER_H

#include <cmath>

#include "GPUCommonDef.h"
#include "Rtypes.h"
#include "TString.h"
#include "AliTPCTransform.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class TPCFastTransform;

///
/// The TPCFastTransformManager class is to initialize TPCFastTransformation object
///

class TPCFastTransformManager
{
 public:
  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  TPCFastTransformManager();

  /// Copy constructor: disabled
  TPCFastTransformManager(const TPCFastTransformManager&) CON_DELETE;

  /// Assignment operator: disabled
  TPCFastTransformManager& operator=(const TPCFastTransformManager&) CON_DELETE;

  /// Destructor
  ~TPCFastTransformManager() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// Initializes TPCFastTransform object
  int create(TPCFastTransform& spline, AliTPCTransform* transform, Long_t TimeStamp);

  /// Updates the transformation with the new time stamp
  Int_t updateCalibration(TPCFastTransform& spline, Long_t TimeStamp);

  /// _______________  Utilities   ________________________

  AliTPCTransform* getOriginalTransform() { return mOrigTransform; }

  ///  Gives error string
  const char* getLastError() const { return mError.Data(); }

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);

  TString mError;                  ///< error string
  AliTPCTransform* mOrigTransform; ///< transient
  int fLastTimeBin;                ///< last calibrated time bin
};

inline int TPCFastTransformManager::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
