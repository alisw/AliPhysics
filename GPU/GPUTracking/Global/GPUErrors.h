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

/// \file GPUErrors.h
/// \author David Rohr

#ifndef GPUERRORS_H
#define GPUERRORS_H

#include "GPUDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class GPUErrors
{
 public:
  enum errorNumbers {
#define GPUCA_ERROR_CODE(num, name) name = num,
#include "GPUErrorCodes.h"
#undef GPUCA_ERROR_CODE
  };

  GPUd() void raiseError(unsigned int code, unsigned int param1 = 0, unsigned int param2 = 0, unsigned int param3 = 0) const;
  GPUd() bool hasError() { return *mErrors > 0; }
  void setMemory(GPUglobalref() unsigned int* m) { mErrors = m; }
  void clear();
  void printErrors();
  static unsigned int getMaxErrors();

 private:
  GPUglobalref() unsigned int* mErrors;
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
