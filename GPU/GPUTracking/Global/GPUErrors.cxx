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

/// \file GPUErrors.cxx
/// \author David Rohr

#include "GPUErrors.h"
#include "GPUDataTypes.h"
#include "GPUCommonMath.h"
#include "GPUDefMacros.h"
#include "GPULogging.h"

using namespace GPUCA_NAMESPACE::gpu;

#define GPUCA_MAX_ERRORS 255u

GPUd() void GPUErrors::raiseError(unsigned int code, unsigned int param1, unsigned int param2, unsigned int param3) const
{
  unsigned int pos = CAMath::AtomicAdd(mErrors, 1u);
  if (pos < GPUCA_MAX_ERRORS) {
    mErrors[4 * pos + 1] = code;
    mErrors[4 * pos + 2] = param1;
    mErrors[4 * pos + 3] = param2;
    mErrors[4 * pos + 4] = param3;
  }
}

#ifndef GPUCA_GPUCODE

#include <cstring>
#include <unordered_map>

unsigned int GPUErrors::getMaxErrors()
{
  return GPUCA_MAX_ERRORS;
}

void GPUErrors::clear()
{
  memset(mErrors, 0, GPUCA_MAX_ERRORS * sizeof(*mErrors));
}

static std::unordered_map<unsigned int, const char*> errorNames = {
#define GPUCA_ERROR_CODE(num, name) {num, GPUCA_M_STR(name)},
#include "GPUErrorCodes.h"
#undef GPUCA_ERROR_CODE
};

void GPUErrors::printErrors()
{
  for (unsigned int i = 0; i < std::min(*mErrors, GPUCA_MAX_ERRORS); i++) {
    const auto& it = errorNames.find(mErrors[4 * i + 1]);
    const char* errorName = it == errorNames.end() ? "INVALID ERROR CODE" : it->second;
    GPUError("GPU Error Code (%u:%u) %s : %u / %u / %u", i, mErrors[4 * i + 1], errorName, mErrors[4 * i + 2], mErrors[4 * i + 3], mErrors[4 * i + 4]);
  }
  if (*mErrors > GPUCA_MAX_ERRORS) {
    GPUError("Additional errors occured (codes not stored)");
  }
}

#endif
