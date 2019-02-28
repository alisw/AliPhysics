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

/// \file GPUCommonFairLogger.h
/// \author David Rohr

#ifndef GPUCOMMONFAIRLOGGER_H
#define GPUCOMMONFAIRLOGGER_H

#if defined(GPUCA_STANDALONE) || defined(GPUCA_ALIROOT_LIB) || defined(GPUCA_GPULIBRARY)

#include <iostream>
#define LOG(type) std::cout
namespace FairLogger
{
static constexpr const char* endl = "\n";
}

#else

#include <FairLogger.h>

#endif

#endif
