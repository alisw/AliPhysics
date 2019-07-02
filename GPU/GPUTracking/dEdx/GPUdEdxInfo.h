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

/// \file GPUdEdxInfo.h
/// \author David Rohr

#ifndef GPUDEDXINFO_H
#define GPUDEDXINFO_H

#ifdef HAVE_O2HEADERS
#include "DataFormatsTPC/dEdxInfo.h"
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
#ifdef HAVE_O2HEADERS
using GPUdEdxInfo = o2::tpc::dEdxInfo;
#else
struct GPUdEdxInfo {
};
#endif
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
