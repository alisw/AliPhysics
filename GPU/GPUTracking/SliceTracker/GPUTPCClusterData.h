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

/// \file GPUTPCClusterData.h
/// \author Matthias Kretz, Sergey Gorbunov, David Rohr

#ifndef GPUTPCCLUSTERDATA_H
#define GPUTPCCLUSTERDATA_H

#include "GPUTPCDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUTPCClusterData {
  int id;
  short row;
  short flags;
  float x;
  float y;
  float z;
  float amp;
#ifdef GPUCA_FULL_CLUSTERDATA
  float pad;
  float time;
  float ampMax;
  float sigmaPad2;
  float sigmaTime2;
#endif
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // CLUSTERDATA_H
