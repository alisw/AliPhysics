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

/// \file GPUTPCMCInfo.h
/// \author David Rohr

#ifndef GPUTPCMCINFO_H
#define GPUTPCMCINFO_H

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUTPCMCInfo {
  int charge;
  char prim;
  char primDaughters;
  int pid;
  float x;
  float y;
  float z;
  float pX;
  float pY;
  float pZ;
  float genRadius;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
