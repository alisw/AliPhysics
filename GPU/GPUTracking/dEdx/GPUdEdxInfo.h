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

#include "GPUDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUdEdxInfo {
  float dEdxTotIROC;
  float dEdxTotOROC1;
  float dEdxTotOROC2;
  float dEdxTotOROC3;
  float dEdxTotTPC;
  float dEdxMaxIROC;
  float dEdxMaxOROC1;
  float dEdxMaxOROC2;
  float dEdxMaxOROC3;
  float dEdxMaxTPC;
  unsigned char NHitsIROC;
  unsigned char NHitsSubThresholdIROC;
  unsigned char NHitsOROC1;
  unsigned char NHitsSubThresholdOROC1;
  unsigned char NHitsOROC2;
  unsigned char NHitsSubThresholdOROC2;
  unsigned char NHitsOROC3;
  unsigned char NHitsSubThresholdOROC3;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
