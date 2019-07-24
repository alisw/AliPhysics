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

/// \file GPUDisplayQuaternion.cxx
/// \author David Rohr

#include "GPUDisplay.h"
#ifdef GPUCA_BUILD_EVENT_DISPLAY

#include <cmath>
using namespace GPUCA_NAMESPACE::gpu;

void GPUDisplay::createQuaternionFromMatrix(float* v, const float* mat)
{
  if (mat[0] > mat[5] && mat[0] > mat[10]) {
    const float S = sqrt(std::max(0.f, 1.0f + mat[0] - mat[5] - mat[10])) * 2;
    v[0] = 0.25 * S;
    v[1] = (mat[4] + mat[1]) / S;
    v[2] = (mat[2] + mat[8]) / S;
    v[3] = (mat[9] - mat[6]) / S;
  } else if (mat[5] > mat[10]) {
    const float S = sqrt(std::max(0.f, 1.0f + mat[5] - mat[0] - mat[10])) * 2;
    v[1] = 0.25 * S;
    v[0] = (mat[4] + mat[1]) / S;
    v[2] = (mat[9] + mat[6]) / S;
    v[3] = (mat[2] - mat[8]) / S;
  } else {
    float S = sqrt(std::max(0.f, 1.0f + mat[10] - mat[0] - mat[5])) * 2;
    v[2] = 0.25 * S;
    if (fabsf(S) < 0.001) {
      S = 1;
    }
    v[0] = (mat[2] + mat[8]) / S;
    v[1] = (mat[9] + mat[6]) / S;
    v[3] = (mat[4] - mat[1]) / S;
  }
  if (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3] < 0.0001) {
    v[3] = 1;
  }
}

#endif
