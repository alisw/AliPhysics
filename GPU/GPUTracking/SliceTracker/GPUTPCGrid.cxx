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

/// \file GPUTPCGrid.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCGrid.h"
#include "GPUCommonMath.h"
using namespace GPUCA_NAMESPACE::gpu;

#ifndef assert
#include <cassert>
#endif

MEM_CLASS_PRE()
void MEM_LG(GPUTPCGrid)::CreateEmpty()
{
  // Create an empty grid
  mYMin = 0.f;
  mYMax = 1.f;
  mZMin = 0.f;
  mZMax = 1.f;

  mNy = 1;
  mNz = 1;
  mN = 1;

  mStepYInv = 1.f;
  mStepZInv = 1.f;
}

MEM_CLASS_PRE()
GPUd() void MEM_LG(GPUTPCGrid)::Create(float yMin, float yMax, float zMin, float zMax, float sy, float sz)
{
  //* Create the grid
  mYMin = yMin;
  mZMin = zMin;

  mStepYInv = 1.f / sy;
  mStepZInv = 1.f / sz;

  mNy = static_cast<unsigned int>((yMax - mYMin) * mStepYInv + 1.f);
  mNz = static_cast<unsigned int>((zMax - mZMin) * mStepZInv + 1.f);

  mN = mNy * mNz;

  mYMax = mYMin + mNy * sy;
  mZMax = mZMin + mNz * sz;
}

MEM_CLASS_PRE()
GPUd() int MEM_LG(GPUTPCGrid)::GetBin(float Y, float Z) const
{
  //* get the bin pointer
  const int yBin = static_cast<int>((Y - mYMin) * mStepYInv);
  const int zBin = static_cast<int>((Z - mZMin) * mStepZInv);
  const int bin = zBin * mNy + yBin;
#ifndef GPUCA_GPUCODE
  assert(bin >= 0);
  assert(bin < static_cast<int>(mN));
#endif
  return bin;
}

MEM_CLASS_PRE()
GPUd() int MEM_LG(GPUTPCGrid)::GetBinBounded(float Y, float Z) const
{
  //* get the bin pointer
  const int yBin = static_cast<int>((Y - mYMin) * mStepYInv);
  const int zBin = static_cast<int>((Z - mZMin) * mStepZInv);
  int bin = zBin * mNy + yBin;
  if (bin >= static_cast<int>(mN)) {
    bin = mN - 1;
  }
  if (bin < 0) {
    bin = 0;
  }
  return bin;
}

MEM_CLASS_PRE()
GPUd() void MEM_LG(GPUTPCGrid)::GetBin(float Y, float Z, int* const bY, int* const bZ) const
{
  //* get the bin pointer

  int bbY = (int)((Y - mYMin) * mStepYInv);
  int bbZ = (int)((Z - mZMin) * mStepZInv);

  if (bbY >= (int)mNy) {
    bbY = mNy - 1;
  }
  if (bbY < 0) {
    bbY = 0;
  }
  if (bbZ >= (int)mNz) {
    bbZ = mNz - 1;
  }
  if (bbZ < 0) {
    bbZ = 0;
  }

  *bY = (unsigned int)bbY;
  *bZ = (unsigned int)bbZ;
}

MEM_CLASS_PRE()
GPUd() void MEM_LG(GPUTPCGrid)::GetBinArea(float Y, float Z, float dy, float dz, int& bin, int& ny, int& nz) const
{
  Y -= mYMin;
  int by = (int)((Y - dy) * mStepYInv);
  ny = (int)((Y + dy) * mStepYInv) - by;
  Z -= mZMin;
  int bz = (int)((Z - dz) * mStepZInv);
  nz = (int)((Z + dz) * mStepZInv) - bz;
  if (by >= (int)mNy) {
    by = mNy - 1;
  }
  if (by < 0) {
    by = 0;
  }
  if (bz >= (int)mNz) {
    bz = mNz - 1;
  }
  if (bz < 0) {
    bz = 0;
  }
  if (by + ny >= (int)mNy) {
    ny = mNy - 1 - by;
  }
  if (bz + nz >= (int)mNz) {
    nz = mNz - 1 - bz;
  }
  bin = bz * mNy + by;
}
