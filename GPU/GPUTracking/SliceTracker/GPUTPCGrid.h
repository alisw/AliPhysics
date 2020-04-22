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

/// \file GPUTPCGrid.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCGRID_H
#define GPUTPCGRID_H

#include "GPUTPCDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCGrid
 *
 * 2-dimensional grid of pointers.
 * pointers to (y,z)-like objects are assigned to the corresponding grid bin
 * used by GPUTPCTracker to speed-up the hit operations
 * grid axis are named Z,Y to be similar to TPC row coordinates.
 */
MEM_CLASS_PRE()
class GPUTPCGrid
{
 public:
  GPUd() void CreateEmpty();
  GPUd() void Create(float yMin, float yMax, float zMin, float zMax, int ny, int nz);

  GPUd() int GetBin(float Y, float Z) const;
  /**
 * returns -1 if the row is empty == no hits
 */
  GPUd() int GetBinBounded(float Y, float Z) const;
  GPUd() void GetBin(float Y, float Z, int* const bY, int* const bZ) const;
  GPUd() void GetBinArea(float Y, float Z, float dy, float dz, int& bin, int& ny, int& nz) const;

  GPUd() unsigned int N() const { return mN; }
  GPUd() unsigned int Ny() const { return mNy; }
  GPUd() unsigned int Nz() const { return mNz; }
  GPUd() float YMin() const { return mYMin; }
  GPUd() float YMax() const { return mYMax; }
  GPUd() float ZMin() const { return mZMin; }
  GPUd() float ZMax() const { return mZMax; }
  GPUd() float StepYInv() const { return mStepYInv; }
  GPUd() float StepZInv() const { return mStepZInv; }

 private:
  friend class GPUTPCNeighboursFinder;

  unsigned int mNy; //* N bins in Y
  unsigned int mNz; //* N bins in Z
  unsigned int mN;  //* total N bins
  float mYMin;      //* minimal Y value
  float mYMax;      //* maximal Y value
  float mZMin;      //* minimal Z value
  float mZMax;      //* maximal Z value
  float mStepYInv;  //* inverse bin size in Y
  float mStepZInv;  //* inverse bin size in Z
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCGRID_H
