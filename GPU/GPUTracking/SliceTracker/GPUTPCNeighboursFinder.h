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

/// \file GPUTPCNeighboursFinder.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCNEIGHBOURSFINDER_H
#define GPUTPCNEIGHBOURSFINDER_H

#include "GPUTPCDef.h"
#include "GPUTPCRow.h"
#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
MEM_CLASS_PRE()
class GPUTPCTracker;

/**
 * @class GPUTPCNeighboursFinder
 *
 */
class GPUTPCNeighboursFinder : public GPUKernelTemplate
{
 public:
  MEM_CLASS_PRE()
  struct GPUSharedMemory {
    int mNHits;   // n hits
    float mUpDx;  // x distance to the next row
    float mDnDx;  // x distance to the previous row
    float mUpTx;  // normalized x distance to the next row
    float mDnTx;  // normalized x distance to the previous row
    int mIRow;    // row number
    int mIRowUp;  // next row number
    int mIRowDn;  // previous row number
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0
    float mA1[GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP][GPUCA_THREAD_COUNT_FINDER];
    float mA2[GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP][GPUCA_THREAD_COUNT_FINDER];
    calink mB[GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP][GPUCA_THREAD_COUNT_FINDER];
#endif
    MEM_LG(GPUTPCRow)
    mRow, mRowUp, mRowDown;
  };

  typedef GPUconstantref() MEM_GLOBAL(GPUTPCTracker) processorType;
  GPUhdi() CONSTEXPRRET static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::TPCSliceTracking; }
  MEM_TEMPLATE()
  GPUhdi() static processorType* Processor(MEM_TYPE(GPUConstantMem) & processors)
  {
    return processors.tpcTrackers;
  }
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & smem, processorType& tracker);
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCNEIGHBOURSFINDER_H
