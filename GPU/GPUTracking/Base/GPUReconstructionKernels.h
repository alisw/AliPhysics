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

/// \file GPUReconstructionKernels.h
/// \author David Rohr

// No header protection, this may be used multiple times
#include "GPUReconstructionKernelMacros.h"

// clang-format off
GPUCA_KRNL((GPUTPCNeighboursFinder                       ), (single, REG, (GPUCA_THREAD_COUNT_FINDER, GPUCA_BLOCK_COUNT_FINDER_MULTIPLIER)), (), ())
GPUCA_KRNL((GPUTPCNeighboursCleaner                      ), (single, REG, (GPUCA_THREAD_COUNT_CLEANER, 1)), (), ())
GPUCA_KRNL((GPUTPCStartHitsFinder                        ), (single, REG, (GPUCA_THREAD_COUNT, 1)), (), ())
GPUCA_KRNL((GPUTPCStartHitsSorter                        ), (single, REG, (GPUCA_THREAD_COUNT, 1)), (), ())
GPUCA_KRNL((GPUTPCTrackletConstructor, singleSlice       ), (single, REG, (GPUCA_THREAD_COUNT_CONSTRUCTOR, GPUCA_BLOCK_COUNT_CONSTRUCTOR_MULTIPLIER)), (), ())
GPUCA_KRNL((GPUTPCTrackletConstructor, allSlices         ), (single, REG, (GPUCA_THREAD_COUNT_CONSTRUCTOR, GPUCA_BLOCK_COUNT_CONSTRUCTOR_MULTIPLIER)), (), ())
GPUCA_KRNL((GPUTPCTrackletSelector                       ), (both, REG, (GPUCA_THREAD_COUNT_SELECTOR, GPUCA_BLOCK_COUNT_SELECTOR_MULTIPLIER)), (), ())
GPUCA_KRNL((GPUMemClean16                                ), (simple, REG, (GPUCA_THREAD_COUNT, 1)), (, GPUPtr1(void*, ptr), unsigned long size), (, GPUPtr2(void*, ptr), size))
#if !defined(GPUCA_OPENCL1) && (!defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE))
GPUCA_KRNL((GPUTPCGMMergerTrackFit                       ), (simple, REG, (GPUCA_THREAD_COUNT_FIT, 1)), (, int mode), (, mode))
#ifdef HAVE_O2HEADERS
GPUCA_KRNL((GPUTRDTrackerGPU                             ), (simple, REG, (GPUCA_THREAD_COUNT_TRD, 1)), (), ())
GPUCA_KRNL((GPUITSFitterKernel                           ), (simple, REG, (GPUCA_THREAD_COUNT_ITS, 1)), (), ())
GPUCA_KRNL((GPUTPCConvertKernel                          ), (simple, REG, (GPUCA_THREAD_COUNT_CONVERTER, 1)), (), ())
GPUCA_KRNL((GPUTPCCompressionKernels,   step0attached    ), (simple, REG, (GPUCA_THREAD_COUNT_COMPRESSION1, 1)), (), ())
GPUCA_KRNL((GPUTPCCompressionKernels,   step1unattached  ), (simple, REG, (GPUCA_THREAD_COUNT_COMPRESSION2, 1)), (), ())
GPUCA_KRNL((GPUTPCCFChargeMapFiller,    fillChargeMap    ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFChargeMapFiller,    resetMaps        ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFPeakFinder                           ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFNoiseSuppression,   noiseSuppression ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFNoiseSuppression,   updatePeaks      ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFDeconvolution                        ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFClusterizer                          ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFMCLabelFlattener,   setRowOffsets    ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (), ())
GPUCA_KRNL((GPUTPCCFMCLabelFlattener,   flatten          ), (single, REG, (GPUCA_THREAD_COUNT_CLUSTERER, 1)), (, unsigned int row, GPUPtr1(GPUTPCLinearLabels*, out)), (, row, GPUPtr2(GPUTPCLinearLabels*, out)))
GPUCA_KRNL((GPUTPCCFStreamCompaction,   nativeScanUpStart), (single, REG, (GPUCA_THREAD_COUNT_SCAN, 1)), (, int iBuf, int stage), (, iBuf, stage))
GPUCA_KRNL((GPUTPCCFStreamCompaction,   nativeScanUp     ), (single, REG, (GPUCA_THREAD_COUNT_SCAN, 1)), (, int iBuf, int nElems), (, iBuf, nElems))
GPUCA_KRNL((GPUTPCCFStreamCompaction,   nativeScanTop    ), (single, REG, (GPUCA_THREAD_COUNT_SCAN, 1)), (, int iBuf, int nElems), (, iBuf, nElems))
GPUCA_KRNL((GPUTPCCFStreamCompaction,   nativeScanDown   ), (single, REG, (GPUCA_THREAD_COUNT_SCAN, 1)), (, int iBuf, unsigned int offset, int nElems), (, iBuf, offset, nElems))
GPUCA_KRNL((GPUTPCCFStreamCompaction,   compactDigit     ), (single, REG, (GPUCA_THREAD_COUNT_SCAN, 1)), (, int iBuf, int stage, GPUPtr1(deprecated::PackedDigit*, in), GPUPtr1(deprecated::PackedDigit*, out)), (, iBuf, stage, GPUPtr2(deprecated::PackedDigit*, in), GPUPtr2(deprecated::PackedDigit*, out)))
GPUCA_KRNL((GPUTPCCFDecodeZS                             ), (single, REG, (GPUCA_THREAD_COUNT_CFDECODE, GPUCA_BLOCK_COUNT_DECODE_MULTIPLIER)), (), ())
#endif
#endif
// clang-format on
