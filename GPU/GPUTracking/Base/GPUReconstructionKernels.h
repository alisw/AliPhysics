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
GPUCA_KRNL((GPUTPCNeighboursFinder                       ), (single), (), ())
GPUCA_KRNL((GPUTPCNeighboursCleaner                      ), (single), (), ())
GPUCA_KRNL((GPUTPCStartHitsFinder                        ), (single), (), ())
GPUCA_KRNL((GPUTPCStartHitsSorter                        ), (single), (), ())
GPUCA_KRNL((GPUTPCTrackletConstructor, singleSlice       ), (single), (), ())
GPUCA_KRNL((GPUTPCTrackletConstructor, allSlices         ), (single), (), ())
GPUCA_KRNL((GPUTPCTrackletSelector                       ), (both),   (), ())
GPUCA_KRNL((GPUMemClean16                                ), (),       (, GPUPtr1(void*, ptr), unsigned long size), (, GPUPtr2(void*, ptr), size))
#ifndef GPUCA_OPENCL1
GPUCA_KRNL((GPUTPCGMMergerTrackFit                       ), (),       (), ())
#ifdef HAVE_O2HEADERS
GPUCA_KRNL((GPUTRDTrackerGPU                             ), (),       (), ())
GPUCA_KRNL((GPUITSFitterKernel                           ), (),       (), ())
GPUCA_KRNL((GPUTPCConvertKernel                          ), (),       (), ())
GPUCA_KRNL((GPUTPCCompressionKernels,   step0attached    ), (),       (), ())
GPUCA_KRNL((GPUTPCCompressionKernels,   step1unattached  ), (),       (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, fillChargeMap    ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, resetMaps        ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, findPeaks        ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, noiseSuppression ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, updatePeaks      ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, countPeaks       ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, computeClusters  ), (single), (), ())
GPUCA_KRNL((GPUTPCClusterFinderKernels, nativeScanUpStart), (single), (, int iBuf, int stage), (, iBuf, stage))
GPUCA_KRNL((GPUTPCClusterFinderKernels, nativeScanUp     ), (single), (, int iBuf, int nElems), (, iBuf, nElems))
GPUCA_KRNL((GPUTPCClusterFinderKernels, nativeScanTop    ), (single), (, int iBuf, int nElems), (, iBuf, nElems))
GPUCA_KRNL((GPUTPCClusterFinderKernels, nativeScanDown   ), (single), (, int iBuf, unsigned int offset, int nElems), (, iBuf, offset, nElems))
GPUCA_KRNL((GPUTPCClusterFinderKernels, compactDigit     ), (single), (, int iBuf, int stage, GPUPtr1(deprecated::PackedDigit*, in), GPUPtr1(deprecated::PackedDigit*, out)), (, iBuf, stage, GPUPtr2(deprecated::PackedDigit*, in), GPUPtr2(deprecated::PackedDigit*, out)))
GPUCA_KRNL((GPUTPCClusterFinderKernels, decodeZS         ), (single), (), ())
#endif
#endif
// clang-format on
