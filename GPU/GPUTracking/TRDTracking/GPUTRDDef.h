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

/// \file GPUTRDDef.h
/// \author David Rohr

#ifndef GPUTRDDEF_H
#define GPUTRDDEF_H

#include "GPUCommonDef.h"

#ifdef GPUCA_ALIROOT_LIB
#define TRD_TRACK_TYPE_ALIROOT
#else
#define TRD_TRACK_TYPE_O2
#endif

#ifdef GPUCA_ALIROOT_LIB
class AliExternalTrackParam;
class AliTrackerBase;
#else
namespace o2
{
namespace dataformats
{
class TrackTPCITS;
} // namespace dataformats
namespace base
{
class Propagator;
} // namespace base
} // namespace o2
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{

#ifdef GPUCA_ALIROOT_LIB
typedef double My_Float;
#else
typedef float My_Float;
#endif

#if defined(TRD_TRACK_TYPE_ALIROOT)
typedef AliExternalTrackParam TRDBaseTrack;
typedef AliExternalTrackParam TRDBaseTrackGPU;
// class GPUTPCGMTrackParam;
// typedef GPUTPCGMTrackParam TRDBaseTrack;
#elif defined(TRD_TRACK_TYPE_O2)
typedef o2::dataformats::TrackTPCITS TRDBaseTrack;
class GPUTPCGMTrackParam;
typedef GPUTPCGMTrackParam TRDBaseTrackGPU;
#endif

#ifdef GPUCA_ALIROOT_LIB
typedef AliTrackerBase TRDBasePropagator;
typedef AliTrackerBase TRDBasePropagatorGPU;
// class GPUTPCGMPropagator;
// typedef GPUTPCGMPropagator TRDBasePropagator;
#else
typedef o2::base::Propagator TRDBasePropagator;
class GPUTPCGMPropagator;
typedef GPUTPCGMPropagator TRDBasePropagatorGPU;
#endif

template <class T>
class trackInterface;
template <class T>
class propagatorInterface;
template <class T>
class GPUTRDTrack_t;
// clang-format off
typedef GPUTRDTrack_t<trackInterface<TRDBaseTrack> > GPUTRDTrack; // Need pre-c++11 compliant formatting
typedef GPUTRDTrack_t<trackInterface<TRDBaseTrackGPU> > GPUTRDTrackGPU;
// clang-foramt on
typedef propagatorInterface<TRDBasePropagator> GPUTRDPropagator;
typedef propagatorInterface<TRDBasePropagatorGPU> GPUTRDPropagatorGPU;

template <class T, class P>
class GPUTRDTracker_t;
typedef GPUTRDTracker_t<GPUTRDTrack, GPUTRDPropagator> GPUTRDTracker;
typedef GPUTRDTracker_t<GPUTRDTrackGPU, GPUTRDPropagatorGPU> GPUTRDTrackerGPU;

#if defined(GPUCA_ALIGPUCODE) && !defined(GPUCA_ALIROOT_LIB) && !defined(__CLING__) && !defined(__ROOTCLING__) && !defined(G__ROOT)
#define Error(...)
#define Warning(...)
#define Info(...)
#endif
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTRDDEF_H
