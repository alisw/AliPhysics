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
// class GPUTPCGMTrackParam;
// typedef GPUTPCGMTrackParam TRDBaseTrack;
#elif defined(TRD_TRACK_TYPE_O2)
class GPUTPCGMTrackParam;
typedef GPUTPCGMTrackParam TRDBaseTrack;
#endif

#ifdef GPUCA_ALIROOT_LIB
typedef AliTrackerBase TRDBasePropagator;
// class GPUTPCGMPropagator;
// typedef GPUTPCGMPropagator TRDBasePropagator;
#else
class GPUTPCGMPropagator;
typedef GPUTPCGMPropagator TRDBasePropagator;
#endif

template <class T>
class trackInterface;
template <class T>
class propagatorInterface;
template <class T>
class GPUTRDTrack_t;
typedef GPUTRDTrack_t<trackInterface<TRDBaseTrack>> GPUTRDTrack;
typedef propagatorInterface<TRDBasePropagator> GPUTRDPropagator;

#if defined(GPUCA_ALIGPUCODE) && !defined(GPUCA_ALIROOT_LIB) && !defined(__CLING__) && !defined(__ROOTCLING__) && !defined(G__ROOT)
#define Error(...)
#define Warning(...)
#define Info(...)
#endif
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTRDDEF_H
