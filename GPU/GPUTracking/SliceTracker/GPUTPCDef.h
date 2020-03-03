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

/// \file GPUTPCDef.h
/// \author David Rohr, Sergey Gorbunov

// clang-format off
#ifndef GPUTPCDEF_H
#define GPUTPCDEF_H

#include "GPUDef.h"

#define CALINK_INVAL ((calink) -1)
namespace GPUCA_NAMESPACE
{
namespace gpu
{
#if defined(GPUCA_O2_LIB) || defined(GPUCA_O2_INTERFACE)
typedef unsigned int calink;
typedef unsigned int cahit;
#else
typedef unsigned int calink;
typedef unsigned int cahit;
#endif
struct cahit2 { cahit x, y; };
}
} // GPUCA_NAMESPACE::GPU

#ifdef GPUCA_TPC_USE_STAT_ERROR
  #define GPUCA_TPC_RAW_PROPAGATE_PAD_ROW_TIME
#endif

#ifdef GPUCA_TPC_RAW_PROPAGATE_PAD_ROW_TIME // Needs full clusterdata
  #define GPUCA_FULL_CLUSTERDATA
#endif

#if defined(GPUCA_STANDALONE) || defined(GPUCA_GPUCODE) // No support for Full Field Propagator or Statistical errors
  #ifdef GPUCA_GM_USE_FULL_FIELD
    #undef GPUCA_GM_USE_FULL_FIELD
  #endif
  #ifdef GPUCA_TPC_USE_STAT_ERROR
    #undef GPUCA_TPC_USE_STAT_ERROR
  #endif
#endif

#ifdef GPUCA_EXTERN_ROW_HITS
  #define CA_GET_ROW_HIT(iRow) tracker.TrackletRowHits()[iRow * s.mNTracklets + r.mItr]
  #define CA_SET_ROW_HIT(iRow, val) tracker.TrackletRowHits()[iRow * s.mNTracklets + r.mItr] = val
#else
  #define CA_GET_ROW_HIT(iRow) tracklet.RowHit(iRow)
  #define CA_SET_ROW_HIT(iRow, val) tracklet.SetRowHit(iRow, val)
#endif

#endif //GPUDTPCEF_H
// clang format on
