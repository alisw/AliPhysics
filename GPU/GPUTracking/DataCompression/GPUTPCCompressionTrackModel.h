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

/// \file GPUTPCCompressionTrackModel.h
/// \author David Rohr

#ifndef GPUTPCCOMPRESSIONTRACKMODEL_H
#define GPUTPCCOMPRESSIONTRACKMODEL_H

// For debugging purposes, we provide means to use other track models
#define GPUCA_COMPRESSION_TRACK_MODEL_MERGER
// #define GPUCA_COMPRESSION_TRACK_MODEL_SLICETRACKER

#include "GPUDef.h"

#ifdef GPUCA_COMPRESSION_TRACK_MODEL_MERGER
#include "GPUTPCGMPropagator.h"
#include "GPUTPCGMTrackParam.h"

#elif defined(GPUCA_COMPRESSION_TRACK_MODEL_SLICETRACKER)
#include "GPUTPCTrackParam.h"

#else // Default internal track model for compression
#error Not yet implemented
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
// ATTENTION! This track model is used for the data compression.
// Changes to the propagation and fit will prevent the decompression of data
// encoded with the old version!!!

struct GPUParam;

class GPUTPCCompressionTrackModel
{
 public:
  GPUd() void Init(float x, float y, float z, float alpha, unsigned char qPt, const GPUParam& proc);
  GPUd() int Propagate(float x, float alpha);
  GPUd() int Filter(float y, float z, int iRow);
  GPUd() int Mirror();

#if defined(GPUCA_COMPRESSION_TRACK_MODEL_MERGER) || defined(GPUCA_COMPRESSION_TRACK_MODEL_SLICETRACKER)
  GPUd() float X() const
  {
    return mTrk.GetX();
  }
  GPUd() float Y() const { return mTrk.GetY(); }
  GPUd() float Z() const { return mTrk.GetZ(); }
  GPUd() float SinPhi() const { return mTrk.GetSinPhi(); }
  GPUd() float DzDs() const { return mTrk.GetDzDs(); }
  GPUd() float QPt() const { return mTrk.GetQPt(); }

#else // Default internal track model for compression

#endif

 protected:
  const GPUParam* mParam;

#ifdef GPUCA_COMPRESSION_TRACK_MODEL_MERGER
  GPUTPCGMPropagator mProp;
  GPUTPCGMTrackParam mTrk;

#elif defined(GPUCA_COMPRESSION_TRACK_MODEL_SLICETRACKER)
  GPUTPCTrackParam mTrk;
  float mAlpha;

#else // Default internal track model for compression

#endif
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
