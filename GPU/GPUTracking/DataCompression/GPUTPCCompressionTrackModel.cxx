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

/// \file GPUTPCCompressionTrackModel.cxx
/// \author David Rohr

#include "GPUTPCCompressionTrackModel.h"
#include "GPUConstantMem.h"
#include "GPUParam.inc"

using namespace GPUCA_NAMESPACE::gpu;

// ATTENTION! This track model is used for the data compression.
// Changes to the propagation and fit will prevent the decompression of data
// encoded with the old version!!!

#ifdef GPUCA_COMPRESSION_TRACK_MODEL_MERGER
GPUd() void GPUTPCCompressionTrackModel::Init(float x, float y, float z, float alpha, unsigned char qPt, const GPUParam& param)
{
  static constexpr float kRho = 1.025e-3f;  // 0.9e-3;
  static constexpr float kRadLen = 29.532f; // 28.94;
  mProp.SetMaterial(kRadLen, kRho);
  mProp.SetMaxSinPhi(GPUCA_MAX_SIN_PHI);
  mProp.SetToyMCEventsFlag(false);
  mProp.SetSeedingErrors(true); // Larger errors for seeds, better since we don't start with good hypothesis
  mProp.SetFitInProjections(true);
  mProp.SetPolynomialField(&param.polynomialField);
  mTrk.X() = x;
  mTrk.Y() = y;
  mTrk.Z() = z;
  mTrk.SinPhi() = 0;
  mTrk.DzDs() = 0;
  mTrk.QPt() = (qPt - 127.f) * (20.f / 127.f);
  mTrk.ResetCovariance();
  mProp.SetTrack(&mTrk, alpha);
  mParam = &param;
  // GPUInfo("Initialized: x %f y %f z %f alpha %f qPt %f", x, y, z, alpha, mTrk.QPt());
}

GPUd() int GPUTPCCompressionTrackModel::Propagate(float x, float alpha)
{
  int retVal = mProp.PropagateToXAlpha(x, alpha, true);
  // GPUInfo("Propagated to: x %f y %f z %f alpha %f qPt %f", x, mTrk.Y(), mTrk.Z(), alpha, mTrk.QPt());
  return retVal;
}

GPUd() int GPUTPCCompressionTrackModel::Filter(float y, float z, int iRow)
{
  mTrk.ConstrainSinPhi();
  int retVal = mProp.Update(y, z, iRow, *mParam, 0, false, false);
  // GPUInfo("Filtered with %f %f: y %f z %f qPt %f", y, z, mTrk.Y(), mTrk.Z(), mTrk.QPt());
  return retVal;
}

GPUd() int GPUTPCCompressionTrackModel::Mirror()
{
  mProp.Mirror(true);
  // GPUInfo("Mirrored: y %f z %f qPt %f", mTrk.Y(), mTrk.Z(), mTrk.QPt());
  return 0;
}

#elif defined(GPUCA_COMPRESSION_TRACK_MODEL_SLICETRACKER)

#else // Default internal track model for compression

#endif
