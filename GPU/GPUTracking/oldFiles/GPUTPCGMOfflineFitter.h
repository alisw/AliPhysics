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

/// \file GPUTPCGMOfflineFitter.h
/// \author Sergey Gorbunov

#ifndef GPUTPCGMOfflineFitter_H
#define GPUTPCGMOfflineFitter_H

#if (defined(GPUCA_ALIROOT_LIB) && !defined(GPUCA_GPUCODE))

#include "GPUParam.h"
#include "AliTPCtracker.h"

class GPUTPCGMMergedTrack;
class GPUTPCGMMergedTrackHit;
class AliTPCclusterMI;
class GPUTPCGMPolynomialField;

class GPUTPCGMOfflineFitter : public AliTPCtracker
{
 public:
  GPUTPCGMOfflineFitter();
  ~GPUTPCGMOfflineFitter();

  void Initialize(const GPUParam& hltParam, Long_t TimeStamp, bool isMC);

  void RefitTrack(GPUTPCGMMergedTrack& track, const GPUTPCGMPolynomialField* field, GPUTPCGMMergedTrackHit* clusters);

  int CreateTPCclusterMI(const GPUTPCGMMergedTrackHit& h, AliTPCclusterMI& c);

  bool FitOffline(const GPUTPCGMPolynomialField* field, GPUTPCGMMergedTrack& gmtrack, GPUTPCGMMergedTrackHit* clusters, int& N);

 private:
  GPUParam fCAParam;
};

#endif

#endif
