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

/// \file genEvents.h
/// \author Sergey Gorbunov

#ifndef GENEVENTS_H
#define GENEVENTS_H

#include "GPUCommonDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUChainTracking;
struct GPUParam;
class GPUTPCGMPhysicalTrackModel;
#if !defined(BUILD_QA) || defined(_WIN32)
class genEvents
{
 public:
  genEvents(GPUChainTracking* rec) {}
  void InitEventGenerator() {}
  int GenerateEvent(const GPUParam& sliceParam, char* filename) { return 1; }
  void FinishEventGenerator() {}

  static void RunEventGenerator(GPUChainTracking* rec){};
};

#else

class genEvents
{
 public:
  genEvents(GPUChainTracking* rec) : mRec(rec) {}
  void InitEventGenerator();
  int GenerateEvent(const GPUParam& sliceParam, char* filename);
  void FinishEventGenerator();

  static void RunEventGenerator(GPUChainTracking* rec);

 private:
  int GetSlice(double GlobalPhi);
  int GetDSlice(double LocalPhi);
  double GetSliceAngle(int iSlice);
  int RecalculateSlice(GPUTPCGMPhysicalTrackModel& t, int& iSlice);
  double GetGaus(double sigma);

  TH1F* mClusterError[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };

  struct GenCluster {
    int sector;
    int row;
    int mcID;
    float x;
    float y;
    float z;
    unsigned int id;
  };

  const double mTwoPi = 2 * M_PI;
  const double mSliceDAngle = mTwoPi / 18.;
  const double mSliceAngleOffset = mSliceDAngle / 2;

  GPUChainTracking* mRec;
};

#endif
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
