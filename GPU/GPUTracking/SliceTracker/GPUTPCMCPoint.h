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

/// \file GPUTPCMCPoint.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCMCPOINT_H
#define GPUTPCMCPOINT_H

#include "GPUTPCDef.h"

/**
 * @class GPUTPCMCPoint
 * store MC point information for GPUTPCPerformance
 */
class GPUTPCMCPoint
{
 public:
  GPUTPCMCPoint();

  float X() const { return fX; }
  float Y() const { return fY; }
  float Z() const { return fZ; }
  float Sx() const { return fSx; }
  float Sy() const { return fSy; }
  float Sz() const { return fSz; }
  float Time() const { return fTime; }
  int ISlice() const { return mISlice; }
  int TrackID() const { return fTrackID; }

  void SetX(float v) { fX = v; }
  void SetY(float v) { fY = v; }
  void SetZ(float v) { fZ = v; }
  void SetSx(float v) { fSx = v; }
  void SetSy(float v) { fSy = v; }
  void SetSz(float v) { fSz = v; }
  void SetTime(float v) { fTime = v; }
  void SetISlice(int v) { mISlice = v; }
  void SetTrackID(int v) { fTrackID = v; }

  static bool Compare(const GPUTPCMCPoint& p1, const GPUTPCMCPoint& p2)
  {
    if (p1.fTrackID != p2.fTrackID) {
      return (p1.fTrackID < p2.fTrackID);
    }
    if (p1.mISlice != p2.mISlice) {
      return (p1.mISlice < p2.mISlice);
    }
    return (p1.Sx() < p2.Sx());
  }

  static bool CompareSlice(const GPUTPCMCPoint& p, int slice) { return (p.ISlice() < slice); }

  static bool CompareX(const GPUTPCMCPoint& p, float X) { return (p.Sx() < X); }

 protected:
  float fX;     //* global X position
  float fY;     //* global Y position
  float fZ;     //* global Z position
  float fSx;    //* slice X position
  float fSy;    //* slice Y position
  float fSz;    //* slice Z position
  float fTime;  //* time
  int mISlice;  //* slice number
  int fTrackID; //* mc track number
};

#endif // GPUTPCMCPOINT_H
