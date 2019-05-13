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

/// \file GPUTRDTrackPoint.h
/// \brief This is a flat data structure for transporting TRD track points via network between the components

/// \author Sergey Gorbunov, Ole Schmidt

#ifndef GPUTRDTRACKPOINT_H
#define GPUTRDTRACKPOINT_H

struct GPUTRDTrackPoint {
  float fX[3];
  short fVolumeId;
};

struct GPUTRDTrackPointData {
  unsigned int fCount; // number of space points
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  GPUTRDTrackPoint fPoints[1]; // array of space points
#else
  GPUTRDTrackPoint fPoints[0]; // array of space points
#endif
};

typedef struct GPUTRDTrackPointData GPUTRDTrackPointData;

#endif // GPUTRDTRACKPOINT_H
