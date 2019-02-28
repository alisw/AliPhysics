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

/// \file GPUTRDTrackData.h
/// \brief This is a flat data structure for transporting TRD tracks via network between the components

/// \author Sergey Gorbunov, Ole Schmidt

#ifndef GPUTRDTRACKDATA_H
#define GPUTRDTRACKDATA_H

struct GPUTRDTrackDataRecord {
  float mAlpha;              // azimuthal angle of reference frame
  float fX;                  // x: radial distance
  float fY;                  // local Y-coordinate of a track (cm)
  float fZ;                  // local Z-coordinate of a track (cm)
  float mSinPhi;             // local sine of the track momentum azimuthal angle
  float fTgl;                // tangent of the track momentum dip angle
  float fq1Pt;               // 1/pt (1/(GeV/c))
  float fC[15];              // covariance matrix
  int fTPCTrackID;           // id of corresponding TPC track
  int fAttachedTracklets[6]; // IDs for attached tracklets sorted by layer

  int GetNTracklets() const
  {
    int n = 0;
    for (int i = 0; i < 6; i++) {
      if (fAttachedTracklets[i] >= 0) {
        n++;
      }
    }
    return n;
  }
};

typedef struct GPUTRDTrackDataRecord GPUTRDTrackDataRecord;

struct GPUTRDTrackData {
  unsigned int fCount; // number of tracklets
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  GPUTRDTrackDataRecord fTracks[1]; // array of tracklets
#else
  GPUTRDTrackDataRecord fTracks[0]; // array of tracklets
#endif
  static size_t GetSize(unsigned int nTracks)
  {
    return sizeof(GPUTRDTrackData) + nTracks * sizeof(GPUTRDTrackDataRecord);
  }
  size_t GetSize() const { return GetSize(fCount); }
};

typedef struct GPUTRDTrackData GPUTRDTrackData;

#endif // GPUTRDTRACKDATA_H
