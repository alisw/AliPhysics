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

/// \file GPUTRDTrackletWord.h
/// \brief TRD Tracklet word for GPU tracker - 32bit tracklet info + half chamber ID + index

/// \author Ole Schmidt

#ifndef GPUTRDTRACKLETWORD_H
#define GPUTRDTRACKLETWORD_H

#include "GPUDef.h"

class AliTRDtrackletWord;
class AliTRDtrackletMCM;

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class GPUTRDTrackletWord
{
 public:
  GPUd() GPUTRDTrackletWord(unsigned int trackletWord = 0);
  GPUd() GPUTRDTrackletWord(unsigned int trackletWord, int hcid, int id);
  GPUdDefault() GPUTRDTrackletWord(const GPUTRDTrackletWord& rhs) CON_DEFAULT;
  GPUdDefault() GPUTRDTrackletWord& operator=(const GPUTRDTrackletWord& rhs) CON_DEFAULT;
  GPUdDefault() ~GPUTRDTrackletWord() CON_DEFAULT;
#ifndef GPUCA_GPUCODE_DEVICE
  GPUTRDTrackletWord(const AliTRDtrackletWord& rhs);
  GPUTRDTrackletWord(const AliTRDtrackletMCM& rhs);
  GPUTRDTrackletWord& operator=(const AliTRDtrackletMCM& rhs);
#endif

  // ----- Override operators < and > to enable tracklet sorting by HCId -----
  GPUd() bool operator<(const GPUTRDTrackletWord& t) const { return (GetHCId() < t.GetHCId()); }
  GPUd() bool operator>(const GPUTRDTrackletWord& t) const { return (GetHCId() > t.GetHCId()); }
  GPUd() bool operator<=(const GPUTRDTrackletWord& t) const { return (GetHCId() < t.GetHCId()) || (GetHCId() == t.GetHCId()); }

  // ----- Getters for contents of tracklet word -----
  GPUd() int GetYbin() const;
  GPUd() int GetdY() const;
  GPUd() int GetZbin() const { return ((mTrackletWord >> 20) & 0xf); }
  GPUd() int GetPID() const { return ((mTrackletWord >> 24) & 0xff); }

  GPUd() int GetId() const { return mId; }

  // ----- Getters for offline corresponding values -----
  GPUd() double GetPID(int /* is */) const { return (double)GetPID() / 256.f; }
  GPUd() int GetDetector() const { return mHCId / 2; }
  GPUd() int GetHCId() const { return mHCId; }
  GPUd() float GetdYdX() const { return (GetdY() * 140e-4f / 3.f); }
  GPUd() float GetY() const { return (GetYbin() * 160e-4f); }
  GPUd() unsigned int GetTrackletWord() const { return mTrackletWord; }

  GPUd() void SetTrackletWord(unsigned int trackletWord) { mTrackletWord = trackletWord; }
  GPUd() void SetDetector(int id) { mHCId = 2 * id + (GetYbin() < 0 ? 0 : 1); }
  GPUd() void SetId(int id) { mId = id; }
  GPUd() void SetHCId(int id) { mHCId = id; }

 protected:
  int mId;                    // index in tracklet array
  int mHCId;                  // half-chamber ID
  unsigned int mTrackletWord; // tracklet word: PID | Z | deflection length | Y
                              //          bits:   8   4            7          13
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTRDTRACKLETWORD_H
