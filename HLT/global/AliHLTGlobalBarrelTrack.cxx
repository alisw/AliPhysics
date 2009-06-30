// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTGlobalBarrelTrack.cxx
    @author Matthias Richter
    @date   2009-06-24
    @brief  An AliKalmanTrack implementation for global HLT barrel tracks.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTGlobalBarrelTrack.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalBarrelTrack)

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack()
: AliKalmanTrack()
  , fPoints()
  , fLastX(0.0)
  , fLastY(0.0)
  , fTrackID(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack(const AliHLTGlobalBarrelTrack& t)
  : AliKalmanTrack(t)
  , fPoints()
  , fLastX(t.GetLastPointX())
  , fLastY(t.GetLastPointY())
  , fTrackID(t.TrackID())
{
  // see header file for class documentation
  fPoints.assign(t.fPoints.begin(), t.fPoints.end());
}

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack(const AliHLTExternalTrackParam& p)
  : AliKalmanTrack()
  , fPoints()
  , fLastX(p.fLastX)
  , fLastY(p.fLastY)
  , fTrackID(p.fTrackID)
{
  // see header file for class documentation

  // the 5 track parameters are named in the AliHLTExternalTrackParam
  // while AliExternalTrackParam just uses an array[5]
  // the members have the some order, fY is the first one
  Set(p.fX, p.fAlpha, &p.fY, p.fC);
  SetPoints(p.fPointIDs, p.fNPoints);
}

template <class c>
AliHLTGlobalBarrelTrack& AliHLTGlobalBarrelTrack::operator=(const c& p)
{
  // see header file for class documentation
  this->~AliHLTGlobalBarrelTrack();
  new (this) AliHLTGlobalBarrelTrack(p);
  return *this;
}

AliHLTGlobalBarrelTrack::~AliHLTGlobalBarrelTrack()
{
  // see header file for class documentation
}

int AliHLTGlobalBarrelTrack::ConvertTrackDataArray(const AliHLTTracksData* pTracks, unsigned sizeInByte, vector<AliHLTGlobalBarrelTrack> &tgtArray)
{
  // see header file for class documentation
  int iResult=0;
  tgtArray.clear();
  if (!pTracks || sizeInByte<sizeof(AliHLTTracksData) || pTracks->fCount==0) return 0;

  const AliHLTUInt8_t* pEnd=reinterpret_cast<const AliHLTUInt8_t*>(pTracks);
  pEnd+=sizeInByte;

  tgtArray.resize(pTracks->fCount);
  const AliHLTUInt8_t* pCurrent=reinterpret_cast<const AliHLTUInt8_t*>(pTracks->fTracklets);
  for (unsigned i=0; i<pTracks->fCount; i++) {
    if (pCurrent+sizeof(AliHLTExternalTrackParam)>pEnd) {
      iResult=-EINVAL; break;
    }
    const AliHLTExternalTrackParam* track=reinterpret_cast<const AliHLTExternalTrackParam*>(pCurrent);
    if (pCurrent+sizeof(AliHLTExternalTrackParam)+track->fNPoints*sizeof(UInt_t)>pEnd) {
      iResult=-EINVAL; break;
    }
    tgtArray[i]=*track;
    pCurrent+=sizeof(AliHLTExternalTrackParam)+track->fNPoints*sizeof(UInt_t);
  }
  if (iResult<0) tgtArray.clear();
  else iResult=tgtArray.size();
  return iResult;
}

UInt_t AliHLTGlobalBarrelTrack::GetNumberOfPoints() const
{
  // see header file for class documentation
  return fPoints.size();
}

const UInt_t* AliHLTGlobalBarrelTrack::GetPoints() const
{
  // see header file for class documentation
  if (fPoints.size()==0) return NULL;
  return &fPoints[0];
}

int AliHLTGlobalBarrelTrack::SetPoints(const UInt_t* pArray, UInt_t arraySize)
{
  // see header file for class documentation
  if (!pArray || arraySize==0) return 0;
  fPoints.resize(arraySize);
  for (unsigned i=0; i<arraySize; i++) fPoints[i]=pArray[i];
  return fPoints.size();
}

void AliHLTGlobalBarrelTrack::Print(Option_t* option) const
{
  // see header file for class documentation
  cout << "********* Track Id: " << fTrackID << " *******************" << endl;
  AliExternalTrackParam::Print(option);
//   cout << "  Alpha "     << GetAlpha();
//   cout << "  X "         << GetX();
//   cout << "  Y "         << GetY();
//   cout << "  Z "         << GetZ() << endl;
//   cout << "  Snp "       << GetSnp();
//   cout << "  Tgl "       << GetTgl();
//   cout << "  Signed1Pt " << GetSigned1Pt() << endl;
}
