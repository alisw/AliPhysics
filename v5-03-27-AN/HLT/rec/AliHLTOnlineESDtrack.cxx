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

/// @file   AliHLTOnlineESDtrack.cxx
/// @author Matthias Richter
/// @date   2010-10-29
/// @brief  A streamlined container class for AliESDtrack.
/// @note   

#include "AliHLTOnlineESDtrack.h"
#include "AliESDtrack.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOnlineESDtrack)

AliHLTOnlineESDtrack::AliHLTOnlineESDtrack()
  : AliExternalTrackParam()
  , fCp(NULL)
  , fIp(NULL)
  , fTPCInner(NULL)
  , fOp(NULL)
  , fFlags(0)
  , fID(0)
  , fLabel(0)
  , fITSLabel(0)
  , fTPCLabel(0)
  , fTRDLabel(0)
  , fTrackLength(0)
  , fdTPC(0),fzTPC(0)
  , fCddTPC(0),fCdzTPC(0),fCzzTPC(0)
  , fCchi2TPC(0)
  , fD(0),fZ(0)
  , fCdd(0),fCdz(0),fCzz(0)
  , fCchi2(0)
  , fITSchi2(0)
  , fTPCchi2(0)
  , fTPCchi2Iter1(0)
  , fTPCncls(0)
  , fTPCnclsF(0)
  , fTPCnclsIter1(0)
  , fTPCnclsFIter1(0)
  , fITSncls(0)
  , fITSClusterMap(0)
  , fTRDncls(0)
  , fTRDncls0(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTOnlineESDtrack::AliHLTOnlineESDtrack(const AliHLTOnlineESDtrack& t)
  : AliExternalTrackParam(t)
  , fCp(t.fCp?(new AliExternalTrackParam(*t.fCp)):NULL)
  , fIp(t.fIp?(new AliExternalTrackParam(*t.fIp)):NULL)
  , fTPCInner(t.fTPCInner?(new AliExternalTrackParam(*t.fTPCInner)):NULL)
  , fOp(t.fOp?(new AliExternalTrackParam(*t.fOp)):NULL)
  , fFlags(t.fFlags)
  , fID(t.fID)
  , fLabel(t.fLabel)
  , fITSLabel(t.fITSLabel)
  , fTPCLabel(t.fTPCLabel)
  , fTRDLabel(t.fTRDLabel)
  , fTrackLength(t.fTrackLength)
  , fdTPC(t.fdTPC),fzTPC(t.fzTPC)
  , fCddTPC(t.fCddTPC),fCdzTPC(t.fCdzTPC),fCzzTPC(t.fCzzTPC)
  , fCchi2TPC(t.fCchi2TPC)
  , fD(t.fD),fZ(t.fZ)
  , fCdd(t.fCdd),fCdz(t.fCdz),fCzz(t.fCzz)
  , fCchi2(t.fCchi2)
  , fITSchi2(t.fITSchi2)
  , fTPCchi2(t.fTPCchi2)
  , fTPCchi2Iter1(t.fTPCchi2Iter1)
  , fTPCncls(t.fTPCncls)
  , fTPCnclsF(t.fTPCnclsF)
  , fTPCnclsIter1(t.fTPCnclsIter1)
  , fTPCnclsFIter1(t.fTPCnclsFIter1)
  , fITSncls(t.fITSncls)
  , fITSClusterMap(t.fITSClusterMap)
  , fTRDncls(t.fTRDncls)
  , fTRDncls0(t.fTRDncls0)
{
  // copy constructor
}

AliHLTOnlineESDtrack::~AliHLTOnlineESDtrack()
{
  // destructor
  if (fCp) delete fCp; fCp=NULL;
  if (fIp) delete fIp; fIp=NULL;
  if (fTPCInner) delete fTPCInner; fTPCInner=NULL;
  if (fOp) delete fOp; fOp=NULL;
}

AliHLTOnlineESDtrack& AliHLTOnlineESDtrack::operator=(const AliHLTOnlineESDtrack& t)
{
  // assignment operator
  if (this==&t) return *this;

  AliExternalTrackParam::operator=(t);

  return *this;
}

AliHLTOnlineESDtrack& AliHLTOnlineESDtrack::operator=(const AliESDtrack& t)
{
  // assignment operator from AliESDtrack

  AliExternalTrackParam::operator=(t);

  CopyInternalParam(fCp, t.GetConstrainedParam());
  CopyInternalParam(fIp, t.GetInnerParam());
  CopyInternalParam(fTPCInner, t.GetTPCInnerParam());
  CopyInternalParam(fOp, t.GetOuterParam());

  fFlags=t.GetStatus();
  fID=t.GetID();
  fLabel=t.GetLabel();
  fITSLabel=t.GetITSLabel();
  fTPCLabel=t.GetTPCLabel();
  fTRDLabel=t.GetTRDLabel();
  fTrackLength=t.GetIntegratedLength();

  Float_t p[2]; Float_t cov[3];

  // copy impact parameters for TPC
  t.GetImpactParametersTPC(p, cov);
  fdTPC=p[0]; fzTPC=p[1]; fCddTPC=cov[0]; fCdzTPC=cov[1]; fCzzTPC=cov[2];
  fCchi2TPC=t.GetConstrainedChi2TPC();

  // copy impact parameters
  t.GetImpactParameters(p, cov);
  fD=p[0]; fZ=p[1]; fCdd=cov[0]; fCdz=cov[1]; fCzz=cov[2];
  fCchi2=t.GetConstrainedChi2();

  fITSchi2=t.GetITSchi2();
  fTPCchi2=t.GetTPCchi2();
  fTPCchi2Iter1=t.GetTPCchi2Iter1();
  fTPCncls=t.GetTPCNcls();
  fTPCnclsF=t.GetTPCNclsF();
  fTPCnclsIter1=t.GetTPCNclsIter1();
  fTPCnclsFIter1=t.GetTPCNclsFIter1();
  fITSncls=t.GetNcls(0);
  fITSClusterMap=t.GetITSClusterMap();
  fTRDncls=t.GetTRDncls();
  fTRDncls0=t.GetTRDncls0();

  return *this;
}

void AliHLTOnlineESDtrack::CopyInternalParam(AliExternalTrackParam* &internalParam, const AliExternalTrackParam* pSrc)
{
  // copy one of the internal AliExternalTrackParam members
  if (pSrc) {
    if (!internalParam) internalParam=new AliExternalTrackParam(*pSrc);
    else (*internalParam)=(*pSrc);
  } else if (internalParam) {
    internalParam->Reset();
  }
}

void AliHLTOnlineESDtrack::Print(const char* options) const
{
  /// overloaded from TObject, print info
  AliExternalTrackParam::Print(options);
}

void AliHLTOnlineESDtrack::Dump() const
{
  /// overloaded from TObject, more crude data dump
  AliExternalTrackParam::Dump();
}

void AliHLTOnlineESDtrack::Clear(Option_t * option)
{
  /// overloaded from TObject, clear object
  
  AliExternalTrackParam::Clear(option);
}

TObject * AliHLTOnlineESDtrack::Clone(const char */*newname*/) const
{
  /// overloaded from TObject, clone object

  AliESDtrack* track=new AliESDtrack(this);

  return track;
}

void AliHLTOnlineESDtrack::Copy(TObject &object) const
{
  /// overloaded from TObject, copy object

  AliESDtrack* pESDTrack=dynamic_cast<AliESDtrack*>(&object);
  if (pESDTrack) {
    pESDTrack->Set(GetX(), GetAlpha(), GetParameter(), GetCovariance());
  }
  AliExternalTrackParam::Copy(object);
}
