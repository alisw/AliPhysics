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

/// @file   AliHLTESDtrack.cxx
/// @author Matthias Richter
/// @date   2010-10-29
/// @brief  An AliESDtrack child class doing the conversion between 
///         AliHLTESDOptTrack and AliESDtrack
/// @note   

#include "AliHLTESDtrack.h"
#include "AliHLTOnlineESDtrack.h"
#include "AliESDtrack.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTESDtrack)

AliHLTESDtrack::AliHLTESDtrack()
  : AliESDtrack()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTESDtrack::AliHLTESDtrack(const AliHLTESDtrack& t)
  : AliESDtrack(t)
{
  // copy constructor
}

AliHLTESDtrack::~AliHLTESDtrack()
{
  // destructor
}

AliHLTESDtrack& AliHLTESDtrack::operator=(const AliHLTESDtrack& t)
{
  // assignment operator
  if (this==&t) return *this;

  AliExternalTrackParam::operator=(t);

  return *this;
}

AliHLTESDtrack& AliHLTESDtrack::operator=(const AliHLTOnlineESDtrack& t)
{
  // assignment operator from AliESDtrack

  AliHLTOnlineESDtrack::CopyInternalParam(fCp, t.GetConstrainedParam());
  AliHLTOnlineESDtrack::CopyInternalParam(fIp, t.GetInnerParam());
  AliHLTOnlineESDtrack::CopyInternalParam(fTPCInner, t.GetTPCInnerParam());
  AliHLTOnlineESDtrack::CopyInternalParam(fOp, t.GetOuterParam());

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
  fITSncls=t.GetITSNcls();
  fITSClusterMap=t.GetITSClusterMap();
  fTRDncls=t.GetTRDncls();
  fTRDncls0=t.GetTRDncls0();

  return *this;
}
