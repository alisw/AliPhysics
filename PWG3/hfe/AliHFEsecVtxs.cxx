 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//  Secondary vertexing container to store secondary vertex characteristics of 
//  2 or 3 particle sec vertex
//  from example, qusi-invariant mass, signed Lxy are stored
//
//  Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//

#include "AliLog.h"
#include "AliHFEsecVtxs.h"

ClassImp(AliHFEsecVtxs)

//_______________________________________________________________________________________________
AliHFEsecVtxs::AliHFEsecVtxs():
  fTrkLabel1(0)
  ,fTrkLabel2(0)
  ,fMCCode(0)
  ,fInvmass(0)
  ,fKFChi2(0)
  ,fSignedLxy(0)
  ,fSignedLxy2(0)
  ,fKFIP(0)
  ,fKFIP2(0)
{ 
  //
  // Default constructor
  //
}

//_______________________________________________________________________________________________
AliHFEsecVtxs::AliHFEsecVtxs(const AliHFEsecVtxs &p):
  TObject(p)
  ,fTrkLabel1(p.fTrkLabel1)
  ,fTrkLabel2(p.fTrkLabel2)
  ,fMCCode(p.fMCCode)
  ,fInvmass(p.fInvmass)
  ,fKFChi2(p.fKFChi2)
  ,fSignedLxy(p.fSignedLxy)
  ,fSignedLxy2(p.fSignedLxy2)
  ,fKFIP(p.fKFIP)
  ,fKFIP2(p.fKFIP2)
{ 
  //
  // Copy constructor
  //
}

//_______________________________________________________________________________________________
AliHFEsecVtxs&
AliHFEsecVtxs::operator=(const AliHFEsecVtxs &)
{
  // Assignment operator
  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEsecVtxs::~AliHFEsecVtxs()
{
  // Destructor
  //cout << "Analysis Done." << endl;
}
