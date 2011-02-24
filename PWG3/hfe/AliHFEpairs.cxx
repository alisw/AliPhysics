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

/* $Id$ */

//
//  Container class to store pair characteristics
//  for secondary vertex analysis
//  from example, qusi-invariant mass, signed Lxy are stored
//
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//

#include "AliLog.h"
#include "AliHFEpairs.h"

ClassImp(AliHFEpairs)

//_______________________________________________________________________________________________
AliHFEpairs::AliHFEpairs():
  fTrkLabel(0)
  ,fPairCode(0)
  ,fInvmass(0)
  ,fKFChi2(0)
  ,fOpenangle(0)
  ,fCosOpenangle(0)
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
AliHFEpairs::AliHFEpairs(const AliHFEpairs &p):
  TObject(p)
  ,fTrkLabel(p.fTrkLabel)
  ,fPairCode(p.fPairCode)
  ,fInvmass(p.fInvmass)
  ,fKFChi2(p.fKFChi2)
  ,fOpenangle(p.fOpenangle)
  ,fCosOpenangle(p.fCosOpenangle)
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
AliHFEpairs&
AliHFEpairs::operator=(const AliHFEpairs &)
{
  //
  // Assignment operator
  //

  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEpairs::~AliHFEpairs()
{
  //
  // Destructor
  //

  //cout << "Analysis Done." << endl;
}
