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

// $Id$
// $MpId: AliMpConnection.cxx,v 1.6 2006/03/17 11:38:06 ivana Exp $
// Category: motif
//
// Class AliMpConnection
// ----------------
// Class that defines a connexion properties.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpConnection.h"

ClassImp(AliMpConnection)

//_____________________________________________________________________________
AliMpConnection::AliMpConnection(Int_t padNum, Int_t bergNum,Int_t kaptonNum,
		                 Int_t gassiNum) 
  : TObject(),
    fPadNum(padNum),
    fBergNum(bergNum),
    fKaptonNum(kaptonNum),
    fGassiNum(gassiNum),
    fOwner(0)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpConnection::AliMpConnection(const AliMpConnection& right) 
  : TObject(right) 
{
/// Protected copy constructor (not provided) 

  Fatal("AliMpConnection", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpConnection::AliMpConnection() 
  : TObject(),
    fPadNum(-1),
    fBergNum(-1),
    fKaptonNum(-1),
    fGassiNum(-1),
    fOwner(0)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpConnection::~AliMpConnection() {
//  
}

//
// operators
//

//_____________________________________________________________________________
AliMpConnection& 
AliMpConnection::operator=(const AliMpConnection& right)
{
/// Protected assignment operator (not provided)

  // check assignment to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignment operator not provided.");
    
  return *this;  
}    

