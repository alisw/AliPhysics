// $Id$
// Category: motif
//
// Class AliMpConnection
// ----------------
// Class that defines a connexion properties.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

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
//
}

//_____________________________________________________________________________
AliMpConnection::AliMpConnection(const AliMpConnection& right) 
  : TObject(right) {
// 
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
//
}

//_____________________________________________________________________________
AliMpConnection::~AliMpConnection() {
//  
}

// operators

//_____________________________________________________________________________
AliMpConnection& 
AliMpConnection::operator=(const AliMpConnection& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

