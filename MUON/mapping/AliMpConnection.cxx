// $Id$
// Category: motif
//
// Class AliMpConnection
// ----------------
// Class that defines a connexion properties.
//
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
//
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
