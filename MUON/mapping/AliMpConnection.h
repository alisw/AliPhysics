/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpConnection.h,v 1.9 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpConnection
/// \brief A connection properties.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_CONNECTION_H
#define ALI_MP_CONNECTION_H

#include <TObject.h>

#include "AliMpMotifType.h"
#include "AliMpIntPair.h"

#include <TString.h>

class AliMpConnection : public TObject
{
  public:
    AliMpConnection();
    AliMpConnection(Int_t padNum,Int_t bergNum,Int_t kaptonNum,Int_t gassiNum);
    virtual ~AliMpConnection();

    // methods

    // accessors
    Int_t GetBergNum()   const {return fBergNum;}
    Int_t GetKaptonNum() const {return fKaptonNum;}
    Int_t GetGassiNum()  const {return fGassiNum;}
    Int_t GetPadNum()  const {return fPadNum;}
    AliMpMotifType *GetOwner() const {return fOwner;}
    
    void SetGassiNum(Int_t n) { fGassiNum = n; }
    
    AliMpIntPair LocalIndices() const;
    TString  PadName() const;
    // modifiers
    void SetOwner(AliMpMotifType *owner) {fOwner=owner;}

  private:
    AliMpConnection(const AliMpConnection& right);
    AliMpConnection& operator=(const AliMpConnection& right);

    // data members
    Int_t fPadNum;    ///< Pad number
    Int_t fBergNum;   ///< Berg connector number
    Int_t fKaptonNum; ///< Kapton connector number
    Int_t fGassiNum;  ///< Gassiplex channel number
    AliMpMotifType *fOwner; ///< The motif type which contains this connection

  ClassDef(AliMpConnection,1)  // Connection description
};

// inline functions

inline TString AliMpConnection::PadName() const 
{ return fOwner->PadName(fPadNum); }

inline AliMpIntPair AliMpConnection::LocalIndices() const
{ return fOwner->FindLocalIndicesByConnection(this);}

#endif //ALI_MP_CONNECTION_H
