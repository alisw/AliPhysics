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
#include "AliMpEncodePair.h"

#include <TString.h>

class AliMpConnection : public TObject
{
  public:
    AliMpConnection(Int_t padNum, 
                    Int_t bergNum,
                    Int_t kaptonNum,
                    Int_t gassiNum,
                    MpPair_t localIndices);
    AliMpConnection(TRootIOCtor* /*ioCtor*/);
    //AliMpConnection();
    virtual ~AliMpConnection();

    //
    // accessors
    //
          /// Return Berg connector number
    Int_t GetBergNum()   const     { return fBergNum; }
          /// Return kapton connector number
    Int_t GetKaptonNum() const     { return fKaptonNum; }
          /// Return manu channel number
    Int_t GetManuChannel() const   { return fGassiNum; }
          /// Return pad number
    Int_t GetPadNum()  const       { return GetUniqueID(); }

          /// Return encoded local indices
    MpPair_t GetLocalIndices() const { return fLocalIndices; }
    Int_t  GetLocalIx() const;
    Int_t  GetLocalIy() const;
    
          /// Return the motif type which contains this connection
    AliMpMotifType *GetOwner() const { return fOwner; }
    
    TString  PadName() const;
    
    //
    // modifiers
    //

          /// Set Gassiplex channel number
    void SetGassiNum(Int_t n)            { fGassiNum = n; }
          /// Set the motif type which contains this connection
    void SetOwner(AliMpMotifType *owner) { fOwner=owner; }
    
  private:
    /// Not implemented
    AliMpConnection();
    /// Not implemented
    AliMpConnection(const AliMpConnection& right);
    /// Not implemented
    AliMpConnection& operator=(const AliMpConnection& right);

    // data members
    Int_t   fBergNum;   ///< Berg connector number
    Int_t   fKaptonNum; ///< Kapton connector number
    Int_t   fGassiNum;  ///< Gassiplex channel number
    MpPair_t  fLocalIndices;  ///< Local indices
    AliMpMotifType *fOwner; ///< The motif type which contains this connection

  ClassDef(AliMpConnection,2)  // Connection description
};

// inline functions

/// Return the pad number converted to a name
inline TString AliMpConnection::PadName() const 
{ return fOwner->PadName(GetUniqueID()); }

#endif //ALI_MP_CONNECTION_H
