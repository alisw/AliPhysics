/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpArrayI.h,v 1.4 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpArrayI
/// \brief Helper class for sorted integer array
///
/// \author Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ARRAY_I_H
#define ALI_MP_ARRAY_I_H

#include "AliMpIntPair.h"

#include <TObject.h>
#include <TArrayI.h>

class TString;

class AliMpArrayI : public TObject
{
  public:
    AliMpArrayI(Bool_t sort = true);
    AliMpArrayI(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpArrayI();
    
    // methods
    Bool_t Add(Int_t value);
    Bool_t Remove(Int_t value);

    // set methods
    void SetSize(Int_t size);

    // get methods
    Int_t   GetSize() const;
    Int_t   GetValue(Int_t index) const;
    Bool_t  HasValue(Int_t value) const;
    
  private:  
    // methods
    Int_t  GetPosition(Int_t value) const;
  
    // static data members
    static const Int_t    fgkDefaultSize; ///< Default initial size

    // data members
    Bool_t   fSort;       ///< Option to sort the values
    Int_t    fNofValues;  ///< Number of values in the array
    TArrayI  fValues;     ///< Array of values 
    AliMpIntPair fLimits; ///< The minimum and maximum values in the array

  ClassDef(AliMpArrayI,1)  // Helper class for sorted integer array
};

#endif //ALI_MP_EX_MAP_H

