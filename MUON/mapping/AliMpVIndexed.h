/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVIndexed.h,v 1.7 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpVIndexed
/// \brief Base class that defines the limits of global pad indices.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_INDEXED_H
#define ALI_MP_V_INDEXED_H

#include <TObject.h>

#include "AliMpEncodePair.h"

class AliMpVPadIterator;

class AliMpVIndexed : public TObject
{
 public:
  AliMpVIndexed();
  virtual ~AliMpVIndexed();

  // methods
          /// Create iterator over this element
  virtual AliMpVPadIterator* CreateIterator() const = 0;

  virtual MpPair_t GlobalIndices(MpPair_t localIndices) const;
  virtual Int_t  GlobalIx(Int_t localIx) const;
  virtual Int_t  GlobalIy(Int_t localIy) const;

  // set methods
  void SetLowIndicesLimit(MpPair_t limit, Bool_t valid = true);
  void SetLowIndicesLimit(Int_t ix, Int_t iy, Bool_t valid = true);

  void SetHighIndicesLimit(MpPair_t limit, Bool_t valid = true);
  void SetHighIndicesLimit(Int_t ix, Int_t iy, Bool_t valid = true);

  // get methods
  Bool_t  HasIndices(MpPair_t indices) const;
  Bool_t  HasIndices(Int_t ix, Int_t iy) const;
  Bool_t  HasValidIndices() const;

  MpPair_t  GetLowIndicesLimit() const;
  Int_t     GetLowLimitIx() const;
  Int_t     GetLowLimitIy() const;
  Bool_t    IsLowLimitValid() const;

  MpPair_t  GetHighIndicesLimit() const;
  Int_t     GetHighLimitIx() const;
  Int_t     GetHighLimitIy() const;
  Bool_t    IsHighLimitValid() const;


 private:
  // data members 
  MpPair_t  fLowLimit;  ///<  the lowest global pad indices 
  MpPair_t  fHighLimit; ///<  the highest global pad indices 
  Bool_t    fLowValid;  ///<  true, if low indices limit is set
  Bool_t    fHighValid; ///<  true, if high indices imit is set

  ClassDef(AliMpVIndexed,2) // A motif position
};


#endif //ALI_MP_V_INDEXED_H
