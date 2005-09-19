/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVIndexed.h,v 1.5 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \class AliMpVIndexed
/// \brief Base class that defines the limits of global pad indices.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_INDEXED_H
#define ALI_MP_V_INDEXED_H

#include <TObject.h>

#include "AliMpIntPair.h"

class AliMpVPadIterator;

class AliMpVIndexed : public TObject
{
 public:
  AliMpVIndexed(const AliMpIntPair& lowLimit, const AliMpIntPair& highLimit);
  AliMpVIndexed();
  virtual ~AliMpVIndexed();

  // methods
  virtual AliMpVPadIterator* CreateIterator() const = 0;
  virtual AliMpIntPair  GlobalIndices(const AliMpIntPair& localIndices) const;

  // set methods
  void SetLowIndicesLimit(const AliMpIntPair& limit);
  void SetHighIndicesLimit(const AliMpIntPair& limit);

  // get methods
  Bool_t    HasIndices(const AliMpIntPair& indices) const;
  Bool_t    HasValidIndices() const;
  AliMpIntPair  GetLowIndicesLimit() const;
  AliMpIntPair  GetHighIndicesLimit() const;

 private:
  // data members 
  AliMpIntPair   fLowIndicesLimit;  // the lowest global pad indices 
  AliMpIntPair   fHighIndicesLimit; // the hihgest global pad indices 

  ClassDef(AliMpVIndexed,1) // A motif position
};

// inline functions

inline void AliMpVIndexed::SetLowIndicesLimit(const AliMpIntPair& limit)
{ fLowIndicesLimit = limit; }
  
inline void AliMpVIndexed::SetHighIndicesLimit(const AliMpIntPair& limit)
{ fHighIndicesLimit = limit; }  
  
inline AliMpIntPair AliMpVIndexed::GetLowIndicesLimit() const
{ return fLowIndicesLimit; }

inline AliMpIntPair AliMpVIndexed::GetHighIndicesLimit() const
{ return fHighIndicesLimit; }

#endif //ALI_MP_V_INDEXED_H
