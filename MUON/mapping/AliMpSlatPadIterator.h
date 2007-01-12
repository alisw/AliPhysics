/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSlatPadIterator.h,v 1.8 2006/05/24 13:58:24 ivana Exp $

/// \ingroup slat
/// \class AliMpSlatPadIterator
/// \brief Iterator for slat pads.
/// 
//  Author: Laurent Aphecetche

#ifndef ALI_MP_SLAT_PAD_ITERATOR_H
#define ALI_MP_SLAT_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "TObjArray.h"

class AliMpSlat;
class AliMpArea;

class AliMpSlatPadIterator : public AliMpVPadIterator
{
 public:
  AliMpSlatPadIterator(); 
  // Area position must be relative to bottom-left of slat.
  AliMpSlatPadIterator(const AliMpSlat* slat, const AliMpArea& area);
  virtual ~AliMpSlatPadIterator();

  void First();
  void Next();
  Bool_t IsDone() const;
  AliMpPad CurrentItem() const;
  void Invalidate();
 
 private:
  AliMpSlatPadIterator(const AliMpSlatPadIterator&);
	AliMpSlatPadIterator& operator=(const AliMpSlatPadIterator&);
  Bool_t Prepare(const AliMpArea& area);
  AliMpArea Intersect(const AliMpArea& a, const AliMpArea& b) const;

 private:
  const AliMpSlat* fkSlat; ///< pointer to the slat being iterated over
  TObjArray fDelegates; ///< iterators we do use (array of AliMpVPadIterator*)
  AliMpVPadIterator* fCurrentDelegate; ///< current iterator
  Int_t fCurrentDelegateIndex; ///< current iterator index

  ClassDef(AliMpSlatPadIterator,2) // Pad iterator for St 345 Slats
};

#endif
