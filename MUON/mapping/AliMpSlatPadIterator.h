/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSlatPadIterator.h,v 1.4 2005/09/19 19:01:09 ivana Exp $

/// \ingroup slat
/// \class AliMpSlatPadIterator
/// \brief Iterator for slat pads.
/// 
/// Author: Laurent Aphecetche

#ifndef ALI_MP_SLAT_PAD_ITERATOR_H
#define ALI_MP_SLAT_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include <vector>

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
  const AliMpSlat* fkSlat;
  std::vector<AliMpVPadIterator*> fDelegates;
  AliMpVPadIterator* fCurrentDelegate;
  UInt_t fCurrentDelegateIndex;

  ClassDef(AliMpSlatPadIterator,1) // Pad iterator for St 345 Slats
};

#endif
