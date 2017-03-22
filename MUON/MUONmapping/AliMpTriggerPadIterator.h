/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpTriggerPadIterator.h,v 1.8 2006/05/24 13:58:24 ivana Exp $

/// \ingroup slat
/// \class AliMpTriggerPadIterator
/// \brief Iterator for trigger slat pads.
/// 
//  Author: Diefo Stocco

#ifndef ALI_MP_TRIGGER_PAD_ITERATOR_H
#define ALI_MP_TRIGGER_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "TObjArray.h"

class AliMpTrigger;

class AliMpTriggerPadIterator : public AliMpVPadIterator
{
 public:
  AliMpTriggerPadIterator(); 
  // Area position must be relative to bottom-left of slat.
  AliMpTriggerPadIterator(const AliMpTrigger* slat);
  virtual ~AliMpTriggerPadIterator();

  void First();
  void Next();
  Bool_t IsDone() const;
  AliMpPad CurrentItem() const;
  void Invalidate();
 
 private:
  /// Not implemented
  AliMpTriggerPadIterator(const AliMpTriggerPadIterator&);
  /// Not implemented
  AliMpTriggerPadIterator& operator=(const AliMpTriggerPadIterator&);

  Bool_t Prepare();

 private:
  const AliMpTrigger* fkTriggerSlat; //!<! pointer to the trigger slat being iterated over
  TObjArray fPadList; //!<! Array of pads
  Int_t fCurrentPadIndex; //!<! current pad index

  ClassDef(AliMpTriggerPadIterator,0) // Pad iterator for Trigger slats
};

#endif
