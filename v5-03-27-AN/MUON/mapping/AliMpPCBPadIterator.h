#ifndef ALI_MP_PCB_PAD_ITERATOR_H
#define ALI_MP_PCB_PAD_ITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId$

/// \ingroup slat
/// \class AliMpPCBPadIterator
/// \brief Iterates over slat pads within a region of constant pad size.
///
//  Author: Laurent Aphecetche

#include "AliMpVPadIterator.h"

class AliMpArea;
class AliMpSlat;
class AliMpSlatSegmentation;

class AliMpPCBPadIterator : public AliMpVPadIterator
{
public:
  AliMpPCBPadIterator(const AliMpSlat* slat, const AliMpArea& area);
  virtual ~AliMpPCBPadIterator();
  
  void First();
  void Next();
  Bool_t IsDone() const;
  AliMpPad CurrentItem() const;
  void Invalidate();
  
  void Print(Option_t* opt="") const;
  
private:
  /// Not implemented
  AliMpPCBPadIterator(const AliMpPCBPadIterator& right);
  /// Not implemented
  AliMpPCBPadIterator&  operator = (const AliMpPCBPadIterator& right);
  
  Bool_t GetNextPosition(Int_t& ix, Int_t& iy) const;
  Bool_t CropArea(const AliMpArea& area);
  void SetPad(AliMpPad& pad, Int_t ix, Int_t iy);
  
private:
  const AliMpSlat*       fkSlat; //!< the slat we're iterating over
  AliMpSlatSegmentation* fSlatSegmentation; //!< segmentation pointer
  MpPair_t fMinIndices; //!< indices of bottom left of region to iterate over
  MpPair_t fMaxIndices; //!< indices of top right of region to iterate over
  MpPair_t fOffset; //!< current position
  AliMpPad fCurrentPad; //!< current pad
  Bool_t fIsDone; //!< whether we've finished or not
  
  ClassDef(AliMpPCBPadIterator,0) // Pad iterator for a zone of constant density, for St345.
};

#endif
