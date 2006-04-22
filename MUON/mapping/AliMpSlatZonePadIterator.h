#ifndef ALI_MP_SLAT_ZONE_PAD_ITERATOR_H
#define ALI_MP_SLAT_ZONE_PAD_ITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSlatZonePadIterator.h,v 1.3 2005/08/26 15:42:32 ivana Exp $

/// \ingroup slat
/// \class AliMpSlatZonePadIterator
/// \brief Iterates over slat pads within a region of constant pad size.
/// \author Laurent Aphecetche

#include "AliMpVPadIterator.h"
#include "AliMpArea.h"
#include "TVector2.h"

class AliMpSlat;
class AliMpSlatSegmentation;

class AliMpSlatZonePadIterator : public AliMpVPadIterator
{
 public:
  AliMpSlatZonePadIterator(const AliMpSlat* slat, const AliMpArea& area);
  virtual ~AliMpSlatZonePadIterator();

  void First();
  void Next();
  Bool_t IsDone() const;
  AliMpPad CurrentItem() const;
  void Invalidate();
 
 protected:
  AliMpSlatZonePadIterator(const AliMpSlatZonePadIterator& right);
  AliMpSlatZonePadIterator&  operator = (const AliMpSlatZonePadIterator& right);
     
 private:
  Bool_t CropArea();
  Bool_t GetNextPosition(Double_t& x, Double_t& y);
  void SetPad(AliMpPad& pad, const TVector2& pos);

 private:
  const AliMpSlat*       fkSlat; //! the slat we're iterating over
  AliMpSlatSegmentation* fSlatSegmentation; //! segmentation pointer
  AliMpArea  fArea; //! area we're iterating over
  TVector2   fOffset; //! current position (relative to bottom-left of area)
  TVector2   fStep; //! step sizes
  AliMpPad   fCurrentPad; //! current pad
  Bool_t     fIsDone; //! whether we've finished or not

  static const Double_t fgkDmax; // maximum double
  static const Double_t fgkEpsilon; // comparison precision 

  ClassDef(AliMpSlatZonePadIterator,1) // Pad iterator for a zone of constant density, for St345.
};

#endif
