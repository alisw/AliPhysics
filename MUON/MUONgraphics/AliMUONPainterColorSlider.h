#ifndef ALIMUONPAINTERCOLORSLIDER_H
#define ALIMUONPAINTERCOLORSLIDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterColorSlider
/// \brief A vertical color palette
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif

class TGNumberEntry;
class TGTextButton;

class AliMUONPainterColorSlider : public TGCompositeFrame
{
public:
  AliMUONPainterColorSlider(const TGWindow* p, UInt_t w, UInt_t h);
  virtual ~AliMUONPainterColorSlider();

  void DefaultButtonWasClicked(); // *SIGNAL*

  void SetDefaultButtonWasClicked(Double_t* range); // *SIGNAL*

  void DataRangeAutoRequested(); // *SIGNAL*

  void DataRangeWasChanged(Double_t* range); // *SIGNAL*

  void LockButtonWasClicked(); 
  
  void SetRange(Double_t min, Double_t max, Bool_t emit=kTRUE);
  
  Bool_t IsLocked() const;
  
private:
  /// Not implemented
  AliMUONPainterColorSlider(const AliMUONPainterColorSlider& rhs);
  /// Not implemented
  AliMUONPainterColorSlider& operator=(const AliMUONPainterColorSlider& rhs);

  void LockDefaultButtons();
  void UnlockDefaultButtons();
  
private:
  TGNumberEntry* fEntryMin; ///< textbox for min value to be represented
  TGNumberEntry* fEntryMax; ///< textbox for max value to be represented
  Double_t fMin; ///< min value to be represented
  Double_t fMax; ///< max value to be represented
  TGTextButton* fAutoButton; ///< to toggle data range computation
  TGTextButton* fLockButton; ///< to toggle locking of range
  TGTextButton* fDefaultButton; ///< to get back to default range (if defined)
  TGTextButton* fSetDefaultButton; ///< to set the default for the current source
  
  ClassDef(AliMUONPainterColorSlider,3) // A painter color palette
};

#endif
