/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONPainterColorSlider.h"
#include "AliMUONPainterHelper.h"
#include "AliLog.h"
#include <TGNumberEntry.h>
#include <TGButton.h>
#include <TMath.h>

///\class AliMUONPainterColorSlider
///
/// A painter color palette
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterColorSlider)
///\endcond

//_____________________________________________________________________________
AliMUONPainterColorSlider::AliMUONPainterColorSlider(const TGWindow* p, 
                                                     UInt_t w, UInt_t h)
: TGCompositeFrame(p,w,h,kVerticalFrame),
  fEntryMin(0x0),
  fEntryMax(0x0),
  fMin(FLT_MAX),
  fMax(-FLT_MAX),
  fAutoButton(new TGTextButton(this,"Auto")),
  fLockButton(new TGTextButton(this,"Lock")),
  fDefaultButton(new TGTextButton(this,"Back to default")),
  fSetDefaultButton(new TGTextButton(this,"Set as default"))
{
    /// ctor
    Int_t ndivisions(20);
  
    Int_t hsize = (h-100)/(ndivisions+2);
  Int_t topBorder(5);
  
  Double_t min(0.0);
  Double_t max(1.0);
  
  Double_t step = (max-min)/ndivisions;
  
  for ( Int_t i = -1; i < ndivisions+1; ++i ) 
  {
    Double_t value = max - (min + step*i);
    
    Int_t color = AliMUONPainterHelper::Instance()->ColorFromValue(value,
                                                                   min,max);
    Pixel_t pixel = gVirtualX->GetPixel(color);
    TGVerticalFrame* frame = new TGVerticalFrame(this,w,hsize,kFixedSize,pixel);
    
    AddFrame(frame,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));
    
    topBorder = 0;
  }
  
  fEntryMax = new TGNumberEntry(this);
  
  AddFrame(fEntryMax,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));
    
  fEntryMin = new TGNumberEntry(this);
  
  AddFrame(fEntryMin,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));
  
  AddFrame(fAutoButton,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));
  
  fAutoButton->Connect("Clicked()","AliMUONPainterColorSlider",this,"DataRangeAutoRequested()");
  
  AddFrame(fLockButton,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));
  
  fLockButton->Connect("Clicked()","AliMUONPainterColorSlider",this,"LockButtonWasClicked()");

  AddFrame(fDefaultButton,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));

  fDefaultButton->Connect("Clicked()","AliMUONPainterColorSlider",this,"DefaultButtonWasClicked()");

  AddFrame(fSetDefaultButton,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));

  fSetDefaultButton->Connect("Clicked()","AliMUONPainterColorSlider",this,"SetDefaultButtonWasClicked(Double_t*)");

  fEntryMax->Connect("ValueSet(Long_t)","AliMUONPainterColorSlider",this,"DataRangeWasChanged(Double_t*)");
  fEntryMin->Connect("ValueSet(Long_t)","AliMUONPainterColorSlider",this,"DataRangeWasChanged(Double_t*)");
}

//_____________________________________________________________________________
AliMUONPainterColorSlider::~AliMUONPainterColorSlider()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONPainterColorSlider::DataRangeAutoRequested()
{
  /// Signal that the "Auto" button was clicked
  Emit("DataRangeAutoRequested()");
}

//_____________________________________________________________________________
void 
AliMUONPainterColorSlider::DataRangeWasChanged(Double_t*)
{
  /// Data range was changed
  
  Double_t values[] = { fEntryMin->GetNumber(), fEntryMax->GetNumber() };
  
  Long_t param[] = { (Long_t)values };
  
  Emit("DataRangeWasChanged(Double_t*)",param);
}

//_____________________________________________________________________________
void AliMUONPainterColorSlider::DefaultButtonWasClicked()
{
  /// Signal that the "Default" button was clicked
  Emit("DefaultButtonWasClicked()");
}

//_____________________________________________________________________________
void AliMUONPainterColorSlider::SetDefaultButtonWasClicked(Double_t*)
{
  /// Signal that the "SetDefault" button was clicked
  
  Double_t values[] = { fEntryMin->GetNumber(), fEntryMax->GetNumber() };
  
  Long_t param[] = { (Long_t)values };
  
  Emit("SetDefaultButtonWasClicked(Double_t*)",param);
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterColorSlider::IsLocked() const
{
  /// Whether our range is locked or not

  return (fLockButton->GetString() == "Unlock"); // if we can unlock it means we are locked...
}

//_____________________________________________________________________________
void AliMUONPainterColorSlider::LockDefaultButtons()
{
  fDefaultButton->SetEnabled(kFALSE);
  fSetDefaultButton->SetEnabled(kFALSE);
}

//_____________________________________________________________________________
void AliMUONPainterColorSlider::UnlockDefaultButtons()
{
  fDefaultButton->SetEnabled(kTRUE);
  fSetDefaultButton->SetEnabled(kTRUE);
}

//_____________________________________________________________________________
void 
AliMUONPainterColorSlider::LockButtonWasClicked()
{
  /// Lock (toggle button) was clicked
  
  if ( IsLocked() )
  {
    // unlock it
    fLockButton->SetText("Lock");
    fEntryMin->SetState(kTRUE);
    fEntryMax->SetState(kTRUE);
    fAutoButton->SetEnabled(kTRUE);
    UnlockDefaultButtons();
  }
  else
  {
    // lock it
    fLockButton->SetText("Unlock");
    fEntryMin->SetState(kFALSE);
    fEntryMax->SetState(kFALSE);
    fAutoButton->SetEnabled(kFALSE);
    LockDefaultButtons();
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterColorSlider::SetRange(Double_t min, Double_t max, Bool_t emit)
{
  /// Set the data range
  
  if ( !IsLocked() )
  {
    fMin = min;
    fMax = max;
  
    fEntryMin->SetNumber(fMin);
    fEntryMax->SetNumber(fMax);
  }
  
  if ( emit ) 
  {
    DataRangeWasChanged(0x0);
  }
}

