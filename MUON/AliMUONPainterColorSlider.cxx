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
  fMax(-FLT_MAX)
{
    /// ctor
    Int_t ndivisions(20);
  
    Int_t hsize = (h-100)/(ndivisions+2);
  Int_t topBorder(5);
  
  Double_t min(0.0);
  Double_t max(1.0);
  
  Double_t step = (max-min)/ndivisions;
  
  for ( Int_t i = -1; i <= ndivisions+1; ++i ) 
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
  
//  fEntryMin->SetFormat(TGNumberFormat::kNESRealOne);
//  fEntryMax->SetFormat(TGNumberFormat::kNESRealOne);
  
  TGTextButton* button = new TGTextButton(this,"Auto");
  
  AddFrame(button,new TGLayoutHints(kLHintsExpandX,0,0,topBorder,0));
  
  button->Connect("Clicked()","AliMUONPainterColorSlider",this,"DataRangeAutoRequested()");
  
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
  AliDebug(1,"");
  Emit("DataRangeAutoRequested()");
}

//_____________________________________________________________________________
void 
AliMUONPainterColorSlider::DataRangeWasChanged(Double_t*)
{
  /// Data range was changed
  
  Double_t values[] = { fEntryMin->GetNumber(), fEntryMax->GetNumber() };
  
  Long_t param[] = { (Long_t)values };
  
  AliDebug(1,Form("double min %e max %e",values[0],values[1]));
  
  Emit("DataRangeWasChanged(Double_t*)",param);
}

//_____________________________________________________________________________
void 
AliMUONPainterColorSlider::SetRange(Double_t min, Double_t max, Bool_t emit)
{
  /// Set the data range
  
  AliDebug(1,Form("min %e max %e emit %d",min,max,emit));
  
  fMin = min;
  fMax = max;
  
  fEntryMin->SetNumber(fMin);
  fEntryMax->SetNumber(fMax);
  
  if ( emit ) 
  {
    DataRangeWasChanged(0x0);
  }
}

