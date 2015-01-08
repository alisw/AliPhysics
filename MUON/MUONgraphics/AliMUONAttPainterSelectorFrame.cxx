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

///\class AliMUONAttPainterSelectorFrame
///
/// Widget to select painter view type
///
///\author Laurent Aphecetche, Subatech

#include "AliMUONAttPainterSelectorFrame.h"

#include "AliMUONAttPainter.h"
#include "AliMUONPainterInterfaceHelper.h"
#include "AliLog.h"
#include <TGButton.h>
#include <TGButtonGroup.h>

///\cond CLASSIMP
ClassImp(AliMUONAttPainterSelectorFrame)
///\endcond

//_____________________________________________________________________________
AliMUONAttPainterSelectorFrame::AliMUONAttPainterSelectorFrame(TGWindow* p, UInt_t w, UInt_t h)
: TGHorizontalFrame(p,w,h),
fCathode(0x0),
fPlane(0x0),
fViewPoint(0x0),
fAttributes()
{
  /// ctor
  
  fCathode = new TGButtonGroup(this,"Cathode",kHorizontalFrame);

  fAttributes.SetCathode(kTRUE,kFALSE);
  AliMUONPainterInterfaceHelper::AddRadioButton(*fCathode,fAttributes.CathodeName(),(void*)(10));
  fAttributes.SetCathode(kFALSE,kTRUE);
  AliMUONPainterInterfaceHelper::AddRadioButton(*fCathode,fAttributes.CathodeName(),(void*)(1));
                                                
  fPlane = new TGButtonGroup(this,"Plane",kHorizontalFrame);
  
  fAttributes.SetPlane(kTRUE,kFALSE);
  AliMUONPainterInterfaceHelper::AddRadioButton(*fPlane,fAttributes.PlaneName(),(void*)(10));
  fAttributes.SetPlane(kFALSE,kTRUE);
  AliMUONPainterInterfaceHelper::AddRadioButton(*fPlane,fAttributes.PlaneName(),(void*)(1));

  fViewPoint = new TGButtonGroup(this,"ViewPoint",kHorizontalFrame);

  fAttributes.SetViewPoint(kTRUE,kFALSE);
  AliMUONPainterInterfaceHelper::AddRadioButton(*fViewPoint,fAttributes.ViewPointName(),(void*)(10));
  fAttributes.SetViewPoint(kFALSE,kTRUE);
  AliMUONPainterInterfaceHelper::AddRadioButton(*fViewPoint,fAttributes.ViewPointName(),(void*)(1));

  fViewPoint->SetState(kFALSE); //FIXME: until we're sure back views are handled correctly
  
  AddFrame(fCathode);
  AddFrame(fPlane);
  AddFrame(fViewPoint);
  
  fCathode->Connect("Clicked(Int_t)","AliMUONAttPainterSelectorFrame",this,"CathodeClicked(Int_t)");
  fPlane->Connect("Clicked(Int_t)","AliMUONAttPainterSelectorFrame",this,"PlaneClicked(Int_t)");
  fViewPoint->Connect("Clicked(Int_t)","AliMUONAttPainterSelectorFrame",this,"ViewClicked(Int_t)");
}

//_____________________________________________________________________________
AliMUONAttPainterSelectorFrame::~AliMUONAttPainterSelectorFrame()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONAttPainterSelectorFrame::CathodeClicked(Int_t buttonId)
{
  /// Cathode button clicked
  
  fAttributes.SetPlane(kFALSE,kFALSE);

  TGButton* button = fCathode->GetButton(buttonId);
  
  Long_t i = reinterpret_cast<Long_t>(button->GetUserData());
  
  if ( i == 10 ) 
  {
    fAttributes.SetCathode(kTRUE,kFALSE);
  }
  else if ( i == 1 ) 
  {
    fAttributes.SetCathode(kFALSE,kTRUE);
  }
  else
  {
    AliFatal("");
  }
  
  Clicked(&fAttributes);
}

//_____________________________________________________________________________
void
AliMUONAttPainterSelectorFrame::PlaneClicked(Int_t buttonId)
{
  /// Plane button clicked
  
  fAttributes.SetCathode(kFALSE,kFALSE);
  
  TGButton* button = fPlane->GetButton(buttonId);
  
  Long_t i = reinterpret_cast<Long_t> (button->GetUserData());
  
  if ( i == 10 ) 
  {
    fAttributes.SetPlane(kTRUE,kFALSE);
  }
  else if ( i == 1 ) 
  {
    fAttributes.SetPlane(kFALSE,kTRUE);
  }
  else
  {
    AliFatal("");
  }
  
  Clicked(&fAttributes);
  
}

//_____________________________________________________________________________
void
AliMUONAttPainterSelectorFrame::ViewClicked(Int_t buttonId)
{
  /// View button clicked

  TGButton* button = fViewPoint->GetButton(buttonId);
  
  Long_t i = reinterpret_cast<Long_t> (button->GetUserData());
  
  if ( i == 10 ) 
  {
    fAttributes.SetViewPoint(kTRUE,kFALSE);
  }
  else if ( i == 1 ) 
  {
    fAttributes.SetViewPoint(kFALSE,kTRUE);
  }
  else
  {
    AliFatal("");
  }
  
  Clicked(&fAttributes);
  
}

//_____________________________________________________________________________
void
AliMUONAttPainterSelectorFrame::Clicked(const AliMUONAttPainter* newValues)
{
  /// Emit a signal
  
  Long_t params[] = { (Long_t)newValues };
  
  Emit("Clicked(AliMUONAttPainter*)",params);
}

//_____________________________________________________________________________
void
AliMUONAttPainterSelectorFrame::Update(const AliMUONAttPainter& att)
{
  /// Update button state from the painter attributes

  AliMUONPainterInterfaceHelper::Unselect(*fCathode,"*");
  AliMUONPainterInterfaceHelper::Unselect(*fPlane,"*");
  AliMUONPainterInterfaceHelper::Unselect(*fViewPoint,"*");

  fAttributes = att;
  
  fCathode->SetState(!fAttributes.IsCathodeAndPlaneDisabled());
  fPlane->SetState(!fAttributes.IsCathodeAndPlaneDisabled());

  if ( fAttributes.IsCathodeDefined() ) 
  {
    AliMUONPainterInterfaceHelper::Select(*fCathode,fAttributes.CathodeName());
  }
  
  if ( fAttributes.IsPlaneDefined() ) 
  {
    AliMUONPainterInterfaceHelper::Select(*fPlane,fAttributes.PlaneName());
  }
  
  AliMUONPainterInterfaceHelper::Select(*fViewPoint,fAttributes.ViewPointName());
  
}
