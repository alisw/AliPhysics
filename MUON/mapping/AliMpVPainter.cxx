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
// $MpId: AliMpVPainter.cxx,v 1.10 2006/05/24 13:58:32 ivana Exp $
// Category: graphics
//
// Class AliMpVPainter
// --------------
// Class for drawing objects into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
  
#include "AliMpVPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpSector.h"
#include "AliMpRow.h"
#include "AliMpZone.h"
#include "AliMpSubZone.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifPosition.h"
#include "AliMpSectorPainter.h"
#include "AliMpRowPainter.h"
#include "AliMpZonePainter.h"
#include "AliMpSubZonePainter.h"
#include "AliMpRowSegmentPainter.h"
#include "AliMpMotifPainter.h"
#include "AliMpMotifType.h"
#include "AliMpPCB.h"
#include "AliMpPCBPainter.h"
#include "AliMpSlat.h"
#include "AliMpSlatPainter.h"

#include <TList.h>
#include <TVirtualX.h>
#include <TPad.h>

/// \cond CLASSIMP
ClassImp(AliMpVPainter)
/// \endcond

//_______________________________________________________________________
AliMpVPainter::AliMpVPainter()
  : TObject(),
    fColor(2),
    fPadPosition(),
    fPadDimensions(),
    fTrashList(0)
{
  /// Default constructor

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  fPadPosition  =  gr->GetPadPosition();
  fPadDimensions =  gr->GetPadDimensions();
  fColor=gr->GetColor();
  fTrashList = new TList;
}

//_______________________________________________________________________
AliMpVPainter::~AliMpVPainter()
{
  /// Destructor

  if (fTrashList){
    fTrashList->Delete();
    delete fTrashList;
  }
}

//_______________________________________________________________________
Bool_t AliMpVPainter::IsInside(const TVector2 &point,const TVector2& pos,const TVector2& dim)
{
  /// Is the point \a point inside the \a area (pos,dim)?

  return ( (TMath::Abs(point.X()-pos.X())<dim.X() ) && (TMath::Abs(point.Y()-pos.Y())<dim.Y() ) );
}

//_______________________________________________________________________
Int_t AliMpVPainter::DistancetoPrimitive(Int_t x, Int_t y)
{
  /// Distance to the center if (x,y) is inside the box defined by (fPadPosition,fPadDimensions)
  /// 9999 otherwise

  TVector2 point = TVector2(gPad->AbsPixeltoX(x), gPad->AbsPixeltoY(y));
  if ( IsInside(point,fPadPosition,fPadDimensions) )
    {
      return (Int_t)(point-fPadPosition).Mod();
    }
  return 9999;
}

//_______________________________________________________________________
void AliMpVPainter::DumpObject() const
{
  /// Dump the painted object
}

//_______________________________________________________________________
TObject* AliMpVPainter::Clone(const char* newname) const
{
  /// Create a clone of this object

  AliMpVPainter *newobj = (AliMpVPainter *)TObject::Clone(newname);
  if (!newobj) return 0;
  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  newobj->fPadPosition  =  gr->GetPadPosition();
  newobj->fPadDimensions =  gr->GetPadDimensions();
  return newobj;
}

//_______________________________________________________________________
TObject* AliMpVPainter::DrawClone(Option_t* option) const
{
  /// Draw the clone object

  TVirtualPad *pad = gROOT->GetSelectedPad();
  TVirtualPad *padsav = gPad;

  TObject *newobj = Clone();

  if (!newobj) return 0;

  if (pad) pad->cd();
  newobj->Draw(option);
  if (padsav) padsav->cd();
  return newobj;
}

//_______________________________________________________________________
AliMpVPainter *AliMpVPainter::CreatePainter(TObject *object)
{
  /// Create a new painter, which correspond to the
  /// class of object

  AliMpVPainter *painter=0;
  if (object->InheritsFrom(AliMpSector::Class()))
    painter = new AliMpSectorPainter((AliMpSector *)object);
  else if (object->InheritsFrom(AliMpZone::Class()))
    painter = new AliMpZonePainter((AliMpZone *)object);
  else if (object->InheritsFrom(AliMpSubZone::Class()))
    painter = new AliMpSubZonePainter((AliMpSubZone *)object);
  else if (object->InheritsFrom(AliMpRow::Class()))
    painter = new AliMpRowPainter((AliMpRow *)object);
  else if (object->InheritsFrom(AliMpVRowSegment::Class()))
    painter = new AliMpRowSegmentPainter((AliMpVRowSegment *)object);
  else if (object->InheritsFrom(AliMpMotifPosition::Class()))
    painter = new AliMpMotifPainter((AliMpMotifPosition *)object);
  else if (object->InheritsFrom(AliMpMotifType::Class()))
    painter = new AliMpMotifPainter((AliMpMotifType *)object);
  else if (object->InheritsFrom(AliMpPCB::Class()))
    painter = new AliMpPCBPainter((AliMpPCB *)object);
  else if (object->InheritsFrom(AliMpSlat::Class()))
    painter = new AliMpSlatPainter((AliMpSlat*)object);
  return painter;
}

//_______________________________________________________________________
void AliMpVPainter::AddPainter(AliMpVPainter *painter)
{
  /// Add a painter to the list of painters (private)

  fTrashList->Add(painter);
}


//_______________________________________________________________________
AliMpVPainter *AliMpVPainter::DrawObject(TObject *object,Option_t *option)
{
  /// Draw the object                                                  \n
  /// Return the AliMpVPainter object created for the drawing

  AliMpVPainter *painter=CreatePainter(object);

  if (painter){
    painter->Draw(option);
    AddPainter(painter);
  }
  return painter;
}

//_______________________________________________________________________
void AliMpVPainter::InitGraphContext()
{
  /// Set the pad and real area of the graphic context to
  /// the one stored in this painter

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  gr->SetPadPosition(fPadPosition);
  gr->SetPadDimensions(fPadDimensions);
  gr->SetRealPosition(GetPosition());
  gr->SetRealDimensions(GetDimensions());
  gVirtualX->SetFillColor(fColor);
  gVirtualX->SetTextColor(1);
  gVirtualX->SetLineColor(1);
}

//_______________________________________________________________________
void AliMpVPainter::PaintWholeBox(Bool_t fill)
{
  /// Paint the box around the total pad area given in this painter
  /// fill it or bnot following the parameter value

  Double_t x1,y1,x2,y2;
  x1 = fPadPosition.X()-fPadDimensions.X();
  y1 = fPadPosition.Y()-fPadDimensions.Y();
  x2 = fPadPosition.X()+fPadDimensions.X();
  y2 = fPadPosition.Y()+fPadDimensions.Y();

  Style_t sty = gVirtualX->GetFillStyle();
  gVirtualX->SetFillStyle(fill?1:0);
  gPad->PaintBox(x1,y1,x2,y2);
  gVirtualX->SetFillStyle(0);
  gPad->PaintBox(x1,y1,x2,y2);
  gVirtualX->SetFillStyle(sty);
}

//_______________________________________________________________________
TVector2 AliMpVPainter::RealToPad(const TVector2& realPos)
{
  /// Transform a real position into its equivalent position in a canvas

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  gr->Push();
  gr->SetPadPosition(fPadPosition);
  gr->SetPadDimensions(fPadDimensions);
  gr->SetRealPosition(GetPosition());
  gr->SetRealDimensions(GetDimensions());


  TVector2 ans = gr->RealToPad(realPos);
  gr->Pop();
  return ans;
}
