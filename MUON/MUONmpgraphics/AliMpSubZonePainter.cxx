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
// $MpId: AliMpSubZonePainter.cxx,v 1.8 2006/05/24 13:58:32 ivana Exp $
// Category: graphics

//-----------------------------------------------------------------------------
// Class AliMpSubZonePainter
// -------------------------
// Class for drawing a subzone into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
//-----------------------------------------------------------------------------
  
#include "AliMpSubZonePainter.h"
#include "AliMpGraphContext.h"
#include "AliMpSubZone.h"
#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"

#include <TVirtualX.h>
#include <TPad.h>

/// \cond CLASSIMP
ClassImp(AliMpSubZonePainter)
/// \endcond

//_______________________________________________________________________
AliMpSubZonePainter::AliMpSubZonePainter()
  : AliMpVPainter(),
    fSubZone(0)
{
  /// Default constructor
}

//_______________________________________________________________________
AliMpSubZonePainter::AliMpSubZonePainter(AliMpSubZone *subZone)
  : AliMpVPainter(),
    fSubZone(subZone)
{
  /// Standard constructor 

}

//_______________________________________________________________________
AliMpSubZonePainter::~AliMpSubZonePainter()
{
  /// Destructor
}

//_______________________________________________________________________
Int_t AliMpSubZonePainter::DistancetoPrimitive(Int_t x, Int_t y)
{
  /// Dist to the nearest segment center if (x,y) is inside the sub-zone
  /// 9999 otherwise
  
  if (fSubZone->GetNofRowSegments()<1) return 9999;
  AliMpGraphContext *gr = AliMpGraphContext::Instance();

  gr->Push();
  InitGraphContext();


  TVector2 point = TVector2(gPad->AbsPixeltoX(x), gPad->AbsPixeltoY(y));

  Double_t res=9999.;
  for (Int_t iseg=0;iseg<fSubZone->GetNofRowSegments();++iseg){
    //for each row segments
    AliMpVRowSegment* seg = fSubZone->GetRowSegment(iseg);

    TVector2 pos,dim;
    gr->RealToPad(TVector2(seg->GetPositionX(), seg->GetPositionY()),
                  TVector2(seg->GetDimensionX(),seg->GetDimensionY()),pos,dim);

    if ( IsInside(point,pos,dim) ){
      Double_t value = (point-pos).Mod();
      if (value<res) res=value;
    }
  }
  gr->Pop();
  return (Int_t)res;
}

//_______________________________________________________________________
void AliMpSubZonePainter::DumpObject()
{
  /// Draw the owned object
  
  fSubZone->Dump();
}

//_______________________________________________________________________
TVector2 AliMpSubZonePainter::GetPosition() const
{
  /// Get the owned object's position

  if (fSubZone->GetNofRowSegments()<1) return TVector2(0.,0.);
  AliMpVRowSegment* seg = fSubZone->GetRowSegment(0);

  // bl = bottom left position;
  TVector2 bl = TVector2(seg->GetPositionX(), seg->GetPositionY()) -
                TVector2(seg->GetDimensionX(),seg->GetDimensionY());
  // ur = upper right position
  TVector2 ur = TVector2(seg->GetPositionX(), seg->GetPositionY()) +
                TVector2(seg->GetDimensionX(),seg->GetDimensionY());

  for (Int_t iseg=1;iseg<fSubZone->GetNofRowSegments();++iseg){
    seg = fSubZone->GetRowSegment(iseg);
    // update the bottom-left corner
    if (bl.X()>seg->GetPositionX()-seg->GetDimensionX())
      bl.Set(seg->GetPositionX()-seg->GetDimensionX(),bl.Y());
    if (bl.Y()>seg->GetPositionY()-seg->GetDimensionY())
      bl.Set(bl.X(),seg->GetPositionY()-seg->GetDimensionY());
    // update the upper-right corner
    if (ur.X()<seg->GetPositionX()+seg->GetDimensionX())
      ur.Set(seg->GetPositionX()+seg->GetDimensionX(),ur.Y());
    if (ur.Y()<seg->GetPositionY()+seg->GetDimensionY())
      ur.Set(ur.X(),seg->GetPositionY()+seg->GetDimensionY());
  }
  return (ur+bl)/2.;
}

//_______________________________________________________________________
TVector2 AliMpSubZonePainter::GetDimensions() const
{
  /// Get the owned object's dimensions

  if (fSubZone->GetNofRowSegments()<1) return TVector2(0.,0.);
  AliMpVRowSegment* seg = fSubZone->GetRowSegment(0);

  // bl = bottom left position;
  TVector2 bl = TVector2(seg->GetPositionX(), seg->GetPositionY()) - 
                TVector2(seg->GetDimensionX(),seg->GetDimensionY());
  // ur = upper right position
  TVector2 ur = TVector2(seg->GetPositionX(), seg->GetPositionY()) + 
                TVector2(seg->GetDimensionX(),seg->GetDimensionY());

  for (Int_t iseg=1;iseg<fSubZone->GetNofRowSegments();++iseg){
    seg = fSubZone->GetRowSegment(iseg);
    // update the bottom-left corner
    if (bl.X()>seg->GetPositionX()-seg->GetDimensionX())
      bl.Set(seg->GetPositionX()-seg->GetDimensionX(),bl.Y());
    if (bl.Y()>seg->GetPositionY()-seg->GetDimensionY())
      bl.Set(bl.X(),seg->GetPositionY()-seg->GetDimensionY());
    // update the upper-right corner
    if (ur.X()<seg->GetPositionX()+seg->GetDimensionX())
      ur.Set(seg->GetPositionX()+seg->GetDimensionX(),ur.Y());
    if (ur.Y()<seg->GetPositionY()+seg->GetDimensionY())
      ur.Set(ur.X(),seg->GetPositionY()+seg->GetDimensionY());
  }
  return (ur-bl)/2.;
}

//_______________________________________________________________________
void AliMpSubZonePainter::Draw(Option_t *option)
{
/// Draw the sector on the current pad
/// The first letter of \a option is treated as follows:
/// - case "S" : each row segments are drawn separately
/// - case ""  : the whole subzone is drawn at once
/// in both cases, the rest of the option is passed
/// as argument to the Draw function of respectively
/// zone or row objects.

  if (!fSubZone) return;
  AliMpGraphContext *gr = AliMpGraphContext::Instance();

  gr->Push();
  InitGraphContext();
  switch (option[0]){
  case 'S':
    {

	for (Int_t iRowSeg=0;iRowSeg<fSubZone->GetNofRowSegments();++iRowSeg){
	  gr->Push();
	  AliMpVRowSegment* rowSegment = fSubZone->GetRowSegment(iRowSeg);

	  gr->SetPadPosForReal(TVector2(rowSegment->GetPositionX(),rowSegment->GetPositionY()),
			       TVector2(rowSegment->GetDimensionX(),rowSegment->GetDimensionY()));
	  gr->SetColor(GetColor());
	  DrawObject(rowSegment,option+1);
      
	  gr->Pop();
	}
    }
    break;
  default: AppendPad(option);
  }
  gr->Pop();
}


//_______________________________________________________________________
void AliMpSubZonePainter::Paint(Option_t *option)
{
/// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fSubZone) return;
  if (fSubZone->GetNofRowSegments()<1) return;
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  Int_t col=gVirtualX->GetFillColor();
  InitGraphContext();
  
  gVirtualX->SetFillColor(GetColor());
  for (Int_t iRowSeg=0;iRowSeg<fSubZone->GetNofRowSegments();++iRowSeg){
    AliMpVRowSegment *rowSegment = fSubZone->GetRowSegment(iRowSeg);
    TVector2 pos,dim;
    gr->RealToPad(TVector2(rowSegment->GetPositionX(),rowSegment->GetPositionY()),
                  TVector2(rowSegment->GetDimensionX(),rowSegment->GetDimensionY()),
		  pos,dim);
    gPad->PaintBox(pos.X()-dim.X(),pos.Y()-dim.Y(),
		   pos.X()+dim.X(),pos.Y()+dim.Y());
    if (option[0]=='T'){
      Float_t textSize =   gVirtualX->GetTextSize();
      gVirtualX->SetTextSize(15);
      gPad->PaintText(pos.X()-0.01,pos.Y()-0.01,
		      fSubZone->GetMotif()->GetID());
      gVirtualX->SetTextSize(textSize);
    }
  }

  gVirtualX->SetFillColor(col);
  gr->Pop();
}
