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
// $MpId: AliMpZonePainter.cxx,v 1.8 2006/05/24 13:58:32 ivana Exp $
// Category: graphics

//-----------------------------------------------------------------------------
// Class AliMpZonePainter
// ----------------------
// Class for drawing a zone into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
//-----------------------------------------------------------------------------
 
#include "AliMpZonePainter.h"
#include "AliMpGraphContext.h"
#include "AliMpZone.h"
#include "AliMpSubZone.h"
#include "AliMpVRowSegment.h"

#include <TVirtualX.h>
#include <TPad.h>

/// \cond CLASSIMP
ClassImp(AliMpZonePainter)
/// \endcond

//_______________________________________________________________________
AliMpZonePainter::AliMpZonePainter()
  : AliMpVPainter(),
    fZone(0)
{
  /// Default constructor
}

//_______________________________________________________________________
AliMpZonePainter::AliMpZonePainter(AliMpZone *zone)
  : AliMpVPainter(),
    fZone(zone)
{
  /// Standard constructor 

}

//_______________________________________________________________________
AliMpZonePainter::~AliMpZonePainter()
{
  /// Destructor 
}

//_______________________________________________________________________
Int_t AliMpZonePainter::DistancetoPrimitive(Int_t x, Int_t y)
{
  /// Distance to the nearest segment center if (x,y) is inside the zone
  /// 9999 otherwise

  if (fZone->GetNofSubZones()<1) return 9999;
  AliMpGraphContext *gr = AliMpGraphContext::Instance();

  gr->Push();
  InitGraphContext();
  
  TVector2 point = TVector2(gPad->AbsPixeltoX(x), gPad->AbsPixeltoY(y));

  Double_t res=9999.;
  for (Int_t isub=0;isub<fZone->GetNofSubZones();++isub){
    // for each sub-zones
    AliMpSubZone* sub = fZone->GetSubZone(isub);
    for (Int_t iseg=0;iseg<sub->GetNofRowSegments();++iseg){
      //for each row segments
      AliMpVRowSegment* seg = sub->GetRowSegment(iseg);

      TVector2 pos,dim;
      gr->RealToPad(TVector2(seg->GetPositionX(), seg->GetPositionY()),
                    TVector2(seg->GetDimensionX(),seg->GetDimensionY()),pos,dim);

      if ( IsInside(point,pos,dim) ){
	Double_t value = (point-pos).Mod();
	if (value<res) res=value;
      }
    }
  }
  gr->Pop();
  return (Int_t)res;
}

//_______________________________________________________________________
void AliMpZonePainter::DumpObject()
{
  /// Dump the owned object

  fZone->Dump();
}

//_______________________________________________________________________
TVector2 AliMpZonePainter::GetPosition() const
{
  //// Get the owned object's position

  if (fZone->GetNofSubZones()<1) return TVector2(0.,0.);

  TVector2 bl(9999,9999),ur(-9999,-9999);

  for (Int_t isub=0;isub<fZone->GetNofSubZones();++isub){
    // for each sub-zones
    AliMpSubZone* sub = fZone->GetSubZone(isub);
    for (Int_t iseg=0;iseg<sub->GetNofRowSegments();++iseg){
      //for each row segments
      AliMpVRowSegment* seg = sub->GetRowSegment(iseg);

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
    } //iseg
  } //isub
  return (ur+bl)/2.;
}

//_______________________________________________________________________
TVector2 AliMpZonePainter::GetDimensions() const
{
  //// Get the owned object's dimensions

  if (fZone->GetNofSubZones()<1) return TVector2(0.,0.);

  TVector2 bl(9999,9999),ur(-9999,-9999);

  for (Int_t isub=0;isub<fZone->GetNofSubZones();++isub){
    // for each sub-zones
    AliMpSubZone* sub = fZone->GetSubZone(isub);
    for (Int_t iseg=0;iseg<sub->GetNofRowSegments();++iseg){
      //for each row segments
      AliMpVRowSegment* seg = sub->GetRowSegment(iseg);

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
    } //iseg
  } //isub
  return (ur-bl)/2.;
}

//_______________________________________________________________________
void AliMpZonePainter::Draw(Option_t *option)
{
/// Draw the sector on the current pad
/// The first letter of \a option is treated as follows:
/// - case "S" : each sub zones are drawn separately
/// - case ""  : the whole zone is drawn at once
/// in both cases, the rest of the option is passed
/// as argument to the Draw function of respectively
/// zone or row objects.

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fZone) return;

  gr->Push();
  InitGraphContext();
  switch (option[0]){
  case 'S':
    {
      for (Int_t iSubZone=0;iSubZone<fZone->GetNofSubZones();++iSubZone){
	AliMpSubZone *subZone = fZone->GetSubZone(iSubZone);
	gr->Push();

	Double_t blx=  9999,  bly=  9999;
	Double_t urx= -9999,  ury= -9999;

	  for (Int_t iRowSeg=0;iRowSeg<subZone->GetNofRowSegments();++iRowSeg){
	    AliMpVRowSegment *rowSegment = subZone->GetRowSegment(iRowSeg);

	    TVector2 bl = TVector2(rowSegment->GetPositionX(),rowSegment->GetPositionY());
	    bl -= TVector2(rowSegment->GetDimensionX(),rowSegment->GetDimensionY());
	    TVector2 ur = TVector2(rowSegment->GetPositionX(),rowSegment->GetPositionY());
	    ur += TVector2(rowSegment->GetDimensionX(),rowSegment->GetDimensionY());
	    
	    if (bl.X()<blx) blx=bl.X();
	    if (bl.Y()<bly) bly=bl.Y();
	    if (ur.X()>urx) urx=ur.X();
	    if (ur.Y()>ury) ury=ur.Y();
	  }
	  TVector2 position ( (urx+blx)/2.,(ury+bly)/2. );
	  TVector2 dimensions( (urx-blx)/2.,(ury-bly)/2. );
      
	  gr->SetPadPosForReal(position,dimensions);
	  gr->SetColor((fZone->GetID()-1)*5+iSubZone+2);
	  DrawObject(subZone,option+1);
      
	  gr->Pop();
      }
    }
    break;
  default: AppendPad(option);
  }
  gr->Pop();
}

//_______________________________________________________________________
void AliMpZonePainter::Paint(Option_t *option)
{
//// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fZone) return;
  if (fZone->GetNofSubZones()<1) return;
  gr->Push();
  gPad->Range(0.,0.,1.,1.);

  Int_t col=gVirtualX->GetFillColor();
  InitGraphContext();
  
  gVirtualX->SetFillColor(GetColor());
  Float_t textSize =   gVirtualX->GetTextSize();
  for (Int_t iSubZone=0;iSubZone<fZone->GetNofSubZones();++iSubZone){
    AliMpSubZone *subZone = fZone->GetSubZone(iSubZone);
    for (Int_t iRowSeg=0;iRowSeg<subZone->GetNofRowSegments();++iRowSeg){
      AliMpVRowSegment *rowSegment = subZone->GetRowSegment(iRowSeg);
      TVector2 pos,dim;
      gr->RealToPad(TVector2(rowSegment->GetPositionX(),rowSegment->GetPositionY()),
                    TVector2(rowSegment->GetDimensionX(),rowSegment->GetDimensionY()),
		    pos,dim);
      gPad->PaintBox(pos.X()-dim.X(),pos.Y()-dim.Y(),
		     pos.X()+dim.X(),pos.Y()+dim.Y());
      if (option[0]=='T'){
	gVirtualX->SetTextSize(15);
	gPad->PaintText(pos.X()-0.01,pos.Y()-0.01,
			Form("%d",fZone->GetID()));
      }
    }
  }
  gVirtualX->SetTextSize(textSize);
  gVirtualX->SetFillColor(col);
  gr->Pop();
}
