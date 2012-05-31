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
// $MpId: AliMpSectorPainter.cxx,v 1.8 2006/05/24 13:58:32 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMpSectorPainter
// ------------------------
// Class for drawing a sector into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
//-----------------------------------------------------------------------------
  
#include "AliMpSectorPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpSector.h"
#include "AliMpZone.h"
#include "AliMpSubZone.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"

#include <TVirtualX.h>
#include <TPad.h>

/// \cond CLASSIMP
ClassImp(AliMpSectorPainter)
/// \endcond

//_______________________________________________________________________
AliMpSectorPainter::AliMpSectorPainter()
  :AliMpVPainter(),
   fSector(0)
{
  /// Default constructor
}
//_______________________________________________________________________
AliMpSectorPainter::AliMpSectorPainter(AliMpSector *sector)
  :AliMpVPainter(),
   fSector(sector)
{
  /// Standard constructor 

}

//_______________________________________________________________________
AliMpSectorPainter::~AliMpSectorPainter()
{
  /// Destructor
}

//_______________________________________________________________________
void AliMpSectorPainter::DumpObject()
{
/// Draw the owned object

  fSector->Dump();
}

//_______________________________________________________________________
TVector2 AliMpSectorPainter::GetPosition() const
{
/// Get the owned object's position

  if (fSector->GetNofRows()<1) return TVector2(0.,0.);
  AliMpRow* row = fSector->GetRow(0);

  // bl = bottom left position;
  TVector2 bl = TVector2(row->GetPositionX(), row->GetPositionY())-
                TVector2(row->GetDimensionX(), row->GetDimensionY());
  // ur = upper right position
  TVector2 ur = TVector2(row->GetPositionX(), row->GetPositionY())+
                TVector2(row->GetDimensionX(), row->GetDimensionY());;

  for (Int_t irow=1;irow<fSector->GetNofRows();++irow){
    row = fSector->GetRow(irow);
    // update the bottom-left corner
    if (bl.X()>row->GetPositionX()-row->GetDimensionX())
      bl.Set(row->GetPositionX()-row->GetPositionX(),bl.Y());
    if (bl.Y()>row->GetPositionY()-row->GetDimensionY())
      bl.Set(bl.X(),row->GetPositionY()-row->GetDimensionY());
    // update the upper-right corner
    if (ur.X()<row->GetPositionX()+row->GetDimensionX())
      ur.Set(row->GetPositionX()+row->GetDimensionX(),ur.Y());
    if (ur.Y()<row->GetPositionY()+row->GetDimensionY())
      ur.Set(ur.X(),row->GetPositionY()+row->GetDimensionY());
  }
  return (ur+bl)/2.;
}

//_______________________________________________________________________
TVector2 AliMpSectorPainter::GetDimensions() const
{
/// Get the owned object's dimensions

  if (fSector->GetNofRows()<1) return TVector2(0.,0.);
  AliMpRow* row = fSector->GetRow(0);


  // bl = bottom left position
  TVector2 bl = TVector2(row->GetPositionX(), row->GetPositionY())-
                TVector2(row->GetDimensionX(), row->GetDimensionY());
  // ur = upper right position
  TVector2 ur = TVector2(row->GetPositionX(), row->GetPositionY())+
                TVector2(row->GetDimensionX(), row->GetDimensionY());

  for (Int_t irow=1;irow<fSector->GetNofRows();++irow){
    row = fSector->GetRow(irow);
    // update the bottom-left corner
    if (bl.X()>row->GetPositionX()-row->GetDimensionX())
      bl.Set(row->GetPositionX()-row->GetDimensionX(),bl.Y());
    if (bl.Y()>row->GetPositionY()-row->GetDimensionY())
      bl.Set(bl.X(),row->GetPositionY()-row->GetDimensionY());
    // update the upper-right corner
    if (ur.X()<row->GetPositionX()+row->GetDimensionX())
      ur.Set(row->GetPositionX()+row->GetDimensionX(),ur.Y());
    if (ur.Y()<row->GetPositionY()+row->GetDimensionY())
      ur.Set(ur.X(),row->GetPositionY()+row->GetDimensionY());
  }
  return (ur-bl)/2.;

}
//_______________________________________________________________________
void AliMpSectorPainter::Draw(Option_t *option)
{
/// Draw the sector on the current pad
/// The first letter of \a option is treated as follows:
/// - case "Z" : each zones are drawn separately
/// - case "R" : each rows are drawn separately
/// - case ""  : the whole sector is drawn at once
/// in both cases, the rest of the option is passed
/// as argument to the Draw function of respectively
/// zone or row objects.

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fSector) return;
  gr->Push();
  InitGraphContext();

  switch (option[0]){
  case 'Z':
    {
      for (Int_t iZone=1;iZone<=fSector->GetNofZones();++iZone){
	AliMpZone *zone = fSector->GetZone(iZone);
	gr->Push();

	Double_t blx=  9999,  bly=  9999;
	Double_t urx= -9999,  ury= -9999;
	
	for (Int_t iSubZone=0;iSubZone<zone->GetNofSubZones();++iSubZone){
	  AliMpSubZone *subZone = zone->GetSubZone(iSubZone);
	  for (Int_t iRowSeg=0;iRowSeg<subZone->GetNofRowSegments();++iRowSeg){
	    AliMpVRowSegment *rowSegment = subZone->GetRowSegment(iRowSeg);
	    
	    TVector2 bl = TVector2(rowSegment->GetPositionX(),
                                   rowSegment->GetPositionY()) -
                          TVector2(rowSegment->GetDimensionX(),
                                   rowSegment->GetDimensionY());
	    TVector2 ur = TVector2(rowSegment->GetPositionX(),
                                   rowSegment->GetPositionY())+
                          TVector2(rowSegment->GetDimensionX(),
                                   rowSegment->GetDimensionY());
	    
	    if (bl.X()<blx) blx=bl.X();
	    if (bl.Y()<bly) bly=bl.Y();
	    if (ur.X()>urx) urx=ur.X();
	    if (ur.Y()>ury) ury=ur.Y();
	  }
	}
	TVector2 position ( (urx+blx)/2.,(ury+bly)/2. );
	TVector2 dimensions( (urx-blx)/2.,(ury-bly)/2. );

	gr->SetPadPosForReal(position,dimensions);
	gr->SetColor(iZone+3);
	DrawObject(zone,option+1);

	gr->Pop();
      }
    }
    break;
  case 'R':
    {
      for (Int_t iRow=0;iRow<fSector->GetNofRows();++iRow){
	AliMpRow *row = fSector->GetRow(iRow);
	gr->Push();
	gr->SetPadPosForReal(TVector2(row->GetPositionX(), row->GetPositionY()),
                             TVector2(row->GetDimensionX(), row->GetDimensionY()));
	DrawObject(row,option+1);
	gr->Pop();
      }
    }
    break;
  default: AppendPad(option);
  }
  gr->Pop();
}

//_______________________________________________________________________
void AliMpSectorPainter::Paint(Option_t* /*option*/)
{
/// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fSector) return;
  if (fSector->GetNofRows()<1) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  InitGraphContext();
  gPad->Range(0.,0.,1.,1.);


  Double_t lx1=0.;
  Double_t lx2=0.;
  Int_t iRow;
  for (iRow=0;iRow<fSector->GetNofRows();++iRow){
    AliMpRow *row = fSector->GetRow(iRow);
    TVector2 pos,dim;
    gr->RealToPad(TVector2(row->GetPositionX(), row->GetPositionY()),
                  TVector2(row->GetDimensionX(), row->GetDimensionY()),pos,dim);
    gPad->PaintBox(pos.X()-dim.X(),pos.Y()-dim.Y(),
		   pos.X()+dim.X(),pos.Y()+dim.Y());
    gPad->PaintLine(pos.X()-dim.X(),pos.Y()-dim.Y(),
		   pos.X()-dim.X(),pos.Y()+dim.Y());
    gPad->PaintLine(pos.X()+dim.X(),pos.Y()-dim.Y(),
		   pos.X()+dim.X(),pos.Y()+dim.Y());
	  
    if (iRow>0){
      gPad->PaintLine(pos.X()-dim.X(),pos.Y()-dim.Y(),lx1,pos.Y()-dim.Y());
      gPad->PaintLine(pos.X()+dim.X(),pos.Y()-dim.Y(),lx2,pos.Y()-dim.Y());
    }
    lx1=pos.X()-dim.X();
    lx2=pos.X()+dim.X();
  }

  // now we draw the lower and upper horizontal lines

  AliMpRow *row = fSector->GetRow(0);
  TVector2 pos,dim;
  gr->RealToPad(TVector2(row->GetPositionX(), row->GetPositionY()),
                TVector2(row->GetDimensionX(), row->GetDimensionY()),pos,dim);
  gPad->PaintLine(pos.X()-dim.X(),pos.Y()-dim.Y(),
		 pos.X()+dim.X(),pos.Y()-dim.Y());
  
  row = fSector->GetRow(fSector->GetNofRows()-1);
  gr->RealToPad(TVector2(row->GetPositionX(), row->GetPositionY()),
                TVector2(row->GetDimensionX(), row->GetDimensionY()),pos,dim);
  gPad->PaintLine(pos.X()-dim.X(),pos.Y()+dim.Y(),
		 pos.X()+dim.X(),pos.Y()+dim.Y());

  gr->Pop();
  gVirtualX->SetFillColor(col);
}

