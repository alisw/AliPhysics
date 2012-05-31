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
// $MpId: AliMpRowPainter.cxx,v 1.8 2006/05/24 13:58:32 ivana Exp $
// Category: graphics

//-----------------------------------------------------------------------------
// Class AliMpRowPainter
// ---------------------
// Class for drawing a row into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
//-----------------------------------------------------------------------------
  
#include "AliMpRowPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpRow.h"
#include "AliMpRowSegment.h"

#include <TVirtualX.h>
#include <TPad.h>
 
/// \cond CLASSIMP
ClassImp(AliMpRowPainter)
/// \endcond

//_______________________________________________________________________
AliMpRowPainter::AliMpRowPainter()
  : AliMpVPainter(),
    fRow(0)
{
  /// Default constructor
}

//_______________________________________________________________________
AliMpRowPainter::AliMpRowPainter(AliMpRow *row)
  : AliMpVPainter(),
    fRow(row)
{
  /// Standard constructor 
}

//_______________________________________________________________________
AliMpRowPainter::~AliMpRowPainter()
{
  /// Destructor
}

//_______________________________________________________________________
void AliMpRowPainter::DumpObject()
{
/// Draw the owned object

  fRow->Dump();
}

//_______________________________________________________________________
TVector2 AliMpRowPainter::GetPosition() const
{
/// Get the owned object's position

  return TVector2(fRow->GetPositionX(), fRow->GetPositionY());
}

//_______________________________________________________________________
TVector2 AliMpRowPainter::GetDimensions() const
{
/// Get the owned object's dimensions

  return TVector2(fRow->GetDimensionX(), fRow->GetDimensionY());
}

//_______________________________________________________________________
void AliMpRowPainter::Draw(Option_t *option)
{
/// Draw the sector on the current pad
/// The first letter of \a option is treated as follows:
/// - case "S" : each row segments are drawn separately
/// - case ""  : the whole row is drawn at once
/// in both cases, the rest of the option is passed
/// as argument to the Draw function of respectively
/// zone or row objects.

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if( !fRow) return;

  gr->Push();
  InitGraphContext();

  switch (option[0]){
  case 'S':
    {
      for (Int_t iRowSeg=0;iRowSeg<fRow->GetNofRowSegments();++iRowSeg){
	AliMpVRowSegment *rowSegment = fRow->GetRowSegment(iRowSeg);
	gr->Push();

	gr->SetPadPosForReal(TVector2(rowSegment->GetPositionX(),rowSegment->GetPositionY()),
                             TVector2(rowSegment->GetDimensionX(),rowSegment->GetDimensionY()));
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
void AliMpRowPainter::Paint(Option_t *option)
{
  /// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if( !fRow) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();  
  PaintWholeBox(kTRUE);

  if (option[0]=='T'){
    Float_t textSize =   gVirtualX->GetTextSize();
    gVirtualX->SetTextSize(12);
    gPad->PaintText(GetPadPosition().X()-0.01,GetPadPosition().Y()-0.01,
		    Form("%d",fRow->GetID()));
    gVirtualX->SetTextSize(textSize);
  }


  gr->Pop();
  gVirtualX->SetFillColor(col);
}
