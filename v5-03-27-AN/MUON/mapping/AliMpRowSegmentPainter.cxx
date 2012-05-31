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
// $MpId: AliMpRowSegmentPainter.cxx,v 1.8 2006/05/24 13:58:32 ivana Exp $
// Category: graphics

//-----------------------------------------------------------------------------
// Class AliMpRowSegmentPainter
// ----------------------------
// Class for drawing a motif into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay
//-----------------------------------------------------------------------------
 
#include "AliMpRowSegmentPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpVRowSegment.h"
#include "AliMpRow.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"

#include <TVirtualX.h>
#include <TPad.h>

/// \cond CLASSIMP
ClassImp(AliMpRowSegmentPainter)
/// \endcond

//_______________________________________________________________________
AliMpRowSegmentPainter::AliMpRowSegmentPainter()
  : AliMpVPainter(),
    fRowSegment(0)
{
  /// Default constructor
}

//_______________________________________________________________________
AliMpRowSegmentPainter::AliMpRowSegmentPainter(AliMpVRowSegment *row)
  : AliMpVPainter(),
    fRowSegment(row)
{
  /// Standard constructor 

}

//_______________________________________________________________________
AliMpRowSegmentPainter::~AliMpRowSegmentPainter()
{
  /// Destructor
}

//_______________________________________________________________________
TVector2 AliMpRowSegmentPainter::GetPosition() const
{
/// Get the owned object's position

  return TVector2(fRowSegment->GetPositionX(), fRowSegment->GetPositionY());
}

//_______________________________________________________________________
TVector2 AliMpRowSegmentPainter::GetDimensions() const
{
/// Get the owned object's dimensions

  return TVector2(fRowSegment->GetDimensionX(),fRowSegment->GetDimensionY());
}

//_______________________________________________________________________
void AliMpRowSegmentPainter::DumpObject()
{
/// Draw the owned object

  fRowSegment->Dump();
}

//_______________________________________________________________________
void AliMpRowSegmentPainter::Draw(Option_t *option)
{
/// Draw the sector on the current pad
/// The first letter of \a option is treated as follows:
/// - case "S" : each row segments are drawn separately
/// - case ""  : the whole row is drawn at once
/// in both cases, the rest of the option is passed
/// as argument to the Draw function of respectively
/// zone or row objects.

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if( !fRowSegment) return;

  gr->Push();
  InitGraphContext();
  switch (option[0]){
  case 'M':
    {

      for (Int_t iMotif=0;iMotif<fRowSegment->GetNofMotifs();++iMotif){
	gr->Push();

	Int_t motifPositionId = fRowSegment->GetMotifPositionId(iMotif);
	AliMpMotifPosition *motifPos = 
	    fRowSegment->GetRow()->GetMotifMap()
	      ->FindMotifPosition(motifPositionId);
				     
	gr->SetPadPosForReal(TVector2(motifPos->GetPositionX(), motifPos->GetPositionY()),
                             TVector2(motifPos->GetDimensionX(),motifPos->GetDimensionX()));
      	gr->SetColor(GetColor());
	DrawObject(motifPos,option+1);

	gr->Pop();
      }
    }
    break;
  default: AppendPad(option);
  }
  gr->Pop();
}


//_______________________________________________________________________
void AliMpRowSegmentPainter::Paint(Option_t* /*option*/)
{
/// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fRowSegment) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();
  PaintWholeBox(kTRUE);

//   Float_t textSize =   gVirtualX->GetTextSize();
//   if (option[0]=='T')
//     gPad->PaintText(GetPadPosition().X()-0.01,GetPadPosition().Y()-0.01,
// 		     Form("%d",fRowSegment->GetMotif()->GetID()));
//   gVirtualX->SetTextSize(textSize);
  gr->Pop();
  gVirtualX->SetFillColor(col);
}
