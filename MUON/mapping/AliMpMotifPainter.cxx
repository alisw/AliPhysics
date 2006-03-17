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
// $MpId: AliMpMotifPainter.cxx,v 1.8 2006/03/17 11:35:29 ivana Exp $
// Category: graphics
//
// Class AliMpMotifPainter
// -----------------------
// Class for drawing a motif into canvas
// Included in AliRoot: 2003/05/02
// Authors: David Guez, IPN Orsay

#include "AliMpMotifPainter.h"
#include "AliMpGraphContext.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"

#include <TVirtualX.h>
#include <TPad.h>
 
ClassImp(AliMpMotifPainter)

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter()
  : AliMpVPainter(),
    fMotifPos(0)
{
  /// Default constructor
}

//_______________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter(AliMpMotifPosition *motifPos)
  : AliMpVPainter(),
    fMotifPos(motifPos)
{
  /// Standard constructor 

}

//_____________________________________________________________________________
AliMpMotifPainter::AliMpMotifPainter(const AliMpMotifPainter& right) 
  : AliMpVPainter(right) 
{
  /// Protected copy constructor (not provided) 

  Fatal("AliMpMotifPainter", "Copy constructor not provided.");
}

//_______________________________________________________________________
AliMpMotifPainter::~AliMpMotifPainter()
{
  /// Default constructor
}

//_____________________________________________________________________________
AliMpMotifPainter& 
AliMpMotifPainter::operator=(const AliMpMotifPainter& right)
{
  /// Assignment operator (not provided)

  // check assignment to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignment operator not provided.");
    
  return *this;  
}    

//_______________________________________________________________________
void AliMpMotifPainter::DumpObject()
{
/// Dump the owned object

  fMotifPos->Dump();
}

//_______________________________________________________________________
TVector2 AliMpMotifPainter::GetPosition() const
{
/// Get the owned object's position

  return fMotifPos->Position();
}

//_______________________________________________________________________
TVector2 AliMpMotifPainter::GetDimensions() const
{
/// Get the owned object's dimensions

  return fMotifPos->Dimensions();
}

//_______________________________________________________________________
void AliMpMotifPainter::Paint(Option_t *option)
{
/// Paint the object

  AliMpGraphContext *gr = AliMpGraphContext::Instance();
  if (!fMotifPos) return;
  Int_t col=gVirtualX->GetFillColor();
  gr->Push();
  gPad->Range(0.,0.,1.,1.);
  InitGraphContext();

  switch (option[0]){
  case 'T':
  case 'I':
  case 'X':
    {
      PaintWholeBox();
      Float_t textSize =   gVirtualX->GetTextSize();
      gVirtualX->SetTextSize(10);
      TString str;
      switch (option[0]) {
        case 'T' : 
          str = Form("%d",fMotifPos->GetID());
          break;
        case 'I':{
	  switch (option[1]){
	    case '+' :
            str = Form("(%d,%d)",fMotifPos->GetHighIndicesLimit().GetFirst(),
                               fMotifPos->GetHighIndicesLimit().GetSecond());
      	    break;
	    default:
            str = Form("(%d,%d)",fMotifPos->GetLowIndicesLimit().GetFirst(),
                               fMotifPos->GetLowIndicesLimit().GetSecond());
      	  }
	}
        break;
        case 'X' :
          str = Form("(%f,%f)",(GetPosition()-GetDimensions()).X(),
                               (GetPosition()-GetDimensions()).Y());
          break;
      }
      gPad->PaintText(GetPadPosition().X()-0.01,GetPadPosition().Y()-0.01,str);
      
      gVirtualX->SetTextSize(textSize);
    }
    break;
  case 'P':
    {
      //PaintWholeBox(kFALSE);
      AliMpMotifType *motifType = fMotifPos->GetMotif()->GetMotifType();
      for (Int_t j=motifType->GetNofPadsY()-1;j>=0;j--){
	   for (Int_t i=0;i<motifType->GetNofPadsX();i++){
	     AliMpIntPair indices = AliMpIntPair(i,j);
               AliMpConnection* connect = 
                 motifType->FindConnectionByLocalIndices(indices);
	     if (connect){
    	       TVector2 realPadPos = 
	        GetPosition()+fMotifPos->GetMotif()->PadPositionLocal(indices);
	       TVector2 padPadPos,padPadDim;
	       gr->RealToPad(realPadPos,
                               fMotifPos->GetMotif()->GetPadDimensions(indices),
	     		 padPadPos,padPadDim);
	       TVector2 bl = padPadPos - padPadDim;
	       TVector2 ur = padPadPos + padPadDim;


	       Style_t sty = gVirtualX->GetFillStyle();
	       gVirtualX->SetFillStyle(1);
	       gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
	       gVirtualX->SetFillStyle(0);
	       gPad->PaintBox(bl.X(),bl.Y(),ur.X(),ur.Y());
	       gVirtualX->SetFillStyle(sty);
	       if (option[1]=='T'){
	         Float_t textSize =   gVirtualX->GetTextSize();
	         gVirtualX->SetTextSize(10);
		 gVirtualX->SetTextAlign(22);
		 //	         gPad->PaintText(padPadPos.X()-0.01,padPadPos.Y()-0.01,
	         gPad->PaintText((bl.X()+ur.X())/2.0,(bl.Y()+ur.Y())/2.0,
			      Form("%d",connect->GetGassiNum()));
	      
	         gVirtualX->SetTextSize(textSize);
                 }
	    }
          }
      }
    }
    break;
  default:
    PaintWholeBox(kFALSE);
  }
  gr->Pop();
  gVirtualX->SetFillColor(col);
}
