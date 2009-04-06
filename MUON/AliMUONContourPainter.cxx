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

/// \class AliMUONContourPainter
///
/// Class to draw AliMUONContour objects (2D)
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONContourPainter.h"

#include "TVirtualX.h"
#include "AliMUONPolygon.h"
#include "AliMUONContour.h"
#include "TObjArray.h"
#include "TVirtualPad.h"

///\cond CLASSIMP
ClassImp(AliMUONContourPainter)
///\endcond

//_____________________________________________________________________________
AliMUONContourPainter::AliMUONContourPainter()
{
  /// Ctor
}

//_____________________________________________________________________________
AliMUONContourPainter::~AliMUONContourPainter()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONContourPainter::Paint(const AliMUONContour& contour, 
                             Int_t lineColor, Int_t lineWidth, 
                             Int_t fillColor, Int_t fillStyle)
{
  /// Paint the given contour. 
  /// If lineColor > 0 the outline is drawn
  /// If fillColor > 0 the contour is filled.
  
  Bool_t outline(lineColor>0);
  Bool_t fill(fillColor>0);
  
  Int_t fc = gVirtualX->GetFillColor();
  Int_t fs = gVirtualX->GetFillStyle();
  Int_t lc = gVirtualX->GetLineColor();
  Int_t lw = gVirtualX->GetLineWidth();
  
  if ( lineColor > 0 ) gVirtualX->SetLineColor(lineColor);
  if ( lineWidth > 0 ) gVirtualX->SetLineWidth(lineWidth);
  if ( fillColor > 0 ) gVirtualX->SetFillColor(fillColor);
  if ( fillStyle > 0 ) gVirtualX->SetFillStyle(fillStyle);
  
  TIter next(contour.Polygons());
  AliMUONPolygon* pol;
  while ( ( pol = static_cast<AliMUONPolygon*>(next()) ) )
  {
    Int_t n = pol->NumberOfVertices();
    Double_t* x = new Double_t[n];
    Double_t* y = new Double_t[n];
    for ( Int_t i = 0; i < n; ++i )
    {
      x[i] = gPad->GetLogx() ? gPad->XtoPad(pol->X(i)) : pol->X(i);
      y[i] = gPad->GetLogy() ? gPad->YtoPad(pol->Y(i)) : pol->Y(i);
    }
    if ( fill ) 
    {
      gPad->PaintFillArea(n,x,y);
    }
    if (outline)
    {
      gPad->PaintPolyLine(n,x,y);
    }
    
    delete[] x;
    delete[] y;
  }
  
  gVirtualX->SetFillColor(fc);
  gVirtualX->SetFillStyle(fs);
  gVirtualX->SetLineColor(lc);
  gVirtualX->SetLineWidth(lw);
  
}
