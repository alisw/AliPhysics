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

#include "AliMUONPainterHighlighter.h"

#include "AliMUONVPainter.h"
#include "AliLog.h"

/// \class AliMUONPainterHighlighter
///
/// A special painter which highlights another one.
/// Highlighting is currently a bold yellow outline of the relevant painter
///
/// \author Laurent Aphecetche, Subatech
///

///\cond CLASSIMP
ClassImp(AliMUONPainterHighlighter)
///\endcond

//_____________________________________________________________________________
AliMUONPainterHighlighter::AliMUONPainterHighlighter()
: TObject(), fPainter(0x0), fX(FLT_MAX), fY(FLT_MAX)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPainterHighlighter::~AliMUONPainterHighlighter()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONPainterHighlighter::SetPainter(AliMUONVPainter* painter, Double_t x, Double_t y)
{
  /// Set the painte we should highlight
  
  fPainter = painter;
  fX = x;
  fY = y;
}

//_____________________________________________________________________________
void 
AliMUONPainterHighlighter::Paint(Option_t*)
{
  /// Actually highlight our painter, if we have one
  if ( fPainter ) 
  {
    fPainter->PaintOutline(5,5,fX,fY);
  }
}
