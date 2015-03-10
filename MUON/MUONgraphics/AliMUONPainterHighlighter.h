#ifndef ALIMUONPAINTERHIGHLIGHTER_H
#define ALIMUONPAINTERHIGHLIGHTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterHighlighter
/// \brief Special painter which highlights (i.e. draws a yellow bold outline) another painter
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#include <float.h>

class AliMUONVPainter;

class AliMUONPainterHighlighter : public TObject
{
public:
  AliMUONPainterHighlighter();
  virtual ~AliMUONPainterHighlighter();

  void SetPainter(AliMUONVPainter* painter, Double_t x=FLT_MAX, Double_t y=FLT_MAX);
  
  void Paint(Option_t* opt="");
  
private:
  /// Not implemented
  AliMUONPainterHighlighter(const AliMUONPainterHighlighter& rhs);
  /// Not implemented
  AliMUONPainterHighlighter& operator=(const AliMUONPainterHighlighter& rhs);
  
private:
  AliMUONVPainter* fPainter; //!<! the painter we should highlight
  Double_t fX; //!<! position within painter to be highlighted
  Double_t fY; //!<! position within painter to be highlighted
  
  ClassDef(AliMUONPainterHighlighter,1) // Painter highlighter
};

#endif
