#ifndef ALIMUONCONTOURPAINTER_H
#define ALIMUONCONTOURPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONContourPainter
/// \brief Class to draw AliMUONContour objects
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONContour;

class AliMUONContourPainter : public TObject
{
public:
  AliMUONContourPainter();
  virtual ~AliMUONContourPainter();
  
  using TObject::Paint;
  
  static void Paint(const AliMUONContour& contour, 
                    Int_t lineColor=1, Int_t lineStyle=1, 
                    Int_t fillColor=-1, Int_t fillStyle=1001);

  ClassDef(AliMUONContourPainter,1) // 
};

#endif
