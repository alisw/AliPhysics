#ifndef ALIMUONPAINTERCONTOUR_H
#define ALIMUONPAINTERCONTOUR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterContour
/// \brief Contour(s) of a painter
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class AliMpArea;
class TObjArray;
class TPolyLine;
class TGeoHMatrix;
class TVector2;

class AliMUONPainterContour : public TNamed
{
public:
  AliMUONPainterContour(const char* name="");
  AliMUONPainterContour(const char* name, const AliMpArea& area);
  AliMUONPainterContour(const AliMUONPainterContour& rhs);
  AliMUONPainterContour& operator=(const AliMUONPainterContour& rhs);
  virtual ~AliMUONPainterContour();
  
  AliMpArea Area() const;
  
  /// Add an offset to all points
  void Offset(const TVector2& offset);
  
  /// Apply a global transformation to all points
  void Transform(const TGeoHMatrix& matrix);
  
  void AdoptPolyLine(TPolyLine* line);

  virtual void Copy(TObject& obj) const;
    
  Bool_t IsInside(Double_t x, Double_t y) const;

  /// Paint the outline
  void Paint(Option_t* ="") { PaintOutline(1,1); }
  
  void PaintOutline(Int_t lineColor, Int_t lineWidth);
  
  void PaintArea(Int_t fillColor, Int_t fillStyle=1001);

  virtual void Print(Option_t* opt="") const;
  
  /// Return as an array of polylines
  const TObjArray* AsPolyLines() const { return fPolyLines; }
  
private:
  TObjArray* fPolyLines; ///< the polylines used to represent to contour
  Double_t fXmin; ///< min x-value
  Double_t fXmax; ///< max x-value
  Double_t fYmin; ///< min y-value
  Double_t fYmax; ///< max y-value
  
  ClassDef(AliMUONPainterContour,1) // Contour for one painter
};

#endif
