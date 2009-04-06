#ifndef ALIMUONCONTOUR_H
#define ALIMUONCONTOUR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONContour
/// \brief 2D contour
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

#ifndef ALI_MP_AREA_H
#  include "AliMpArea.h"
#endif

class AliMUONPolygon;
class TGeoHMatrix;
class TObjArray;

class AliMUONContour : public TNamed
{
public:
  AliMUONContour(const char* name="");
  AliMUONContour(const char* name, const AliMpArea& area);
  AliMUONContour(const AliMUONContour& rhs);
  AliMUONContour& operator=(const AliMUONContour& rhs);
  virtual ~AliMUONContour();
  
  AliMpArea Area() const;
  
  /// Get a full copy of this object.
  virtual TObject* Clone(const char* /*newname*/="") const { return new AliMUONContour(*this); }
  
  /// Add an offset to all points
  void Offset(Double_t x, Double_t y);
  
  /// Apply a global transformation to all points
  void Transform(const TGeoHMatrix& matrix);
  
  void Add(const AliMUONPolygon& polygon);
  
  virtual void Copy(TObject& obj) const;
    
  Bool_t IsInside(Double_t x, Double_t y) const;

  virtual void Print(Option_t* opt="") const;
  
  /// Get the number of vertices of this contour
  Int_t NumberOfVertices() const { return fNofVertices; }
  
  Bool_t IsValid() const;
  
  /// Get the list of polygons we have
  const TObjArray* Polygons() const { return fPolygons; }
  
  void AssertOrientation(Bool_t autoCorrect=kFALSE);
  
private:
  TObjArray* fPolygons; ///< the polygons that this contour is made of
  Double_t fXmin; ///< min x-value
  Double_t fXmax; ///< max x-value
  Double_t fYmin; ///< min y-value
  Double_t fYmax; ///< max y-value
  Int_t fNofVertices; ///< total number of vertices
  
  ClassDef(AliMUONContour,1) // 2D-contour of an object
};

#endif
