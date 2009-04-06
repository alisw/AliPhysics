#ifndef ALIMUONPOLYGON_H
#define ALIMUONPOLYGON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONPolygon
/// \brief A planar polygon
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONPolygon : public TObject
{
public:
  AliMUONPolygon(Int_t nvertices=5);
  AliMUONPolygon(Double_t xpos, Double_t ypos, Double_t halfsizex, Double_t halfsizey);
  AliMUONPolygon(const AliMUONPolygon& rhs);
  AliMUONPolygon& operator=(const AliMUONPolygon& rhs);
  virtual ~AliMUONPolygon();

  virtual TObject* Clone(const char* /*newname*/="") const { return new AliMUONPolygon(*this); }
  
  Bool_t Contains(Double_t x, Double_t y) const;
  
  Double_t SignedArea() const;
  
  /// Whether this polygon is oriented counter clockwise
  Bool_t IsCounterClockwiseOriented() const { return SignedArea() > 0.0; }
  
  void ReverseOrientation();
  
  void SetVertex(Int_t i, Double_t x, Double_t y);

  Double_t X(Int_t i) const { return fX[i]; }

  Double_t Y(Int_t i) const { return fY[i]; }

  Int_t NumberOfVertices() const { return fN; }
  
  void Print(Option_t* opt="") const;
  
  void Copy(TObject& obj) const;
  
  void Close();
  
private:
  Int_t fN; /// Number of vertices 
  Double_t* fX; //[fN]
  Double_t* fY; //[fN]
  
  ClassDef(AliMUONPolygon,1) // A simple polygon
};

#endif
