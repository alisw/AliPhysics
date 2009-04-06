#ifndef ALIMUONSEGMENT_H
#define ALIMUONSEGMENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONSegment
/// \brief A basic line segment, used for contour making algorithm(s)
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONSegment : public TObject
{
public:
  AliMUONSegment();
  AliMUONSegment(Double_t xstart, Double_t ystart, Double_t xend, Double_t yend);
  virtual ~AliMUONSegment() {}
  
  virtual Int_t	Compare(const TObject* obj) const;

  /// We are sortable
  virtual Bool_t IsSortable() const { return kTRUE; }
  
  /// Return the x-coordinate of our starting point
  Double_t StartX() const { return fStartX; }
  /// Return the y-coordinate of our starting point
  Double_t StartY() const { return fStartY; }  
  /// Return the x-coordinate of our ending point
  Double_t EndX() const { return fEndX; }
  /// Return the y-coordinate of our ending point
  Double_t EndY() const { return fEndY; }
    
  /// Return our smallest y (of starting or ending point)
  double SmallerY() const { return fSmallerY; }

  /// Whether we are a horizontal segment
  Bool_t IsHorizontal() const { return fIsHorizontal; }
  
  /// Whethere we are a vertical segment
  Bool_t IsVertical() const { return fIsVertical; }
  
  /// Whether we are a left edge
  Bool_t IsLeftEdge() const { return fIsLeftEdge; }

  /// Whether we are a right edge
  Bool_t IsRightEdge() const { return fIsRightEdge; }
  
  /// Return our bottom y
  double Bottom() const { return SmallerY(); }
  
  double Top() const;
  
  double Distance() const;
  
  /// Whether we're just a point
  Bool_t IsAPoint() const { return fIsAPoint; }
  
  const char* AsString() const;
  
  static Bool_t AreEqual(double a, double b);

  void Print(Option_t* opt="") const;
  
  void Set(Double_t xstart, Double_t ystart, Double_t xend, Double_t yend);
  
private:
  Double_t fStartX; /// x of start point
  Double_t fStartY; /// y of start point
  Double_t fEndX; /// x of end point
  Double_t fEndY; /// y of end point
  Double_t fSmallerY; /// Either StartY or EndY
  Bool_t fIsHorizontal; /// Whether the segment is horizontal
  Bool_t fIsVertical; /// Whether the segment is vertical
  Bool_t fIsLeftEdge; /// Whether the segment is a left edge 
  Bool_t fIsRightEdge; /// Whether the segment is a right edge
  Bool_t fIsAPoint; /// Whether start==end
  
  static const Double_t fgkPrecision; /// Floating point precision used in comparisons
  
  ClassDef(AliMUONSegment,1) // A basic line segment
};


#endif
