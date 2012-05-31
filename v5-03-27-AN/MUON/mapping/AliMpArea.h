/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpArea.h,v 1.9 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpArea
/// \brief A rectangle area positioned in plane..
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_AREA_H
#define ALI_MP_AREA_H

#include <TObject.h>

class AliMpArea : public TObject
{
 public:
  AliMpArea(Double_t x, Double_t y, 
            Double_t dx, Double_t dy);
  AliMpArea(const AliMpArea& rhs);
  AliMpArea();
  virtual ~AliMpArea();

  // operators
  AliMpArea& operator = (const AliMpArea& right);

  // methods
  Double_t LeftBorder() const;
  Double_t RightBorder() const;
  Double_t UpBorder() const;
  Double_t DownBorder() const;

  void LeftDownCorner(Double_t& x, Double_t& y) const;
  void LeftUpCorner(Double_t& x, Double_t& y) const;
  void RightDownCorner(Double_t& x, Double_t& y) const;
  void RightUpCorner(Double_t& x, Double_t& y) const;

  AliMpArea Intersect(const AliMpArea& area) const;
  Bool_t    Overlap(const AliMpArea& area) const;
  Bool_t    Contains(const AliMpArea& area) const;
  
  void Print(Option_t* opt="") const;

  // get methods
  void      GetParameters(Double_t& x, Double_t& y,
                          Double_t& dx, Double_t& dy) const;
  Double_t  GetPositionX() const;
  Double_t  GetPositionY() const;
  Double_t  GetDimensionX() const;    
  Double_t  GetDimensionY() const;    
  Bool_t    IsValid() const;
  
  
 private:
  // data members
  Double_t  fPositionX;  ///<  x position
  Double_t  fPositionY;  ///<  y position
  Double_t  fDimensionX; ///<   x dimension (half lengths)
  Double_t  fDimensionY; ///<   y dimension (half lengths)
  Bool_t    fValidity;   ///<  validity

  ClassDef(AliMpArea,2) //utility class for area iterators
};

ostream& operator << (ostream &stream,const AliMpArea& area);

// inline functions

                 /// Return x position
inline Double_t  AliMpArea::GetPositionX() const   { return fPositionX; }
                 /// Return y position
inline Double_t  AliMpArea::GetPositionY() const   { return fPositionY; }
                 /// Return x dimensions
inline Double_t  AliMpArea::GetDimensionX() const { return fDimensionX; }    
                 /// Return y dimensions
inline Double_t  AliMpArea::GetDimensionY() const { return fDimensionY; }    
                 /// Return validity
inline Bool_t    AliMpArea::IsValid() const    { return fValidity; }

#endif //ALI_MP_AREA_H
