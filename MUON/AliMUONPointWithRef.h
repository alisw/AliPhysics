#ifndef ALIMUONPOINTWITHREF_H
#define ALIMUONPOINTWITHREF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONPointWithRef
/// \brief A TVector2 with an integer ref, and a specific Compare
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONPointWithRef : public TObject
{
public:
  AliMUONPointWithRef(Double_t x, Double_t y, Int_t ref);
  AliMUONPointWithRef();
  virtual ~AliMUONPointWithRef() {}
  
  /// We are sortable
  virtual Bool_t IsSortable() const { return kTRUE; }

  virtual Int_t	Compare(const TObject* obj) const;

  Double_t X() const { return fX; }
  
  Double_t Y() const { return fY; }
  
  Int_t Ref() const { return fRef; }
  
  void Print(Option_t* opt="") const;
  
private:
  Double_t fX; //< x value
  Double_t fY; //< y value
  Int_t fRef; //< index of the original point in some array
  
  ClassDef(AliMUONPointWithRef,1) // A point with an external integer reference
};

#endif
