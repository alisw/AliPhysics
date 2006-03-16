#ifndef ALIMUONPIXEL_H
#define ALIMUONPIXEL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONPixel
/// \brief Pixel for MLEM method of cluster finding

#include <TObject.h>

class AliMUONPixel : public TObject {

 public:

  AliMUONPixel();
  AliMUONPixel(Double_t xc, Double_t yc, Double_t wx, Double_t wy, Double_t charge); // constructor
  virtual ~AliMUONPixel(); // Destructor

  Double_t Charge(void) const { return fCharge; } // pixel charge
  Double_t Size(Int_t ixy) const { return fSize[ixy]; } // pixel size
  Double_t Coord(Int_t ixy) const { return fXY[ixy]; } // pixel coordinate
  Int_t Flag(void) const { return fFlag; } // pixel flag
  
  void SetCharge(Double_t Charge) { fCharge = Charge; } // set charge
  void SetSize(Int_t ixy, Double_t Size) { fSize[ixy] = Size; } // set size
  void SetCoord(Int_t ixy, Double_t Coord) { fXY[ixy] = Coord; } // set coordinate
  void SetFlag(Int_t flag) { fFlag = flag; } // set flag
  void Shift(Int_t ixy, Double_t shift) { fXY[ixy] += shift; }
  void Print(const char* opt="") const;
  // What is necessary for sorting TObjArray's
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* pixel) const; // "Compare" function for sorting

 protected:

 private:
 
  Double_t fCharge; // ! pixel charge
  Double_t fSize[2]; // ! pixel size
  Double_t fXY[2]; // ! pixel coordinates
  Int_t fFlag; // ! pixel flag

  // Functions

  ClassDef(AliMUONPixel,0) // pixel for MLEM method of cluster finding
    };
#endif
