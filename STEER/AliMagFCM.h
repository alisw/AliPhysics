#ifndef ALIMAGFCM_H
#define ALIMAGFCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMagF.h"
class TVector;

class AliMagFCM : public AliMagF
{
  //Alice Magnetic Field with constan mesh

public:
  AliMagFCM(){}
  AliMagFCM(const char *name, const char *title, const Int_t integ,
	   const Float_t factor, const Float_t fmax);
  AliMagFCM(const AliMagFCM &mag);
  virtual ~AliMagFCM() {delete fB;}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField();
  void Copy(AliMagFCM &magf) const;
  virtual AliMagFCM & operator=(const AliMagFCM &magf);

  Float_t Bx(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix));
  }
  Float_t By(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+1);
  }
  Float_t Bz(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+2);
  }

protected:

  Float_t    fXbeg;  // Start of mesh in x
  Float_t    fYbeg;  // Start of mesh in y
  Float_t    fZbeg;  // Start of mesh in z
  Float_t    fXdel;  // Mesh step in x
  Float_t    fYdel;  // Mesh step in y
  Float_t    fZdel;  // Mesh step in z
  Double_t   fXdeli; // Inverse of Mesh step in x
  Double_t   fYdeli; // Inverse of Mesh step in y
  Double_t   fZdeli; // Inverse of Mesh step in z
  Int_t      fXn;    // Number of mesh points in x
  Int_t      fYn;    // Number of mesh points in y
  Int_t      fZn;    // Number of mesh points in z
  TVector   *fB;     // Field map
  
  ClassDef(AliMagFCM,1)  //Class for all Alice MagField with Constant Mesh
};

#endif
