#ifndef ALIMAGFCM_H
#define ALIMAGFCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------
//  Class for Alice magnetic field with constant mesh
//  Used in the configuration macros (macros/Config.C, etc.)
//  Author:
//-----------------------------------------------------------------------

#include "AliMagFC.h"
#include <TVector.h>

class AliMagFCM : public AliMagFC
{
  //Alice Magnetic Field with constant mesh

public:
  AliMagFCM();
  AliMagFCM(const char *name, const char *title, Int_t integ,
	   Float_t factor, Float_t fmax);
  AliMagFCM(const AliMagFCM &mag);
  virtual ~AliMagFCM() {delete fB;}
  virtual void Field(float *x, float *b) const;
  virtual void Field(double *x, double *b) const;
  virtual void ReadField();
  virtual void    SetSolenoidField(Float_t field = 2.) {fSolenoid = field;}
  virtual Float_t SolenoidField() const {
    return -Factor()*fSolenoid;
  }
  
  void Copy(TObject &magf) const;
  virtual AliMagFCM & operator=(const AliMagFCM &magf)
    {magf.Copy(*this); return *this;}

  Float_t Bx(Int_t ix, Int_t iy, Int_t iz) const {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix));
  }
  Float_t By(Int_t ix, Int_t iy, Int_t iz) const {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+1);
  }
  Float_t Bz(Int_t ix, Int_t iy, Int_t iz) const {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+2);
  }

protected:

  Float_t    fXbeg;     // Start of mesh in x
  Float_t    fYbeg;     // Start of mesh in y
  Float_t    fZbeg;     // Start of mesh in z
  Float_t    fXdel;     // Mesh step in x
  Float_t    fYdel;     // Mesh step in y
  Float_t    fZdel;     // Mesh step in z
  Float_t    fSolenoid; // Solenoid Field Strength
  Double_t   fXdeli;    // Inverse of Mesh step in x
  Double_t   fYdeli;    // Inverse of Mesh step in y
  Double_t   fZdeli;    // Inverse of Mesh step in z
  Int_t      fXn;       // Number of mesh points in x
  Int_t      fYn;       // Number of mesh points in y
  Int_t      fZn;       // Number of mesh points in z
  TVector   *fB;        // Field map
  
  ClassDef(AliMagFCM,1)  //Class for all Alice MagField with Constant Mesh
};

#endif
