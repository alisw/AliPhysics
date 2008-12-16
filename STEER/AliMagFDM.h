#ifndef ALIMAGFDM_H
#define ALIMAGFDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//   Field with Magnetic Field map
//   Used by AliRun class
//   Author:
//-------------------------------------------------------------------------

#include "AliMagFC.h"
//
class AliMagFDM : public AliMagFC
{
//Alice Magnetic Field:Magnetic field map from IP to muon filter for Muon arm

public:
  AliMagFDM();
  AliMagFDM(const char *name, const char *title, Int_t integ,
	    Float_t factor, Float_t fmax);
  virtual ~AliMagFDM(){} 
  virtual void Field(float *x, float *b) const;
  virtual void Field(double *x, double *b) const;
  virtual void ReadField(); 
  virtual void SetSolenoidField(Float_t field = 2.) {fSolenoid = field;}
  virtual Float_t SolenoidField() const {
     return -Factor()*fSolenoid;
  }
  Int_t FZ(Double_t u, const Float_t *Ar, Float_t du, Int_t ki, Int_t nu) const;
  Double_t Ba(Int_t kai, Double_t za1, Double_t za2, Double_t al1, Double_t al2, Double_t al3, Int_t ka, Int_t ma) const;
  Double_t Bb(Double_t z1, Double_t z2, Double_t y1, Double_t y2, Double_t x1, Double_t x2, Int_t kvr, Int_t k, Int_t l, Int_t m) const; 


protected:

//
  Float_t fSolenoid; // Solenoid Field Strength
  Int_t   fInd;      // Character number of validity Map region

  Float_t fZmin;  // Start of the cartesian  part  of MAP in z
  Float_t fZmax;  // End of Map in z   
  Float_t fYmax;  // Start of the cartesian  part  of MAP in y
  Float_t fYmin;  // End  of the cartesian  part  of MAP in y
  Float_t fZpmx;  // End of the polar  part  of MAP in z
  Float_t fZpmn;  // Start of the polar  part  of MAP in z
  Float_t fRmax;  // Maximal radius of the polar  part  of MAP 
  Float_t fRmin;  // Minimal radius of the polar  part  of MAP  
              

  Float_t    fXdel;  //  step in x - cartesian  part  of MAP
  Float_t    fYdel;  //  step in y - cartesian  part  of MAP
  Float_t    fZdel;  //  step in z - cartesian  part  of MAP
  
  Float_t    fRdel;  //  step in r - polar  part  of MAP
  Float_t    fPhid;  //  step in Phi - polar  part  of MAP
  Float_t    fZpdl;  //  step in z - polar  part  of MAP 
  
  Float_t    fCx1;   // Field constant
  Float_t    fCx2;   // Field constant
  Float_t    fAx1;   // Field constant
  Float_t    fAx2;   // Field constant
   
  Float_t fZc[81];  // z coordinates in cartesian  part
  Float_t fY[81];   // y coordinates in cartesian  part 
  Float_t fBcx[81][81][44]; // Bx array for cartesian  part
  Float_t fBcy[81][81][44]; // By array for cartesian  part
  Float_t fBcz[81][81][44]; // Bz array for cartesian  part

  Float_t  fZp[51];  // z coordinates in polar  part
  Float_t  fR[10];   // r coordinates in polar  part  
  Float_t  fPhi[33]; // Phi coordinates in polar  part

  Float_t  fBpx[51][10][33]; // Bx array for polar  part
  Float_t  fBpy[51][10][33]; // By array for polar  part
  Float_t  fBpz[51][10][33]; // Bx array for polar  part 
  Float_t  fB[2][2][32];     // Limits of field
  
  Int_t      fXl;    // Number steps in x for cartesian  part
  Int_t      fYl;    // Number steps in y  for cartesian  par
  Int_t      fZl;    // Number steps in z  for cartesian  part
    
  Int_t      fRn;    // Number steps in r for polar  part
  Int_t      fPhin;  // Number steps in Phi for polar  part
  Int_t      fZpl;   // Number steps in z for polar  part 
  

  ClassDef(AliMagFDM,1) //Class Magnetic field map from IP till muon filter
};

#endif
