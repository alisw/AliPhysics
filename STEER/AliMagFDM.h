#ifndef ALIMAGFDM_H
#define ALIMAGFDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMagF.h"
//
class AliMagFDM : public AliMagF
{
//Alice Magnetic Field:Magnetic field map from IP to muon filter for Muon arm

public:
  AliMagFDM(){}
  AliMagFDM(const char *name, const char *title, const Int_t integ, const Int_t
  map, const Float_t factor, const Float_t fmax);
  virtual ~AliMagFDM(){} 
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField(); 
  

  void FZ(Double_t *u, Float_t *Ar, Float_t *du, Int_t *ki, Int_t *kf, Double_t *a1, Double_t *a2 , Int_t *nu);
  void FRfuncBi(Int_t *kai, Double_t *za1, Double_t *za2, Double_t *al1, Double_t *al2, Double_t *al3, Int_t *ka, Int_t *ma,Double_t  *ba);
  void FGfuncBi(Double_t *z1, Double_t *z2, Double_t *y1, Double_t *y2, Double_t *x1, Double_t *x2, Int_t *kvr, Int_t *k, Int_t *l, Int_t *m, Double_t *bb); 


protected:

//

  Int_t      fdInd;   // Character number of validity Map region

  Float_t fdZmin;  // Start of the cartesian  part  of MAP in z
  Float_t fdZmax;  // End of Map in z   
  Float_t fdYmax;  // Start of the cartesian  part  of MAP in y
  Float_t fdYmin;  // End  of the cartesian  part  of MAP in y
  Float_t fdZpmx;  // End of the polar  part  of MAP in z
  Float_t fdZpmn;  // Start of the polar  part  of MAP in z
  Float_t fdRmax;  // Maximal radius of the polar  part  of MAP 
  Float_t fdRmin;  // Minimal radius of the polar  part  of MAP  
              

  Float_t    fdXdel;  //  step in x - cartesian  part  of MAP
  Float_t    fdYdel;  //  step in y - cartesian  part  of MAP
  Float_t    fdZdel;  //  step in z - cartesian  part  of MAP
  
  Float_t    fdRdel;  //  step in r - polar  part  of MAP
  Float_t    fdPhid;  //  step in Phi - polar  part  of MAP
  Float_t    fdZpdl;  //  step in z - polar  part  of MAP 
  
  Float_t    fdCx1;   // Field constant
  Float_t    fdCx2;   // Field constant
  Float_t    fdAx1;   // Field constant
  Float_t    fdAx2;   // Field constant
   
  Float_t fdZc[81];  // z coordinates in cartesian  part
  Float_t fdY[81];   // y coordinates in cartesian  part 
  Float_t fdBcx[81][81][44]; // Bx array for cartesian  part
  Float_t fdBcy[81][81][44]; // By array for cartesian  part
  Float_t fdBcz[81][81][44]; // Bz array for cartesian  part

  Float_t  fdZp[51];  // z coordinates in polar  part
  Float_t  fdR[10];   // r coordinates in polar  part  
  Float_t  fdPhi[33]; // Phi coordinates in polar  part

  Float_t  fdBpx[51][10][33]; // Bx array for polar  part
  Float_t  fdBpy[51][10][33]; // By array for polar  part
  Float_t  fdBpz[51][10][33]; // Bx array for polar  part 
  Float_t  fdB[2][2][32];     // Limits of field
  
  Int_t      fdXl;    // Number steps in x for cartesian  part
  Int_t      fdYl;    // Number steps in y  for cartesian  par
  Int_t      fdZl;    // Number steps in z  for cartesian  part
    
  Int_t      fdRn;    // Number steps in r for polar  part
  Int_t      fdPhin;  // Number steps in Phi for polar  part
  Int_t      fdZpl;   // Number steps in z for polar  part 
  

  ClassDef(AliMagFDM,1) //Class Magnetic field map from IP till muon filter
};

#endif
