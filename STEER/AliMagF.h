#ifndef AliMagF_H
#define AliMagF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNamed.h"
#include "TVector.h"

enum Field_t {Undef=1, Const=1, ConMesh=2, DipoMap=3};

class AliMagF : public TNamed {

protected:
  Int_t     fMap;    // Field Map identifier
  Int_t     fType;   // Mag Field type
  Int_t     fInteg;  // Integration method as indicated in Geant
  Float_t   fFactor; // Multiplicative factor
  Float_t   fMax;    // Max Field as indicated in Geant

public:
  AliMagF(){}
  AliMagF(const char *name, const char *title, const Int_t integ, const Int_t map, 
	  const Float_t factor, const Float_t fmax);
  virtual ~AliMagF() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual Int_t Type() {return fType;}
  virtual Float_t Max() const {return fMax;}
  virtual Int_t Map() const {return fMap;}
  virtual Int_t Integ() const {return fInteg;}
  virtual Float_t Factor() const {return fFactor;}
  virtual void ReadField() {}
  
  ClassDef(AliMagF,1)  //Base class for all Alice MagField
};

class AliMagFC  : public AliMagF
{
  //Alice Constant Magnetic Field
private:

public:
  AliMagFC(){}
  AliMagFC(const char *name, const char *title, const Int_t integ, const Int_t map, 
	   const Float_t factor, const Float_t fmax);
  virtual ~AliMagFC() {}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField() {}
  
  ClassDef(AliMagFC,1)  //Class for all Alice Constant MagField 
};

class AliMagFCM : public AliMagF
{
  //Alice Magnetic Field with constan mesh
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
public:
  AliMagFCM(){}
  AliMagFCM(const char *name, const char *title, const Int_t integ, const Int_t map, 
	   const Float_t factor, const Float_t fmax);
  virtual ~AliMagFCM() {delete fB;}
  virtual void Field(Float_t *x, Float_t *b);
  virtual void ReadField();

  inline Float_t Bx(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix));
  }
  inline Float_t By(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+1);
  }
  inline Float_t Bz(const Int_t ix, const Int_t iy, const Int_t iz) {
    return (*fB)(3*(iz*(fXn*fYn)+iy*fXn+ix)+2);
  }
  
  ClassDef(AliMagFCM,1)  //Class for all Alice MagField with Constant Mesh
};
//************************************
//
class AliMagFDM : public AliMagF
{
//Alice Magnetic Field:Magnetic field map from IP to muon filter for Muon arm

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
  
  Float_t    fdCx1, fdCx2;
  Float_t    fdAx1, fdAx2; 
   
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
  Float_t  fdB[2][2][32]; 
  
  Int_t      fdXl;    // Number steps in x for cartesian  part
  Int_t      fdYl;    // Number steps in y  for cartesian  par
  Int_t      fdZl;    // Number steps in z  for cartesian  part
    
  Int_t      fdRn;    // Number steps in r for polar  part
  Int_t      fdPhin;  // Number steps in Phi for polar  part
  Int_t      fdZpl;   // Number steps in z for polar  part 
  
  Float_t  rrtes;

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
//_________________________________________

  ClassDef(AliMagFDM,1) //Class Magnetic field map from IP till muon filter
};


#endif
