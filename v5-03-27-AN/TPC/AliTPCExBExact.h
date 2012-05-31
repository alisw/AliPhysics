/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
// An exact implementation of the ExB effect.
////

#ifndef ALITPCEXBEXACT_H
#define ALITPCEXBEXACT_H

#include "AliTPCExB.h"

class AliMagF;

class AliTPCExBExact:public AliTPCExB  {
public:
  AliTPCExBExact(); // just for the I/O stuff
  //AliTPCExBExact(const AliMagF *bFieldMap,Double_t driftVelocity,Int_t n=100);
  AliTPCExBExact(const AliMagF *bField,Double_t driftVelocity,Int_t n=100,
		 Int_t nx=30,Int_t ny=30,Int_t nz=100);
  virtual ~AliTPCExBExact();
  virtual void Correct(const Double_t *position,Double_t *corrected);
  //void TestThisBeautifulObject(const AliFieldMap *bFieldMap,const char* fileName);
  void TestThisBeautifulObject(const AliMagF *bField,const char* fileName);
protected:
  Double_t fDriftVelocity; // The electron drift velocity.
private:
  AliTPCExBExact& operator=(const AliTPCExBExact&); // don't assign me
  AliTPCExBExact(const AliTPCExBExact&); // don't copy me
  void TestThisBeautifulObjectGeneric(const char* fileName);
  void CreateLookupTable();
  void GetE(Double_t *e,const Double_t *x) const;
  void GetB(Double_t *b,const Double_t *x) const;
  void Motion(const Double_t *x,Double_t t,Double_t *dxdt) const;
  void CalculateDistortion(const Double_t *x,Double_t *dist) const;
  void DGLStep(Double_t *x,Double_t t,Double_t h) const;
  //const AliFieldMap *fkMap; //! the magnetic field map as supplied by the user
  const AliMagF *fkField;   //! the magnetic field as supplied by the user
  Int_t fkN;        // max number of integration steps
  Int_t fkNX;       // field mesh points in x direction
  Int_t fkNY;       // field mesh points in y direction
  Int_t fkNZ;       // field mesh points in z direction
  Double_t fkXMin;  // the first grid point in x direction
  Double_t fkXMax;  // the last grid point in x direction
  Double_t fkYMin;  // the first grid point in y direction
  Double_t fkYMax;  // the last grid point in y direction
  Double_t fkZMin;  // the first grid point in z direction
  Double_t fkZMax;  // the last grid point in z direction
  Int_t fkNLook;    // size of the lookup table
  Double_t *fkLook; //[fkNLook] the great lookup table
  static const Double_t fgkEM; //! elementary charge over electron mass (C/kg)
  static const Double_t fgkDriftField; //! the TPC drift field (V/m) (modulus)

  ClassDef(AliTPCExBExact,1)
};

#endif
