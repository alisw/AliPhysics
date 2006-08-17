#ifndef ALIEMCALSHISHKEBABMODULE_H
#define ALIEMCALSHISHKEBABMODULE_H

/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 */

/* $Id$ */

//_________________________________________________________________________
// Main class for "twist" geometry of Shish-Kebab case.
// Author: Aleksei Pavlinov(WSU).
// Sep 2004.

#include <TNamed.h>
#include <TVector2.h>

class AliEMCALGeometry;

class AliEMCALShishKebabModule : public TNamed {
 public:
  AliEMCALShishKebabModule();
  AliEMCALShishKebabModule(AliEMCALShishKebabModule &leftNeighbor);
  void Init(Double_t A, Double_t B);
  AliEMCALShishKebabModule(const AliEMCALShishKebabModule& mod);

  AliEMCALShishKebabModule & operator = (const AliEMCALShishKebabModule& /*rvalue*/)  {
    Fatal("operator =", "not implemented") ;  
    return *this ; 
  }

  virtual ~AliEMCALShishKebabModule(void) {}
  Bool_t GetParameters();
  void DefineName(Double_t theta);
  void DefineFirstModule();
  void DefineSecondModuleFirstAssumption(); // need for testing

  Double_t Solve(Double_t (*fcn)(Double_t*, Double_t*), Double_t xmin=0., Double_t xmax=1.,
  Int_t npar=0, Double_t *par=0, Double_t eps=1.0e-8, Int_t maxIter=1000);

  static Double_t Y2(Double_t *x, Double_t *par);
  static Double_t YALL(Double_t *x, Double_t *par);

  Double_t GetTheta() const {return fTheta;}
  Double_t GetThetaInDegree() const;

  Double_t GetPosX() const {return fOK.Y();}
  Double_t GetPosZ() const {return fOK.X();}
  Double_t GetPosXfromR() const {return fOK.Y() - fgr;}
  Double_t GetA() const {return fA;}
  Double_t GetB() const {return fB;}
  // service methods
  void PrintShish(Int_t pri=1) const;  // *MENU*

 protected:
  // geometry info
  static AliEMCALGeometry *fgGeometry; //!
  static Double_t fga; // x size of module; default 11.2cm 
  static Double_t fgb; // y size of module;
  static Double_t fgr; // radius to IP

  TVector2 fOK; // position the module center x->y; z->x;
  Double_t fA;  // parameters of line; y = A*z + B
  Double_t fB;  // parameters of line; y = A*z + B
  Double_t fTheta; // theta for SK module

  //public:
  ClassDef(AliEMCALShishKebabModule,1) // Turned Shish-Kebab module 
};

#endif
/* To do
 1. Insert position the center of towers - 2 additional TVector2
 */
