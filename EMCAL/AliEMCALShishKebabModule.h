#ifndef ALIEMCALSHISHKEBABMODULE_H
#define ALIEMCALSHISHKEBABMODULE_H

/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "TNamed.h"
#include "TMath.h"
#include "TVector2.h"

class AliEMCALGeometry;

class AliEMCALShishKebabModule : public TNamed {
 public:
  AliEMCALShishKebabModule(const double theta=TMath::Pi()/2.);
  AliEMCALShishKebabModule(AliEMCALShishKebabModule &leftNeighbor);
  void Init(const double A,const double B);

  virtual ~AliEMCALShishKebabModule(void) {}
  Bool_t GetParameters();
  void DefineName(const double theta);
  void DefineFirstModule();
  void DefineSecondModuleFirstAssumption(); // need for testing

  Double_t Solve(Double_t (*fcn)(Double_t*, Double_t*), Double_t xmin=0., Double_t xmax=1.,
  Int_t npar=0, Double_t *par=0, Double_t eps=1.0e-8, Int_t maxIter=1000);

  static Double_t Y2(double *x, double *par);
  static Double_t YALL(double *x, double *par);

  Double_t GetTheta() const {return fTheta;}
  Double_t GetThetaInDegree() const {return fTheta*180./TMath::Pi();}

  Double_t GetPosX() {return fOK.Y();}
  Double_t GetPosZ() {return fOK.X();}
  Double_t GetPosXfromR() {return fOK.Y() - fgr;}
  Double_t GetA() {return fA;}
  Double_t GetB() {return fB;}

  // geometry info
  static AliEMCALGeometry *fgGeometry; //!
  static Double_t fga; // default 11.2cm 
  static Double_t fgb; // 
  // radius to IP
  static Double_t fgr;

  TVector2 fOK; // position the module center x->y; z->x;
  Double_t fA;  // parameters of line = y = A*z + B
  Double_t fB;  // 
  // service methods
  void Print(const int pri=1) const;  // *MENU*
 protected:
  // size of SK module
  Double_t fTheta; // theta for SK module

  ClassDef(AliEMCALShishKebabModule,0) // Turned Shish-Kebab module 
};

#endif
/* To do
 1. Insert position the center of towers - 2 additional TVector2
 */
