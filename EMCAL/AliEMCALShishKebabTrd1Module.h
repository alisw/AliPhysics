#ifndef ALIEMCALSHISHKEBABTRD1MODULE_H
#define ALIEMCALSHISHKEBABTRD1MODULE_H

/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//*-- Author: Aleksei Pavlinov (WSU)
// TO DO : create class for Super module Geometry - 4-nov-04 

#include "TNamed.h"
#include "TMath.h"
#include "TVector2.h"

class AliEMCALGeometry;

class AliEMCALShishKebabTrd1Module : public TNamed {
 public:
  AliEMCALShishKebabTrd1Module(double theta=TMath::Pi()/2.);
  AliEMCALShishKebabTrd1Module(AliEMCALShishKebabTrd1Module &leftNeighbor);
  void Init(double A, double B);

  virtual ~AliEMCALShishKebabTrd1Module(void) {}
  Bool_t GetParameters();
  void DefineName(double theta);
  void DefineFirstModule();

  Double_t GetTheta() const{return fTheta;}
  Double_t GetThetaInDegree() const {return fTheta*180./TMath::Pi();}
  TVector2& GetCenterOfModule() {return fOK;}
  Double_t  GetEtaOfCenterOfModule(){return -TMath::Log(TMath::Tan(fOK.Phi()/2.));}

  Double_t  GetPosX() {return fOK.Y();}
  Double_t  GetPosZ() {return fOK.X();}
  Double_t  GetPosXfromR() {return fOK.Y() - fgr;}
  Double_t  GetA() {return fA;}
  Double_t  GetB() {return fB;}
  //  Additional offline staff 
  TVector2& GetCenterOfCell(Int_t ieta)
  { if(ieta<=1) return fOK1;
    else        return fOK2;}
  // 
  Double_t GetTanBetta() {return fgtanBetta;}
  Double_t Getb()        {return fgb;}
  // service methods
  void Print(int pri=1) const;  // *MENU*

  // geometry info
  static AliEMCALGeometry *fgGeometry; //!
  static Double_t fga;        // 2*dx1=2*dy1
  static Double_t fga2;       // 2*dx2
  static Double_t fgb;        // 2*dz1
  static Double_t fgangle;    // ~1 degree
  static Double_t fgtanBetta; // tan(fgangle/2.)
  // radius to IP
  static Double_t fgr;

 protected:
  TVector2 fOK;     // position the module center x->y; z->x;
  Double_t fA;      // parameters of right line : y = A*z + B
  Double_t fB;      // system where zero point is IP.
  Double_t fThetaA; // angle coresponding fA - for convinience
  Double_t fTheta;  // theta angle of perependicular to SK module
  // position of towers with differents ieta (1 or 2) -  4-nov-04
  TVector2 fOK1;
  TVector2 fOK2;

 public:
  ClassDef(AliEMCALShishKebabTrd1Module,0) // Turned Shish-Kebab module 
};

#endif
/* To do
 1. Insert position the center of towers - 2 additional TVector2
 */
