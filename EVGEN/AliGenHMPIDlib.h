#ifndef ALIGENHMPIDLIB_H
#define ALIGENHMPIDLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

// Library class for particle pt and y distributions used for 
// HMPIDsimulations.
// To be used with AliGenParam.
//
// Author: Annalisa Mastroserio <Annalisa.Mastroserio@ba.infn.it>

#include "AliGenLib.h"

class TRandom;

class AliGenHMPIDlib :public AliGenLib {

 public:
  enum EPartId {kPhi=333};

  //Getters
    
  GenFunc   GetPt(Int_t iPID, const char * sForm=0) const;
  GenFunc   GetY (Int_t iPID, const char * sForm=0) const;
  GenFuncIp GetIp(Int_t iPID, const char * sForm=0) const;    
 private:

//Pi+
  static Int_t IpPiPlus(TRandom *ran);
  static Double_t PtPiPlusFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtPiPlusExp (const Double_t *px, const Double_t *dummy);
  static Double_t YPiPlusFlat (const Double_t *py, const Double_t *dummy);

//Pi-
  static Int_t IpPiMinus(TRandom *ran);
  static Double_t PtPiMinusFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtPiMinusExp (const Double_t *px, const Double_t *dummy);
  static Double_t YPiMinusFlat (const Double_t *py, const Double_t *dummy);

//K+
  static Int_t IpKPlus(TRandom *ran);
  static Double_t PtKPlusFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtKPlusExp (const Double_t *px, const Double_t *dummy);
  static Double_t YKPlusFlat (const Double_t *py, const Double_t *dummy);

//K-
  static Int_t IpKMinus(TRandom *ran);
  static Double_t PtKMinusFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtKMinusExp (const Double_t *px, const Double_t *dummy);
  static Double_t YKMinusFlat (const Double_t *py, const Double_t *dummy);

// K0_s
  static Int_t IpK0s(TRandom *ran);
  static Double_t PtK0sFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtK0sExp (const Double_t *px, const Double_t *dummy);
  static Double_t YK0sFlat (const Double_t *py, const Double_t *dummy);

// Phi(1020)
  static Int_t IpPhi(TRandom *ran);
  static Double_t PtPhiFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtPhiExp (const Double_t *px, const Double_t *dummy);
  static Double_t YPhiFlat (const Double_t *py, const Double_t *dummy);

//Proton
  static Int_t IpProton(TRandom *ran);
  static Double_t PtProtonFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtProtonExp (const Double_t *px, const Double_t *dummy);
  static Double_t YProtonFlat (const Double_t *py, const Double_t *dummy);

//ProtonBar
  static Int_t IpProtonBar(TRandom *ran);
  static Double_t PtProtonBarFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtProtonBarExp (const Double_t *px, const Double_t *dummy);
  static Double_t YProtonBarFlat (const Double_t *py, const Double_t *dummy);

// Lambda
  static Int_t IpLambda(TRandom *ran);
  static Double_t PtLambdaFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtLambdaExp (const Double_t *px, const Double_t *dummy);
  static Double_t YLambdaFlat (const Double_t *py, const Double_t *dummy);

// LambdaBar
  static Int_t IpLambdaBar(TRandom *ran);
  static Double_t PtLambdaBarFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtLambdaBarExp (const Double_t *px, const Double_t *dummy);
  static Double_t YLambdaBarFlat (const Double_t *py, const Double_t *dummy);

  ClassDef(AliGenHMPIDlib,0)
};

#endif







