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
  static Double_t PtPiPlusFlat(Double_t *px, Double_t *dummy);
  static Double_t PtPiPlusExp (Double_t *px, Double_t *dummy);
  static Double_t YPiPlusFlat (Double_t *py, Double_t *dummy);

//Pi-
  static Int_t IpPiMinus(TRandom *ran);
  static Double_t PtPiMinusFlat(Double_t *px, Double_t *dummy);
  static Double_t PtPiMinusExp (Double_t *px, Double_t *dummy);
  static Double_t YPiMinusFlat (Double_t *py, Double_t *dummy);

//K+
  static Int_t IpKPlus(TRandom *ran);
  static Double_t PtKPlusFlat(Double_t *px, Double_t *dummy);
  static Double_t PtKPlusExp (Double_t *px, Double_t *dummy);
  static Double_t YKPlusFlat (Double_t *py, Double_t *dummy);

//K-
  static Int_t IpKMinus(TRandom *ran);
  static Double_t PtKMinusFlat(Double_t *px, Double_t *dummy);
  static Double_t PtKMinusExp (Double_t *px, Double_t *dummy);
  static Double_t YKMinusFlat (Double_t *py, Double_t *dummy);

// K0_s
  static Int_t IpK0s(TRandom *ran);
  static Double_t PtK0sFlat(Double_t *px, Double_t *dummy);
  static Double_t PtK0sExp (Double_t *px, Double_t *dummy);
  static Double_t YK0sFlat (Double_t *py, Double_t *dummy);

// Phi(1020)
  static Int_t IpPhi(TRandom *ran);
  static Double_t PtPhiFlat(Double_t *px, Double_t *dummy);
  static Double_t PtPhiExp (Double_t *px, Double_t *dummy);
  static Double_t YPhiFlat (Double_t *py, Double_t *dummy);

//Proton
  static Int_t IpProton(TRandom *ran);
  static Double_t PtProtonFlat(Double_t *px, Double_t *dummy);
  static Double_t PtProtonExp (Double_t *px, Double_t *dummy);
  static Double_t YProtonFlat (Double_t *py, Double_t *dummy);

//ProtonBar
  static Int_t IpProtonBar(TRandom *ran);
  static Double_t PtProtonBarFlat(Double_t *px, Double_t *dummy);
  static Double_t PtProtonBarExp (Double_t *px, Double_t *dummy);
  static Double_t YProtonBarFlat (Double_t *py, Double_t *dummy);

// Lambda
  static Int_t IpLambda(TRandom *ran);
  static Double_t PtLambdaFlat(Double_t *px, Double_t *dummy);
  static Double_t PtLambdaExp (Double_t *px, Double_t *dummy);
  static Double_t YLambdaFlat (Double_t *py, Double_t *dummy);

// LambdaBar
  static Int_t IpLambdaBar(TRandom *ran);
  static Double_t PtLambdaBarFlat(Double_t *px, Double_t *dummy);
  static Double_t PtLambdaBarExp (Double_t *px, Double_t *dummy);
  static Double_t YLambdaBarFlat (Double_t *py, Double_t *dummy);

  ClassDef(AliGenHMPIDlib,0)
};

#endif







