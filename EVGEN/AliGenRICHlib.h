#ifndef ALIGENRICHLIB_H
#define ALIGENRICHLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TPDGCode.h>
#include <AliGenLib.h>

class TRandom;

class AliGenRICHlib :public AliGenLib {

 public:
   enum EPartId {kPhi=333};

// Phi(1020)
  static Int_t IpPhi(TRandom *ran);
  static Double_t PtPhiFlat(Double_t *px, Double_t *dummy);
  static Double_t PtPhiExp (Double_t *px, Double_t *dummy);
  static Double_t YPhiFlat (Double_t *py, Double_t *dummy);

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


// K0_s
  static Int_t IpK0s(TRandom *ran);
  static Double_t PtK0sFlat(Double_t *px, Double_t *dummy);
  static Double_t PtK0sExp (Double_t *px, Double_t *dummy);
  static Double_t YK0sFlat (Double_t *py, Double_t *dummy);


  typedef Double_t (*GenFunc)  (Double_t *, Double_t *);
  typedef Int_t    (*GenFuncIp)(TRandom *ran);

  //Getters
    
  GenFunc   GetPt(Int_t iPID, const char * sForm=0) const;
  GenFunc   GetY (Int_t iPID, const char * sForm=0) const;
  GenFuncIp GetIp(Int_t iPID, const char * sForm=0) const;    

  ClassDef(AliGenRICHlib,0)
};

#endif







