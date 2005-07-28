#ifndef ALIGENSTRANGELIB_H
#define ALIGENSTRANGELIB_H
#include "AliGenLib.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//======================================================================
//  AliGenSTRANGElib class contains parameterizations of the
//  kaon, phi and hyperon (Lambda, Anti-Lambda, Xi, anti-Xi, Omega,
//  anti-Omega)  for the PPR study of the strange particle production. 
//
//  Rocco CALIANDRO. Rosa Anna FINI, Tiziano VIRGILI
//======================================================================

class TRandom;

class AliGenSTRANGElib :
public AliGenLib
{
 public:
    enum constants{kKaon, kPhi, kLambda, kXiMinus, kOmegaMinus, kLambda1520};
    GenFunc   GetPt(Int_t param, const char* tname=0) const;
    GenFunc   GetY (Int_t param, const char* tname=0) const;
    GenFuncIp GetIp(Int_t param, const char* tname=0) const;    
 private:
// pions
    static Double_t PtPion(Double_t *px, Double_t *dummy);
    static Double_t PtScal(Double_t pt, Int_t np);
// kaons
    static Double_t PtKaon(Double_t *px, Double_t *dummy);
    static Double_t YKaon( Double_t *py, Double_t *dummy);
    static Int_t    IpKaon(TRandom* ran);
// phis
    static Double_t PtPhi(Double_t *px, Double_t *dummy);
    static Double_t YPhi( Double_t *py, Double_t *dummy);
    static Int_t    IpPhi(TRandom* ran);
// lambda
    static Double_t PtLambda(Double_t *px, Double_t *dummy);
    static Double_t YLambda( Double_t *py, Double_t *dummy);
    static Int_t    IpLambda(TRandom *ran);
// Ximinus
    static Double_t PtXiMinus(Double_t *px, Double_t *dummy);
    static Double_t YXiMinus( Double_t *py, Double_t *dummy);
    static Int_t    IpXiMinus(TRandom *ran);
// Omegaminus
    static Double_t PtOmegaMinus(Double_t *px, Double_t *dummy);
    static Double_t YOmegaMinus( Double_t *py, Double_t *dummy);
    static Int_t    IpOmegaMinus(TRandom *ran);
// Lambda(1520)
    static Double_t PtLambda1520(Double_t *px, Double_t *dummy);
    static Double_t YLambda1520(Double_t *py, Double_t *dummy);
    static Int_t    IpLambda1520(TRandom *ran);
    
    ClassDef(AliGenSTRANGElib,0) // Library providing y and pT parameterisations
};
#endif














