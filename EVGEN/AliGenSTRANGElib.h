#ifndef ALIGENSTRANGELIB_H
#define ALIGENSTRANGELIB_H
#include "AliGenLib.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

class TRandom;

class AliGenSTRANGElib :
public AliGenLib
{
 public:
    enum constants{kKaon, kPhi, kLambda, kXiMinus, kOmegaMinus};
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
    
    GenFunc   GetPt(Int_t param, const char* tname=0);
    GenFunc   GetY (Int_t param, const char* tname=0);
    GenFuncIp GetIp(Int_t param, const char* tname=0);    
    ClassDef(AliGenSTRANGElib,0) // Library providing y and pT parameterisations
};
#endif







