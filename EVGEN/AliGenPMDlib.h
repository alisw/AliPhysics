#ifndef ALIGENPMDLIB_H
#define ALIGENPMDLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenLib.h"

class AliGenPMDlib :
public AliGenLib
{
 public:
// Neutral pions
    static Double_t PtPi0(Double_t *px, Double_t *dummy);
    static Double_t PtScal(Double_t pt, Int_t np);
    static Double_t YPi0( Double_t *py, Double_t *dummy);
    static Int_t    IpPi0();
// Etas
    static Double_t PtEta(Double_t *px, Double_t *dummy);
    static Double_t YEta( Double_t *py, Double_t *dummy);
    static Int_t    IpEta();
//
    GenFunc   GetPt(Param_t param, const char* tname=0);
    GenFunc   GetY (Param_t param, const char* tname=0);
    GenFuncIp GetIp(Param_t param, const char* tname=0);    
    ClassDef(AliGenPMDlib,1) // Library providing y and pT parameterisations
};
#endif







