#ifndef ALIGENMUONLIB_H
#define ALIGENMUONLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenLib.h"

class AliGenMUONlib :
public AliGenLib
{
 public:
// pions
    static Double_t PtPion(Double_t *px, Double_t *dummy);
    static Double_t PtScal(Double_t pt, Int_t np);
    static Double_t YPion( Double_t *py, Double_t *dummy);
    static Int_t    IpPion();
// kaons
    static Double_t PtKaon(Double_t *px, Double_t *dummy);
    static Double_t YKaon( Double_t *py, Double_t *dummy);
    static Int_t    IpKaon();
// Phi
    static Double_t PtPhi( Double_t *px, Double_t *dummy);
    static Double_t YPhi( Double_t *px, Double_t *dummy);
    static Int_t    IpPhi();
// J/Psi     
    static Double_t PtJpsi( Double_t *px, Double_t *dummy);
    static Double_t YJpsi(Double_t *py, Double_t *dummy);
    static Int_t    IpJpsi();
// Upsilon    
    static Double_t PtUpsilon( Double_t *px, Double_t *dummy );
    static Double_t YUpsilon(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilon();
//
// Charm    
    static Double_t PtCharm( Double_t *px, Double_t *dummy );
    static Double_t YCharm(Double_t *py, Double_t *dummy);
    static Int_t    IpCharm();
//
// Beauty
    static Double_t PtBeauty( Double_t *px, Double_t *dummy );
    static Double_t YBeauty(Double_t *py, Double_t *dummy);
    static Int_t    IpBeauty();
//
    GenFunc   GetPt(Param_t param, const char* tname=0);
    GenFunc   GetY (Param_t param, const char* tname=0);
    GenFuncIp GetIp(Param_t param, const char* tname=0);    
    ClassDef(AliGenMUONlib,0) // Library providing y and pT parameterisations
};
#endif







