#ifndef ALIGENPHOSLIB_H
#define ALIGENPHOSLIB_H
#include "AliGenLib.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliGenPHOSlib :
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
// etas
    static Double_t PtEta(Double_t *px, Double_t *dummy);
    static Double_t YEta( Double_t *py, Double_t *dummy);
    static Int_t    IpEta();
// omegas
    static Double_t PtOmega(Double_t *px, Double_t *dummy);
    static Double_t YOmega( Double_t *py, Double_t *dummy);
    static Int_t    IpOmega();
// etaprime
    static Double_t PtEtaprime(Double_t *px, Double_t *dummy);
    static Double_t YEtaprime( Double_t *py, Double_t *dummy);
    static Int_t    IpEtaprime();
// phis
    static Double_t PtPhi(Double_t *px, Double_t *dummy);
    static Double_t YPhi( Double_t *py, Double_t *dummy);
    static Int_t    IpPhi();
// baryons
    static Double_t PtBaryon(Double_t *px, Double_t *dummy);
    static Double_t YBaryon( Double_t *py, Double_t *dummy);
    static Int_t    IpBaryon();
    
    GenFunc   GetPt(Param_t param, const char* tname=0);
    GenFunc   GetY (Param_t param, const char* tname=0);
    GenFuncIp GetIp(Param_t param, const char* tname=0);    
    ClassDef(AliGenPHOSlib,0) // Library providing y and pT parameterisations
};
#endif







