#ifndef ALIGENPHOSLIB_H
#define ALIGENPHOSLIB_H
#include <TObject.h>
#include "GenTypeDefs.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliGenPHOSlib :
public TObject
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
    
    typedef Double_t (*GenFunc)  (Double_t *, Double_t *dummy);
    typedef Int_t    (*GenFuncIp)();    
    static GenFunc   GetPt(Param_t param);
    static GenFunc   GetY(Param_t param);
    static GenFuncIp GetIp(Param_t param);    
    ClassDef(AliGenPHOSlib,1) // Library providing y and pT parameterisations
};
#endif







