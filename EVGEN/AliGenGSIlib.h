#ifndef ALIGENGSILIB_H
#define ALIGENGSILIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Implementation of AliGenLib 
// using GSI specific parameterisations.
// Responsible: Andres.Sandoval@cern.ch

#include "AliGenLib.h"
class TRandom;

class AliGenGSIlib :public AliGenLib {
 public:
    enum constants{kUpsilon};
// Upsilon RITMAN   
    static Double_t PtUpsilonRitman( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonRitman(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilonRitman(TRandom *ran);
// Upsilon Karel
    static Double_t PtUpsilonKarel( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonKarel(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilonKarel(TRandom *ran);
// YpsMUON
    static Double_t PtUpsilonMUON( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonMUON(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilonMUON(TRandom *ran);

//
    typedef Double_t (*GenFunc)  (Double_t *, Double_t *);
    typedef Int_t    (*GenFuncIp)(TRandom *ran);
    
    GenFunc   GetPt(Int_t param, const char * tname=0);
    GenFunc   GetY(Int_t param,  const char * tname=0);
    GenFuncIp GetIp(Int_t param, const char * tname=0);    
    static void SetDebug(Bool_t debug){fgDebug=debug;}
private:
    static Bool_t fgDebug;  // Debug flag 
  ClassDef(AliGenGSIlib,0)
};

#endif







