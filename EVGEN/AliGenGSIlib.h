#ifndef ALIGENGSILIB_H
#define ALIGENGSILIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TROOT.h>
#include "GenTypeDefs.h"

class AliGenGSIlib :public TObject{
 public:

// Upsilon RITMAN   
    static Double_t PtUpsilonRitman( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonRitman(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilonRitman();
// Upsilon Karel
    static Double_t PtUpsilonKarel( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonKarel(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilonKarel();
// YpsMUON
    static Double_t PtUpsilonMUON( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonMUON(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilonMUON();

//
    typedef Double_t (*GenFunc)  (Double_t *, Double_t *);
    typedef Int_t    (*GenFuncIp)();    
    static GenFunc   GetPt(Param_t param,const char * tname=0);
    static GenFunc   GetY(Param_t param,const char * tname=0);
    static GenFuncIp GetIp(Param_t param,const char *tname=0);    
    static void SetDebug(Bool_t debug){fgDebug=debug;}
private:
    static Bool_t fgDebug;  // Debug flag 
  ClassDef(AliGenGSIlib,1)
};

#endif







