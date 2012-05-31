#ifndef ALIGENLCLIB_H
#define ALIGENLCLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Library class for particle pt and y distributions used for 
// LambdaC simulations.
// To be used with AliGenParam.
//
// Author: Annalisa Mastroserio <Annalisa.Mastroserio@cern.ch>

#include "AliGenLib.h"

class TRandom;

class AliGenLcLib :public AliGenLib {

 public:
  enum EPartId {kLcPlus=4122,kLcMinus=-4122};

  //Getters
    
  GenFunc   GetPt(Int_t iPID, const char * sForm=0) const;
  GenFunc   GetY (Int_t iPID, const char * sForm=0) const;
  GenFuncIp GetIp(Int_t iPID, const char * sForm=0) const;    
 private:


  static Int_t IpLcPlus(TRandom *ran);
  static Int_t IpLcMinus(TRandom *ran);
  static Double_t PtLcFlat(const Double_t *px, const Double_t *dummy);
  static Double_t PtLcExp (const Double_t *px, const Double_t *dummy);
  static Double_t YLcFlat (const Double_t *py, const Double_t *dummy);

  ClassDef(AliGenLcLib,0)
    };

#endif







