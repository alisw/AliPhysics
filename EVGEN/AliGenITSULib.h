#ifndef ALIGENITSULIB_H
#define ALIGENITSULIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Library class for particle pt and y distributions used for 
// LambdaB simulations.
// To be used with AliGenParam.
//
// Author: Annalisa Mastroserio <Annalisa.Mastroserio@cern.ch>

#include "AliGenLib.h"

class TRandom;

class AliGenITSULib :public AliGenLib {

 public:

  enum EPartId {kLb=5122,kLc=4122,kXi_c = 4232, kB = 521};

  //Getters
    
  GenFunc   GetPt(Int_t iPID, const char * sForm=0) const;
  GenFunc   GetY (Int_t iPID, const char * sForm=0) const;
  GenFuncIp GetIp(Int_t iPID, const char * sForm=0) const;    

 private:

  static Int_t IpLcPlus(TRandom * /*ran*/)  {return   (int)kLc;}
  static Int_t IpLcMinus(TRandom * /*ran*/) {return  -(int)kLc;}
  static Int_t IpLb(TRandom * /*ran*/)      {return   (int)kLb;}
  static Int_t IpLbBar(TRandom * /*ran*/)   {return  -(int)kLb;}
  static Int_t IpXic(TRandom * /*ran*/)     {return (int)kXi_c;}

  static Double_t PtFlat(const Double_t * /*px*/, const Double_t * /*dummy*/) {return 1;}
  static Double_t YFlat (const Double_t * /*py*/, const Double_t * /*dummy*/) {return 1;}

  static Double_t PtLbDist (const Double_t *px, const Double_t *dummy);
  static Double_t PtLcDist (const Double_t *px, const Double_t *dummy);
  static Double_t PtBDist  (const Double_t *px, const Double_t *dummy);



  ClassDef(AliGenITSULib,0)
    };

#endif







