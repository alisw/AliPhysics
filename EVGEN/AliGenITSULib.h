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

  enum EPartId {kLb=5122,kLc=4122,kXi_c = 4232,kBplus = 521,kBzero = 511,kDs=431,kDplus=411, kOmega_ccc=4444,kXi_czero = 4132, kOmega_c=4332, kBs=531};

  //Getters
    
  GenFunc   GetPt(Int_t iPID, const char * sForm=0) const;
  GenFunc   GetY (Int_t iPID, const char * sForm=0) const;
  GenFuncIp GetIp(Int_t iPID, const char * sForm=0) const;    

 private:

  static Int_t IpLcPlus(TRandom * /*ran*/)  {return     (int)kLc;}
  static Int_t IpLcMinus(TRandom * /*ran*/) {return    -(int)kLc;}
  static Int_t IpLb(TRandom * /*ran*/)      {return     (int)kLb;}
  static Int_t IpLbBar(TRandom * /*ran*/)   {return    -(int)kLb;}
  static Int_t IpXicZero(TRandom * /*ran*/)     {return   (int)kXi_czero;}
  static Int_t IpXicZeroBar(TRandom * /*ran*/)  {return  -(int)kXi_czero;}
  static Int_t IpXic(TRandom * /*ran*/)     {return   (int)kXi_c;}
  static Int_t IpXicBar(TRandom * /*ran*/)  {return  -(int)kXi_c;}
  static Int_t IpOmegac(TRandom * /*ran*/)     {return   (int)kOmega_c;}
  static Int_t IpOmegacBar(TRandom * /*ran*/)  {return  -(int)kOmega_c;}
  static Int_t IpBPlus(TRandom * /*ran*/)   {return      (int)kBplus;}
  static Int_t IpBMinus(TRandom * /*ran*/)  {return     -(int)kBplus;}
  static Int_t IpBs(TRandom * /*ran*/)  {return         (int)kBs;}
  static Int_t IpBsBar(TRandom * /*ran*/)  {return         -(int)kBs;}
  static Int_t IpB0(TRandom * /*ran*/)  {return         (int)kBzero;}
  static Int_t IpB0Bar(TRandom * /*ran*/)  {return         -(int)kBzero;}
  static Int_t IpDsPlus(TRandom * /*ran*/)  {return     (int)kDs;}
  static Int_t IpDsMinus(TRandom * /*ran*/) {return    -(int)kDs;}
  static Int_t IpDPlus(TRandom * /*ran*/)   {return  (int)kDplus;}
  static Int_t IpDMinus(TRandom * /*ran*/)  {return  -(int)kDplus;}

  static Int_t IpOmegaccc(TRandom * /*ran*/);
  
  static Double_t PtFlat(const Double_t * /*px*/, const Double_t * /*dummy*/) {return 1;}
  static Double_t YFlat (const Double_t * /*py*/, const Double_t * /*dummy*/) {return 1;}

  static Double_t PtLbDist (const Double_t *px, const Double_t *dummy);
  static Double_t PtLcDist (const Double_t *px, const Double_t *dummy);
  static Double_t PtBDist  (const Double_t *px, const Double_t *dummy);



  ClassDef(AliGenITSULib,1)
    };

#endif







