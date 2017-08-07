#ifndef ALIMAGFAST_H
#define ALIMAGFAST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Fast polynomial parametrization of Alice magnetic field, to be used for reconstruction.
// Solenoid part fitted by Shuto Yamasaki from AliMagWrapCheb in the |Z|<260Interface and R<500 cm 
// Dipole part: to do
//
// Author: ruben.shahoyan@cern.ch
//
#include <TObject.h>

class AliMagFast : public TObject
{

 public:
  enum {kNSolRRanges=5, kNSolZRanges=22, kNQuadrants=4};
  enum {kX,kY,kZ};
  
  struct SolParam { float mParBxyz[3][20];};
  typedef SolParam SolParam_t;

  AliMagFast(const char* inpFName=0);
  AliMagFast(Float_t factor, Int_t nomField = 5, const char* inpFmt="$(ALICE_ROOT)/data/maps/sol%dk.txt");
  AliMagFast(const AliMagFast& src);
  AliMagFast& operator=(const AliMagFast& src);
  
  Bool_t LoadData(const char* inpFName);
  virtual ~AliMagFast() {}

  Bool_t Field(const double xyz[3], double bxyz[3]) const;
  Bool_t GetBz(const double xyz[3], double& bz)     const;
  Bool_t Field(const float  xyz[3], float bxyz[3])  const;
  Bool_t GetBz(const float  xyz[3], float& bz)      const;

  void    SetFactorSol(float v=1.f)                       {fFactorSol = v;}
  Float_t GetFactorSol()                            const {return fFactorSol;}
  
 protected:

  Bool_t GetSegment(const float xyz[3], int& zSeg,int &rSeg, int &quadrant) const;  
  static const float fgkSolR2Max[kNSolRRanges];       // Rmax2 of each range
  static const float fgkSolZMax;                      // max |Z| for solenoid parametrization


  int GetQuadrant(float x,float y) const
  {
    // get point quadrant
    return y>0 ? (x>0 ? 0:1) : (x>0 ? 3:2);
  }
  
  float CalcPol(const float* cf, float x,float y, float z) const;

  Float_t fFactorSol; // scaling factor
  SolParam_t fSolPar[kNSolRRanges][kNSolZRanges][kNQuadrants];
  
  ClassDef(AliMagFast,1)
};

inline float AliMagFast::CalcPol(const float* cf, float x,float y, float z) const
{

  /** calculate polynomial
   *   cf[0] + cf[1]*x + cf[2]*y + cf[3]*z + cf[4]*xx + cf[5]*xy + cf[6]*xz + cf[7]*yy + cf[8]*yz + cf[9]*zz +
   *   cf[10]*xxx + cf[11]*xxy + cf[12]*xxz + cf[13]*xyy + cf[14]*xyz + cf[15]*xzz + cf[16]*yyy + cf[17]*yyz + cf[18]*yzz + cf[19]*zzz
  **/

    float val = cf[0] +
    x*(cf[1] + x*(cf[4] + x*cf[10] + y*cf[11] + z*cf[12]) + y*(cf[5]+z*cf[14]) ) +
    y*(cf[2] + y*(cf[7] + x*cf[13] + y*cf[16] + z*cf[17]) + z*(cf[8]) ) +
    z*(cf[3] + z*(cf[9] + x*cf[15] + y*cf[18] + z*cf[19]) + x*(cf[6]) );
    
  return val;
}


#endif
