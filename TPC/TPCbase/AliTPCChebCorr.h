#ifndef ALITPCCHEBCORR_H
#define ALITPCCHEBCORR_H


/*****************************************************************************
 *         Wrapper for 2D->ND Chebishev parameterizations in TPC volume      *
 *                                                                           *
 *  Creation requires                                                        *
 *  1) Optional number of patches in Z direction and Z limits                *
 *  2) Optional number of patches per sector                                 *
 *  Once the constructor is called and optional parameters are set, use      *
 *  Parameterize(stFun_t fun,int dimOut,const int np[2],const float* prec)   *
 *  or                                                                       *
 *  Parameterize(stFun_t fun,int dimOut,const int np[][2],const float* prec) *
 *  methods, where                                                           *
 *  1) fun: user function with signature                                     *
 *  void (*fun)(int row, float* tz,float* valND)                             *
 *  where row is 0-158, tz is x2y within the sector corrdinates and          *
 *  valND is returned ND-dim. array of values to parameterize.               *
 *  The function should memorize internally selected sector when called as   *
 *  fun(sector,0,0)                                                          *
 *  2) dimOut: dimensionality of the output                                  *
 *  3) np: number of training points to use in each input dimension,         *
 *  common for all output dimensions (np[2]) or for unique                   *
 *  for each output dimension (np[dimOut][2])                                *
 *  4) optional array prec[dimOut] with requested absolute tolerances        *
 *                                                                           *
 *         Author: ruben.shahoyan@cern.ch                                    *
 *****************************************************************************/

#include <TNamed.h>
#include "AliCheb2DStack.h"

class AliTPCChebCorr : public TNamed
{
 public:
  enum {kNSectors=18,kNRows=159};
  enum {kParamDone=BIT(14), // parameterization done
	kUseParF=BIT(15),   // if ON - internal FLOAT representation, otherwise - SHORT
	kUseZ2R=BIT(16),    // 2nd dimension parameterizes Z/R
	kTimeDependent=BIT(17) // flag to signal time-dependent objects to follow
  };
  //
 public:
  //
  AliTPCChebCorr();
  AliTPCChebCorr(const char* name, const char* title, int nps=1,int nz=2, float zmin=-250, float zmax=250);
  virtual ~AliTPCChebCorr();
  void     Parameterize(stFun_t fun,int dimOut,const int np[2], const float *prec=0);
  void     Parameterize(stFun_t fun,int dimOut,const int np[][2], const float *prec=0);
  void     SetBinning(int nps=1,int nz=2,float zmn=-250,float zmx=250);
  Bool_t   GetUseFloatPrec()                     const {return TestBit(kUseParF);}
  Bool_t   GetUseShortPrec()                     const {return !TestBit(kUseParF);}
  Bool_t   SetUseFloatPrec(Bool_t v)                   {SetBit(kUseParF,v);}
  Bool_t   GetUseZ2R()                           const {return TestBit(kUseZ2R);}
  void     SetUseZ2R(Bool_t v=kTRUE)                   {SetBit(kUseZ2R,v);}
  Bool_t   GetTimeDependent()                    const {return TestBit(kTimeDependent);}
  void     SetTimeDependent(Bool_t v=kTRUE)            {SetBit(kTimeDependent,v);}
  Float_t  GetZMin()                             const {return fZMin;}
  Float_t  GetZMax()                             const {return fZMax;}
  Int_t    GetNStacksZ()                         const {return fNStacksZ;}
  Int_t    GetNStacksSector()                    const {return fNStacksSect;}
  const AliCheb2DStack* GetParam(int id)         const {return (const AliCheb2DStack*) fParams ?  fParams[id] : 0;}
  //
  UInt_t   GetTimeStampStart()                   const {return fTimeStampStart;}
  UInt_t   GetTimeStampEnd()                     const {return fTimeStampEnd;}
  UInt_t   GetTimeStampCenter()                  const {return (fTimeStampStart>>1)+(fTimeStampEnd>>1);}
  void     SetTimeStampStart(UInt_t t=0)               {fTimeStampStart = t;}
  void     SetTimeStampEnd(UInt_t t=0xffffffff)        {fTimeStampEnd = t;}
  //
  void     Print(const Option_t* opt="")         const;
  void     Eval(int sector, int row, float y2x, float z,float *corr) const;
  void     Eval(int sector, int row, float tz[2], float *corr) const;
  static   float GetMaxY2X()                    {return fgkY2XHSpan;}
  //
 protected:
  int      GetParID(int iz,int isect,int istack) const {return (iz*kNSectors+isect)*fNStacksSect+istack;}
  //
 protected:
  Int_t    fNStacksSect;            // number of stacks per sector in phi
  Int_t    fNStacksZ;               // number of stacks in Z
  Int_t    fNStacks;                // total number of stacks
  Float_t  fZMin;                   // zmin
  Float_t  fZMax;                   // zmax
  //
  UInt_t   fTimeStampStart;         // time stamp for start of validity
  UInt_t   fTimeStampEnd;           // time stamp for end of validity
  //
  // precalculated parameters
  Float_t  fZCen;                   // Z center
  Float_t  fZScaleI;                // 1/Zspan of single stack
  Float_t  fY2XScaleI;              // 1/Y2Xspan of single stack
  //
  AliCheb2DStack** fParams;         //[fNStacks] set of AliCheb2DStack parameterizations
  //
  static const float fgkY2XHSpan;   // half span of sector
 protected:
  //
  AliTPCChebCorr(const AliTPCChebCorr& src);            // dummy
  AliTPCChebCorr& operator=(const AliTPCChebCorr& rhs); // dummy
  //
  ClassDef(AliTPCChebCorr,1)
};

#endif
