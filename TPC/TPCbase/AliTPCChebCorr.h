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
#include <time.h>
#include "AliCheb2DStack.h"

class AliTPCChebCorr : public TNamed
{
 public:
  enum {kFieldAny, kFieldPos, kFieldNeg, kFieldZero};
  enum {kNSectors=18,kNSectorsIROC=2*kNSectors,kNRows=159,kNRowsIROC=63,kMaxIROCSector=kNSectorsIROC-1};
  enum {kParamDone=BIT(14), // parameterization done
	kUseParF=BIT(15),   // if ON - internal FLOAT representation, otherwise - SHORT
	kUseZ2R=BIT(16),    // 2nd dimension parameterizes Z/R
	kTimeDependent=BIT(17) // flag to signal time-dependent objects to follow
  };
  //
 public:
  //
  AliTPCChebCorr();
  AliTPCChebCorr(const char* name, const char* title, int nps=1,int nzs=1, float zmaxAbs=250, float deadZone=1.5, const float *xi=0);
  virtual ~AliTPCChebCorr();
  Int_t    GetFieldType()                        const {return fFieldType;}
  void     SetFieldType(Char_t t=kFieldAny)            {fFieldType = t;}
  void     Parameterize(stFun_t fun,int dimOut,const int np[2], const float *prec=0);
  void     Parameterize(stFun_t fun,int dimOut,const int np[][2], const float *prec=0);
  void     SetBinning(int nps=1,int nzs=1, float zmxAbs=250);
  Bool_t   GetUseFloatPrec()                     const {return TestBit(kUseParF);}
  Bool_t   GetUseShortPrec()                     const {return !TestBit(kUseParF);}
  Bool_t   SetUseFloatPrec(Bool_t v)                   {SetBit(kUseParF,v); return v;}
  Bool_t   GetUseZ2R()                           const {return TestBit(kUseZ2R);}
  void     SetUseZ2R(Bool_t v=kTRUE)                   {SetBit(kUseZ2R,v);}
  Bool_t   GetTimeDependent()                    const {return TestBit(kTimeDependent);}
  void     SetTimeDependent(Bool_t v=kTRUE)            {SetBit(kTimeDependent,v);}
  Float_t  GetZMin()                             const {return -fZMaxAbs;}
  Float_t  GetZMax()                             const {return fZMaxAbs;}
  Int_t    GetNStacksZ()                         const {return fNStacksZ;}
  Int_t    GetNStacksSector()                    const {return fNStacksSect;}
  Int_t    GetNRows()                            const {return fNRows;}
  Float_t  GetDeadZone()                         const {return fDeadZone;}
  //
  const AliCheb2DStack* GetParam(int id)         const {return (const AliCheb2DStack*) fParams ?  fParams[id] : 0;}
  const AliCheb2DStack* GetParam(int sector, float y2x, float z) const;
  //
  time_t   GetTimeStampStart()                   const {return fTimeStampStart;}
  time_t   GetTimeStampEnd()                     const {return fTimeStampEnd;}
  time_t   GetTimeStampCenter()                  const {return (fTimeStampStart+fTimeStampEnd)/2;}
  void     SetTimeStampStart(time_t t=0)               {fTimeStampStart = t;}
  void     SetTimeStampEnd(time_t t=0xffffffff)        {fTimeStampEnd = t;}
  //
  void     Print(const Option_t* opt="")         const;
  void     Eval(int sector, int row, float y2x, float z,float *corr) const;
  void     Eval(int sector, int row, float tz[2], float *corr)       const;
  Float_t  Eval(int sector, int row, float y2x, float z, int dimOut) const;
  Float_t  Eval(int sector, int row, float tz[2], int dimOut)        const;
  void     Init();
  static   float GetMaxY2X()                    {return fgkY2XHSpan;}
  static const float* GetPadRowX()              {return fgkPadRowX;}
  //
  virtual  Bool_t   IsCorrection()               const {return kTRUE;}
  virtual  Bool_t   IsDistortion()               const {return kFALSE;}
  //
 protected:
  //
  int      GetParID(int iz,int isect,int istack) const {return (iz*kNSectors+isect)*fNStacksSect+istack;}
  //
 protected:
  Char_t   fFieldType;              // info about the field type
  Int_t    fNRows;                  // number of slices along the radius (e.g. rows)
  Int_t    fNStacksSect;            // number of stacks per sector in phi
  Int_t    fNStacksZSect;           // number of stacks per sector (side) in Z 
  Int_t    fNStacksZ;               // number of stacks in Z
  Int_t    fNStacks;                // total number of stacks
  Float_t  fZMaxAbs;                // zmax abs
  //
  time_t   fTimeStampStart;         // time stamp for start of validity
  time_t   fTimeStampEnd;           // time stamp for end of validity
  //
  // precalculated parameters
  Float_t  fZScaleI;                // 1/Zspan of single stack
  Float_t  fY2XScaleI;              // 1/Y2Xspan of single stack
  //
  Float_t  fDeadZone;               // dead zone in cm
  Float_t* fRowXI;                  //[fNRows]  // 1/X of each row is dead zone to be used
  //
  AliCheb2DStack** fParams;         //[fNStacks] set of AliCheb2DStack parameterizations
  //
  static const float fgkY2XHSpan;   // half span of sector
  static const float fgkPadRowX[];  // nominal rows
  static const char* fgkFieldTypeName[]; // names of field types
 protected:
  //
  AliTPCChebCorr(const AliTPCChebCorr& src);            // dummy
  AliTPCChebCorr& operator=(const AliTPCChebCorr& rhs); // dummy
  //
  ClassDef(AliTPCChebCorr,4)
};

//_________________________________________________________________
inline const AliCheb2DStack* AliTPCChebCorr::GetParam(int sector, float y2x, float z) const
{
  // Find appropriate param. Sector is in ROC0-71 conventions
  int iz = (z+fZMaxAbs)*fZScaleI, side = (sector/kNSectors)&0x1;
  // correct for eventual Z calculated in wrong ROC
  if (side)  {if (iz>=fNStacksZSect) iz = fNStacksZSect-1;} // C side
  else       {if (iz<fNStacksZSect)  iz = fNStacksZSect;}   // A side
  if (iz<0) iz=0; else if (iz>=fNStacksZ) iz=fNStacksZ-1;
  int is = (y2x+fgkY2XHSpan)*fY2XScaleI;
  if (is<0) is=0; else if (is>=fNStacksSect) is=fNStacksSect-1;
  return GetParam(GetParID(iz,sector%kNSectors,is));
  //
}

//____________________________________________________________________
inline void AliTPCChebCorr::Eval(int sector, int row, float y2x, float z, float *corr) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // Sector/row is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  const AliCheb2DStack* par = GetParam(sector,y2x,z);
  if (par) par->Eval(row, tz, corr);
  //
}

//____________________________________________________________________
inline void AliTPCChebCorr::Eval(int sector, int row, float tz[2], float *corr) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // Sector is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  const AliCheb2DStack* par = GetParam(sector,tz[0],tz[1]);
  if (par) par->Eval(row, tz, corr);
  //
}

//____________________________________________________________________
inline Float_t AliTPCChebCorr::Eval(int sector, int row, float y2x, float z, int dimOut) const
{
  // Calculate dimOut-th correction for point with x,y,z sector corrdinates
  // Sector/row is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  const AliCheb2DStack* par = GetParam(sector,y2x,z);
  return par ? par->Eval(row, dimOut, tz) : 0;
  //
}

//____________________________________________________________________
inline Float_t AliTPCChebCorr::Eval(int sector, int row, float tz[2], int dimOut) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // Sector is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  const AliCheb2DStack* par = GetParam(sector,tz[0],tz[1]);
  return par ? par->Eval(row, dimOut, tz) : 0;
  //
}


#endif
