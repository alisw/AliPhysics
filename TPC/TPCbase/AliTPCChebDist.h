#ifndef ALITPCCHEBDIST_H
#define ALITPCCHEBDIST_H


/*****************************************************************************
 *         Wrapper for 2D->ND Chebishev parameterizations in TPC volume      *
 *                                                                           *
 *  Class is similar to its base class AliTPCChebCorr, but the evaluation    *
 *  is obtained as a linear interpolation of the evaluations on the slices   *
 *  encompassing the queried X (in sector coordinates). These slices are     *
 *  are not the TPC rows but just a grid for evaluation, their number        *
 *  should be in general > kNRows of the TPC                                 *
 *                                                                           *
 *         Author: ruben.shahoyan@cern.ch                                    *
 *****************************************************************************/

#include <TNamed.h>
#include "AliTPCChebCorr.h"

class AliTPCChebDist : public AliTPCChebCorr
{
 public:
  //
 public:
  //
  AliTPCChebDist();
  AliTPCChebDist(const char* name, const char* title, int nps=1,int nzs=1, float zmaxAbs=250);
  virtual ~AliTPCChebDist() {}
  //
  Float_t  GetXMin()                             const {return fXMin;}
  Float_t  GetXMax()                             const {return fXMax;}  
  //
  void     Eval(int sector, float x, float y2x, float z,float *distortion) const;
  void     Eval(int sector, float xtz[3], float *distortion)               const;
  //
  virtual  Bool_t   IsCorrection()               const {return kFALSE;}
  virtual  Bool_t   IsDistorttion()              const {return kTRUE;}
  //
 protected:
  Int_t    X2Slice(float x) const;
  Float_t  Slice2X(int ix)  const;
  //
 protected:
  //
  Float_t  fXMin;                                       // min X
  Float_t  fXMax;                                       // max X
  Float_t  fDX;                                         // X step
  Float_t  fDXInv;                                      // inverse of X step
  //
  static Float_t fgRMinTPC;                             // def. min radius
  static Float_t fgRMaxTPC;                             // def. max radius
  static Int_t   fgNSlices;                             // def. number of slices in X
 private:
  AliTPCChebDist(const AliTPCChebDist& src);            // dummy
  AliTPCChebDist& operator=(const AliTPCChebDist& rhs); // dummy
  //
  ClassDef(AliTPCChebDist,1)
};

//_________________________________________________________________
inline Float_t AliTPCChebDist::Slice2X(Int_t ix) const
{
  // get the lower slice encompacing given X, except if the X is outside of the fid. range
  return fXMin + ix*fDX;
}

//_________________________________________________________________
inline Int_t AliTPCChebDist::X2Slice(float x) const
{
  // get the lower slize covering given X
  int ix = (x-fXMin)*fDXInv;
  if      (ix<0)         ix = 0;          
  else if (ix>=fNRows)   ix = fNRows-1;
  return ix;
}

//____________________________________________________________________
inline void AliTPCChebDist::Eval(int sector, float xtz[3], float *distortion) const
{
  // Calculate distortion for point with x,y,z sector corrdinates
  // Sector is in 0-71 ROC convention, to check Zs outlying from the sector
  Eval(xtz[0],&xtz[1],distortion);
  //
}


#endif
