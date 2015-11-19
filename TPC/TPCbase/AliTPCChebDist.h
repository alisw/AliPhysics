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
#include "AliTPCChecDist.h"

class AliTPCChebDist : public AliTPCChecDist
{
 public:
  //
 public:
  //
  AliTPCChebDist();
  AliTPCChebDist(const char* name, const char* title, int nps=1,int nzs=1, float zmaxAbs=250);
  virtual ~AliTPCChebDist();
  //
  void     Print(const Option_t* opt="")         const;
  void     Eval(int sector, float x, float y2x, float z,float *distortion) const;
  void     Eval(int sector, float xtz[2], float *distortion)       const;
  //
 protected:
  //
 protected:
  //
  AliTPCChebDist(const AliTPCChebDist& src);            // dummy
  AliTPCChebDist& operator=(const AliTPCChebDist& rhs); // dummy
  //
  ClassDef(AliTPCChebDist,2)
};

//_________________________________________________________________
inline const AliCheb2DStack* AliTPCChebDist::GetParam(int sector, float y2x, float z) const
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
inline void AliTPCChebDist::Eval(int sector, float x, float y2x, float z, float *distortion) const
{
  // Calculate distortion for point with x,y,z sector corrdinates
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  GetParam(sector,y2x,z)->Eval(row, tz, corr);
  //
}

//____________________________________________________________________
inline void AliTPCChebDist::Eval(int sector, int row, float tz[2], float *corr) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // Sector is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  GetParam(sector,tz[0],tz[1])->Eval(row, tz, corr);
  //
}

//____________________________________________________________________
inline Float_t AliTPCChebDist::Eval(int sector, int row, float y2x, float z, int dimOut) const
{
  // Calculate dimOut-th correction for point with x,y,z sector corrdinates
  // Sector/row is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  return GetParam(sector,y2x,z)->Eval(row, dimOut, tz);
  //
}

//____________________________________________________________________
inline Float_t AliTPCChebDist::Eval(int sector, int row, float tz[2], int dimOut) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // Sector is in 0-71 ROC convention, to check Zs outlying from the sector
  if (sector>kMaxIROCSector) row += kNRowsIROC;   // we are in OROC
  return GetParam(sector,tz[0],tz[1])->Eval(row, dimOut, tz);
  //
}


#endif
