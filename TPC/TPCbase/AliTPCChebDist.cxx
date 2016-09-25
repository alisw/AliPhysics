/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliCheb2DStackF.h"
#include "AliCheb2DStackS.h"
#include "AliTPCChebDist.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliTPCChebDist)

Float_t AliTPCChebDist::fgRMinTPC = 83.65; //RS Make sure these are correct radii to use
Float_t AliTPCChebDist::fgRMaxTPC = 247.7;
Int_t   AliTPCChebDist::fgNSlices = AliTPCChebCorr::kNRows*2+2; // check

//____________________________________________________________________
AliTPCChebDist::AliTPCChebDist()
  : AliTPCChebCorr()
  ,fXMin(fgRMinTPC)
  ,fXMax(fgRMaxTPC)
  ,fDX(0)
  ,fDXInv(0)
  ,fCacheValid(kFALSE)
  ,fCacheY2X(0),fCacheZ2X(0),fCacheSector(0),fCacheIXLow(0)
{
// def. c-tor  
  for (int i=4;i--;) fCacheDistLow[i]=fCacheDistUp[i] = 0;
}

//____________________________________________________________________
AliTPCChebDist::AliTPCChebDist(const char* name, const char* title, 
			       int nps, int nzs, float zmaxAbs)
  :AliTPCChebCorr(name,title,nps,nzs,zmaxAbs,0,0)
  ,fXMin(fgRMinTPC)
  ,fXMax(fgRMaxTPC)
  ,fDX(0)
  ,fDXInv(0)
  ,fCacheValid(kFALSE)
  ,fCacheY2X(0),fCacheZ2X(0),fCacheSector(0),fCacheIXLow(0)
{
  // c-tor
  fNRows = fgNSlices;
  if (fNRows<2) AliFatalF("Number of rows=%d cannot be<2",fNRows);
  fDX = (fXMax-fXMin)/(fNRows-1);
  if (fDX<=0)  AliFatalF("X boundaries are not increasing: %f %f",fXMin,fXMax);
  fDXInv = 1./fDX;
  for (int i=4;i--;) fCacheDistLow[i]=fCacheDistUp[i] = 0;
  //
}

//____________________________________________________________________
void AliTPCChebDist::Eval(int sector, float x, float y2x, float z, float *distortion) const
{
  // Calculate distortion for point with x,y,z sector corrdinates (sector id in 0-71 format)
  int ixLow = X2Slice(x);
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  const AliCheb2DStack* chpar = GetParam(sector,y2x,z);
  //
  const float kMaxY2XDiff=3.e-3; // y2x spans +-1.7632e-01, take ~1% precision
  const float kMaxZ2XDiff=1.e-2; // z2x spans 0:1, take ~1% precision  
  // check if cache is valid
  if (fCacheValid && 
      (sector!=fCacheSector || ixLow!=fCacheIXLow || 
       TMath::Abs(y2x-fCacheY2X)>kMaxY2XDiff  || 
       TMath::Abs(z-fCacheZ2X)>kMaxZ2XDiff))
    fCacheValid = kFALSE;
  //
  float scl = (x-Slice2X(ixLow))*fDXInv; // lever arm for interpolation
  if (!fCacheValid) chpar->Eval(ixLow  , tz, fCacheDistLow);
  if (ixLow<fNRows-1) {
    if (!fCacheValid) chpar->Eval(ixLow+1, tz, fCacheDistUp);
    for (int i=4;i--;) distortion[i] = fCacheDistLow[i]+scl*(fCacheDistUp[i]-fCacheDistLow[i]); // linear interpolation
  }
  else { // we are at the last slice, extrapolate
    if (!fCacheValid) chpar->Eval(ixLow-1, tz, fCacheDistUp);
    for (int i=4;i--;) distortion[i] = fCacheDistLow[i]+scl*(fCacheDistLow[i]-fCacheDistUp[i]); // linear extrapolation   
  }
  // protection against extrapolation of dispersion making it negative
  if (distortion[3]<1e-3) distortion[3] = 1e-3;
  //
  if (!fCacheValid) { // memorize
    fCacheValid = kTRUE;
    fCacheY2X = y2x;
    fCacheZ2X = z;
    fCacheSector = sector;
    fCacheIXLow = ixLow;
  }
  //
}

//____________________________________________________________________
void AliTPCChebDist::Init()
{
  // mark cach invalid and proceed to initialization
  fCacheValid = kFALSE;
  AliTPCChebCorr::Init();
}
