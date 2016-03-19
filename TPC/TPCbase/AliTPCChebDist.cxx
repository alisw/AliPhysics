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
Int_t   AliTPCChebDist::fgNSlices = AliTPCChebCorr::kNRows+10; // check

//____________________________________________________________________
AliTPCChebDist::AliTPCChebDist()
  : AliTPCChebCorr()
  ,fXMin(fgRMinTPC)
  ,fXMax(fgRMaxTPC)
  ,fDX(0)
  ,fDXInv(0)
{
// def. c-tor  
}

//____________________________________________________________________
AliTPCChebDist::AliTPCChebDist(const char* name, const char* title, 
			       int nps, int nzs, float zmaxAbs)
  :AliTPCChebCorr(name,title,nps,nzs,zmaxAbs)
  ,fXMin(fgRMinTPC)
  ,fXMax(fgRMaxTPC)
{
  // c-tor
  fNRows = fgNSlices;
  if (fNRows<2) AliFatalF("Number of rows=%d cannot be<2",fNRows);
  fDX = (fXMax-fXMin)/(fNRows-1);
  if (fDX<=0)  AliFatalF("X boundaries are not increasing: %f %f",fXMin,fXMax);
  fDXInv = 1./fDX;
  //
}

//____________________________________________________________________
void AliTPCChebDist::Eval(int sector, float x, float y2x, float z, float *distortion) const
{
  // Calculate distortion for point with x,y,z sector corrdinates (sector id in 0-71 format)
  int ixLow = X2Slice(x);
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  const AliCheb2DStack* chpar = GetParam(sector,y2x,z);
  float distUp[3], scl = (x-Slice2X(ixLow))*fDXInv; // lever arm for interpolation
  chpar->Eval(ixLow  , tz, distortion);
  if (ixLow<fNRows-1) {
    chpar->Eval(ixLow+1, tz, distUp);
    for (int i=3;i--;) distortion[i] += scl*(distUp[i]-distortion[i]); // linear interpolation
  }
  else { // we are at the last slice, extrapolate
    chpar->Eval(ixLow-1, tz, distUp);
    for (int i=3;i--;) distortion[i] += scl*(distortion[i]-distUp[i]); // linear extrapolation   
  }
  //
}
