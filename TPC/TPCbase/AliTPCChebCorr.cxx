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
#include "AliTPCChebCorr.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliTPCChebCorr)


const float AliTPCChebCorr::fgkAngGap = 5e-4; // allow some empty ang. space between sectors
const float AliTPCChebCorr::fgkY2XHSpan = TMath::Tan(TMath::Pi()/AliTPCChebCorr::kNSectors - AliTPCChebCorr::fgkAngGap);

//____________________________________________________________________
AliTPCChebCorr::AliTPCChebCorr()
  : TNamed()
  ,fNStacksSect(0)
  ,fNStacksZSect(0)
  ,fNStacksZ(0)
  ,fNStacks(0)
  ,fZMaxAbs(-1)
  ,fTimeStampStart(0)
  ,fTimeStampEnd(0xffffffff)
  ,fZScaleI(0)
  ,fY2XScaleI(0)
  ,fParams(0)
{
// def. c-tor  
}

//____________________________________________________________________
AliTPCChebCorr::AliTPCChebCorr(const char* name, const char* title, 
			       int nps, int nzs, float zmaxAbs)
  : TNamed(name,title)
  ,fNStacksSect(0)
  ,fNStacksZSect(0)
  ,fNStacksZ(0)
  ,fNStacks(0)
  ,fZMaxAbs(-1)
  ,fTimeStampStart(0)
  ,fTimeStampEnd(0xffffffff)
  ,fZScaleI(0)
  ,fY2XScaleI(0)
  ,fParams(0)
{
  // c-tor
  SetBinning(nps,nzs,zmaxAbs);
  //
}

//____________________________________________________________________
AliTPCChebCorr::~AliTPCChebCorr()
{
  // d-tor
  if (fParams) for (int i=fNStacks;i--;) delete fParams[i];
  delete[] fParams;
}

//____________________________________________________________________
void AliTPCChebCorr::Parameterize(stFun_t fun,int dimOut,const int np[2],const float* prec)
{
  // build parameterizations for 2->dimout, on the same grid of np[0]xnp[1] points for
  // every output dimension, optionally with prec[i] precision for i-th dimension
  //
  if (TestBit(kParamDone)) {
    AliError("Parameterization is already done");
    return;
  }
  //
  if (fZMaxAbs<0) AliFatal("First the binning and Z limits should be set");
  //
  float bmn[2],bmx[2];
  Bool_t useS = GetUseShortPrec();  // float or short representation for param
  fParams = new AliCheb2DStack*[fNStacks];
  //
  for (int iz=0;iz<fNStacksZ;iz++) {
    bmn[1] = iz/fZScaleI - fZMaxAbs;     // boundaries in Z
    bmx[1] = iz<fNStacksZ-1 ? (iz+1)/fZScaleI - fZMaxAbs : fZMaxAbs;
    for (int isc=0;isc<kNSectors;isc++) {
      fun(isc,0,0);
      for (int isl=0;isl<fNStacksSect;isl++) {
	bmn[0] = -fgkY2XHSpan+isl/fY2XScaleI; // boundaries in phi
	bmx[0] = isl<fNStacksSect-1 ?  -fgkY2XHSpan+(isl+1)/fY2XScaleI : fgkY2XHSpan;
	int id = GetParID(iz,isc,isl);
	AliInfoF("Doing param #%03d Iz:%d Sect:%02d Slice:%d | %+.1f<Z<%+.1f %+.3f<y2x<%+.3f",
		 id,iz,isc,isl,bmn[1],bmx[1],bmn[0],bmx[0]);
	if (useS) fParams[id] = new  AliCheb2DStackS(fun,kNRows,dimOut,bmn,bmx,np,prec);
	else      fParams[id] = new  AliCheb2DStackF(fun,kNRows,dimOut,bmn,bmx,np,prec);
      }
    }
  }
  //
  SetBit(kParamDone);
  //
}

//____________________________________________________________________
void AliTPCChebCorr::Parameterize(stFun_t fun,int dimOut,const int np[][2],const float* prec)
{
  // build parameterizations for 2->dimout, on the same grid of np[0]xnp[1] points for
  // every output dimension, optionally with prec[i] precision for i-th dimension
  //
  if (TestBit(kParamDone)) {
    AliError("Parameterization is already done");
    return;
  }
  //
  float bmn[2],bmx[2];
  float y2xMax = TMath::Tan(TMath::Pi()/kNSectors); // half-sector span
  //
  fParams = new AliCheb2DStack*[fNStacks];
  Bool_t useS = GetUseShortPrec();  // float or short representation for param
  //
  for (int iz=0;iz<fNStacksZ;iz++) {
    bmn[1] = iz/fZScaleI-fZMaxAbs;     // boundaries in Z
    bmx[1] = iz<fNStacksZ-1 ? (iz+1)/fZScaleI-fZMaxAbs : fZMaxAbs;
    for (int isc=0;isc<kNSectors;isc++) {
      fun(isc,0,0);
      for (int isl=0;isl<fNStacksSect;isl++) {
	bmn[0] = -y2xMax+isl/fY2XScaleI; // boundaries in phi
	bmx[0] = isl<fNStacksSect-1 ?  -fgkY2XHSpan+(isl+1)/fY2XScaleI : fgkY2XHSpan;
	int id = GetParID(iz,isc,isl);
	AliInfoF("Doing param #%03d Iz:%d Sect:%02d Slice:%d | %+.1f<Z<%+.1f %+.3f<y2x<%+.3f",
		 id,iz,isc,isl,bmn[1],bmx[1],bmn[0],bmx[0]);
	if (useS) fParams[id] = new  AliCheb2DStackS(fun,kNRows,dimOut,bmn,bmx,np,prec);
	else      fParams[id] = new  AliCheb2DStackF(fun,kNRows,dimOut,bmn,bmx,np,prec);
      }
    }
  }
  //
  SetBit(kParamDone);
  //
}

//____________________________________________________________________
void AliTPCChebCorr::Eval(int sector, int row, float y2x, float z, float *corr) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // If sector in 0-71 ROC convention, possible Z-outlying is checked
  int iz = (z+fZMaxAbs)*fZScaleI, side = sector/kNSectors;
  sector %= kNSectors;
  // correct for eventual Z calculated in wrong ROC
  if (side)  {if (iz>=fNStacksZSect) iz = fNStacksZSect-1;} // C side
  else       {if (iz<fNStacksZSect)  iz = fNStacksZSect;}   // A side
  if (iz<0) iz=0; else if (iz>=fNStacksZ) iz=fNStacksZ-1;
  int is = (y2x+fgkY2XHSpan)*fY2XScaleI;
  if (is<0) is=0; else if (is>=fNStacksSect) is=fNStacksSect-1;
  float tz[2] = {y2x,z}; // params use row, Y/X, Z
  GetParam(GetParID(iz,sector,is))->Eval(row,tz,corr);
  //
}

//____________________________________________________________________
void AliTPCChebCorr::Eval(int sector, int row, float tz[2], float *corr) const
{
  // Calculate correction for point with x,y,z sector corrdinates
  // If sector in 0-71 ROC convention, possible Z-outlying is checked
  int iz = (tz[1]+fZMaxAbs)*fZScaleI,  side = sector/kNSectors;
  sector %= kNSectors;
  // correct for eventual Z calculated in wrong ROC
  if (side)  {if (iz>=fNStacksZSect) iz = fNStacksZSect-1;} // C side
  else       {if (iz<fNStacksZSect)  iz = fNStacksZSect;}   // A side
  if (iz<0) iz=0; else if (iz>=fNStacksZ) iz=fNStacksZ-1;
  int is = (tz[0]+fgkY2XHSpan)*fY2XScaleI;
  if (is<0) is=0; else if (is>=fNStacksSect) is=fNStacksSect-1;
  GetParam(GetParID(iz,sector,is))->Eval(row,tz,corr);
  //
}

//____________________________________________________________________
void AliTPCChebCorr::Print(const Option_t* opt) const
{
  // print itself
  printf("%s:%s Cheb2D[%c] Param: %d slices in %+.1f<%s<%+.1f %d per sector\n",
	 GetName(),GetTitle(),GetUseFloatPrec()?'F':'S',
	 fNStacksZ,-fZMaxAbs,GetUseZ2R() ? "Z/R":"Z",fZMaxAbs,fNStacksSect);
  printf("Time span: %10d:%10d TimeDependent flag: %s\n",fTimeStampStart,fTimeStampEnd,
	 GetTimeDependent() ? "ON":"OFF");
  TString opts = opt; opts.ToLower();
  if (opts.Contains("p") && TestBit(kParamDone)) {
    for (int iz=0;iz<fNStacksZ;iz++) {
      for (int isc=0;isc<kNSectors;isc++) {
	for (int isl=0;isl<fNStacksSect;isl++) {
	  int id = GetParID(iz,isc,isl);
	  printf("Z%d Sector%02d Slice%d: ",iz,isc,isl);
	  GetParam(id)->Print(opt);
	}
      }
    }
  }
}

//____________________________________________________________________
void AliTPCChebCorr::SetBinning(int nps,int nzs, float zmxAbs)
{
  // set binning, limits
  fNStacksSect = nps;
  fNStacksZSect = nzs;
  fNStacksZ = nzs*2; // nzs is the number of bins per side!
  fZMaxAbs  = zmxAbs;
  fNStacks  = fNStacksZ*fNStacksSect*kNSectors;
  //  
  if (zmxAbs<1e-8 || nzs<1 || nps<1) AliFatalF("Wrong settings: |Zmax|:%f Nz:%d Nphi:%d",
					       zmxAbs,nzs,nps);
  fZScaleI   = nzs/zmxAbs;
  fY2XScaleI = nps/(2*TMath::Tan(TMath::Pi()/kNSectors));
  //
}

