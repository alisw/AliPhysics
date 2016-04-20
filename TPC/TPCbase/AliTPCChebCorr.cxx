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

const char* AliTPCChebCorr::fgkFieldTypeName[4] = {"Any","B>0"," B<0","B=0"};

const float AliTPCChebCorr::fgkY2XHSpan = TMath::Tan(TMath::Pi()/18);

const float AliTPCChebCorr::fgkPadRowX[AliTPCChebCorr::kNRows] = {
  85.225, 85.975, 86.725, 87.475, 88.225, 88.975, 89.725, 90.475, 91.225, 91.975, 92.725, 93.475, 94.225, 94.975, 95.725,
  96.475, 97.225, 97.975, 98.725, 99.475,100.225,100.975,101.725,102.475,103.225,103.975,104.725,105.475,106.225,106.975,
  107.725,108.475,109.225,109.975,110.725,111.475,112.225,112.975,113.725,114.475,115.225,115.975,116.725,117.475,118.225,
  118.975,119.725,120.475,121.225,121.975,122.725,123.475,124.225,124.975,125.725,126.475,127.225,127.975,128.725,129.475,
  130.225,130.975,131.725,135.100,136.100,137.100,138.100,139.100,140.100,141.100,142.100,143.100,144.100,145.100,146.100,
  147.100,148.100,149.100,150.100,151.100,152.100,153.100,154.100,155.100,156.100,157.100,158.100,159.100,160.100,161.100,
  162.100,163.100,164.100,165.100,166.100,167.100,168.100,169.100,170.100,171.100,172.100,173.100,174.100,175.100,176.100,
  177.100,178.100,179.100,180.100,181.100,182.100,183.100,184.100,185.100,186.100,187.100,188.100,189.100,190.100,191.100,
  192.100,193.100,194.100,195.100,196.100,197.100,198.100,199.350,200.850,202.350,203.850,205.350,206.850,208.350,209.850,
  211.350,212.850,214.350,215.850,217.350,218.850,220.350,221.850,223.350,224.850,226.350,227.850,229.350,230.850,232.350,
  233.850,235.350,236.850,238.350,239.850,241.350,242.850,244.350,245.850
};

//____________________________________________________________________
AliTPCChebCorr::AliTPCChebCorr()
  : TNamed()
  ,fFieldType(kFieldAny)
  ,fNRows(0)
  ,fNStacksSect(0)
  ,fNStacksZSect(0)
  ,fNStacksZ(0)
  ,fNStacks(0)
  ,fZMaxAbs(-1)
  ,fTimeStampStart(0)
  ,fTimeStampEnd(0xffffffff)
  ,fZScaleI(0)
  ,fY2XScaleI(0)
  ,fDeadZone(0)
  ,fRowXI(0)
  ,fParams(0)
{
// def. c-tor  
}

//____________________________________________________________________
AliTPCChebCorr::AliTPCChebCorr(const char* name, const char* title, 
			       int nps, int nzs, float zmaxAbs, float deadZone, const float *xrow)
  : TNamed(name,title)
  ,fFieldType(kFieldAny)
  ,fNRows(kNRows)
  ,fNStacksSect(0)
  ,fNStacksZSect(0)
  ,fNStacksZ(0)
  ,fNStacks(0)
  ,fZMaxAbs(-1)
  ,fTimeStampStart(0)
  ,fTimeStampEnd(0xffffffff)
  ,fZScaleI(0)
  ,fY2XScaleI(0)
  ,fDeadZone(deadZone)
  ,fRowXI(0)
  ,fParams(0)
{
  // c-tor, optional parameter xrow provides the X of every slice if deadZone is requested 
  // (if xrow is absent, default will be used)
  //
  if (deadZone>0) {
    fRowXI = new Float_t[kNRows];
    if (xrow) for (int i=kNRows;i--;) fRowXI[i] = 1./xrow[i];
    else      for (int i=kNRows;i--;) fRowXI[i] = 1.0f/fgkPadRowX[i];
  }
  //
  SetBinning(nps,nzs,zmaxAbs);
  //
}

//____________________________________________________________________
AliTPCChebCorr::~AliTPCChebCorr()
{
  // d-tor
  if (fParams) for (int i=fNStacks;i--;) delete fParams[i];
  delete[] fRowXI;
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
  float y2xMax = GetMaxY2X(); // half-sector span
  float bmn[2],bmx[2];
  Bool_t useS = GetUseShortPrec();  // float or short representation for param
  fParams = new AliCheb2DStack*[fNStacks];
  //
  for (int iz=0;iz<fNStacksZ;iz++) {
    bmn[1] = iz/fZScaleI - fZMaxAbs;     // boundaries in Z
    bmx[1] = iz<fNStacksZ-1 ? (iz+1)/fZScaleI - fZMaxAbs : fZMaxAbs;
    for (int isc=0;isc<kNSectors;isc++) {
      int isc72 = isc;
      if (iz<fNStacksZSect) isc72 += 18; // change side
      fun(isc72,0,0);
      for (int isl=0;isl<fNStacksSect;isl++) {
	float dead[2] = {0.,0.};
	float *deadP = 0;
	float *xinv = 0;
	if (fDeadZone>0 && isl==0) { // assign dead zone
	  dead[0] = fDeadZone;
	  deadP = dead;
	  xinv = fRowXI;
	}
	if (fDeadZone>0 && isl==fNStacksSect-1) { // assign dead zone
	  dead[1] = fDeadZone;
	  deadP = dead;
	  xinv = fRowXI;
	}
	//
	bmn[0] = -y2xMax+isl/fY2XScaleI; // boundaries in phi
	bmx[0] = isl<fNStacksSect-1 ?  -y2xMax+(isl+1)/fY2XScaleI : y2xMax;
	int id = GetParID(iz,isc,isl);
	AliInfoF("Doing param #%03d Iz:%d Sect:%02d Slice:%d | %+.1f<Z<%+.1f %+.3f<y2x<%+.3f",
		 id,iz,isc,isl,bmn[1],bmx[1],bmn[0],bmx[0]);
	if (useS) fParams[id] = new  AliCheb2DStackS(fun,fNRows,dimOut,bmn,bmx,np,dead,xinv,prec);
	else      fParams[id] = new  AliCheb2DStackF(fun,fNRows,dimOut,bmn,bmx,np,dead,xinv,prec);
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
  float y2xMax =  GetMaxY2X(); // half-sector span
  //
  fParams = new AliCheb2DStack*[fNStacks];
  Bool_t useS = GetUseShortPrec();  // float or short representation for param
  //
  for (int iz=0;iz<fNStacksZ;iz++) {
    bmn[1] = iz/fZScaleI-fZMaxAbs;     // boundaries in Z
    bmx[1] = iz<fNStacksZ-1 ? (iz+1)/fZScaleI-fZMaxAbs : fZMaxAbs;
    for (int isc=0;isc<kNSectors;isc++) {
      int isc72 = isc;
      if (iz<fNStacksZSect) isc72 += 18; // change side
      fun(isc72,0,0);  // this will set sector ingo to function
      //
      for (int isl=0;isl<fNStacksSect;isl++) {
	float dead[2] = {0.,0.};
	float *deadP = 0;
	float *xinv = 0;
	if (fDeadZone>0 && isl==0) { // assign dead zone
	  dead[0] = fDeadZone;
	  deadP = dead;
	  xinv = fRowXI;
	}
	if (fDeadZone>0 && isl==fNStacksSect-1) { // assign dead zone
	  dead[1] = fDeadZone;
	  deadP = dead;
	  xinv = fRowXI;
	}
	//
	bmn[0] = -y2xMax+isl/fY2XScaleI; // boundaries in phi
	bmx[0] = isl<fNStacksSect-1 ?  -y2xMax+(isl+1)/fY2XScaleI : y2xMax;
	int id = GetParID(iz,isc,isl);
	AliInfoF("Doing param #%03d Iz:%d Sect:%02d Slice:%d | %+.1f<Z<%+.1f %+.3f<y2x<%+.3f",
		 id,iz,isc,isl,bmn[1],bmx[1],bmn[0],bmx[0]);
	if (useS) fParams[id] = new  AliCheb2DStackS(fun,fNRows,dimOut,bmn,bmx,np,dead,xinv,prec);
	else      fParams[id] = new  AliCheb2DStackF(fun,fNRows,dimOut,bmn,bmx,np,dead,xinv,prec);
      }
    }
  }
  //
  SetBit(kParamDone);
  //
}

//____________________________________________________________________
void AliTPCChebCorr::Print(const Option_t* opt) const
{
  // print itself
  printf("%s:%s Cheb2D[%c] Param: %d slices in %+.1f<%s<%+.1f %d per sector. DeadZone: %.1fcm\n",
	 GetName(),GetTitle(),GetUseFloatPrec()?'F':'S',
	 fNStacksZ,-fZMaxAbs,GetUseZ2R() ? "Z/R":"Z",fZMaxAbs,fNStacksSect,fDeadZone);
  printf("Time span: %ld:%ld TimeDependent flag: %s Field type: %s\n",fTimeStampStart,fTimeStampEnd,
	 GetTimeDependent() ? "ON ":"OFF", fgkFieldTypeName[fFieldType]);
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

//____________________________________________________________________
void AliTPCChebCorr::Init()
{
  // make necessary initializations
  if (fRowXI) {
    for (int i=fNStacks;i--;) {
      if (fParams[i] && !fParams[i]->GetXRowInv()) fParams[i]->SetXRowInv(fRowXI);
    }
  }
}

//____________________________________________________________________
Int_t AliTPCChebCorr::GetDimOut() const
{
  // get number of output dimensions
  const AliCheb2DStack* par = GetParam(0);
  return par ? par->GetDimOut() : 0;
}
