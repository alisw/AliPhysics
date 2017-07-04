/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Fast polynomial parametrization of Alice magnetic field, to be used for reconstruction.
// Solenoid part fitted by Shuto Yamasaki from AliMagWrapCheb in the |Z|<260Interface and R<500 cm 
// Dipole part: to do
//
// Author: ruben.shahoyan@cern.ch
//
#include "AliMagFast.h"
#include "AliLog.h"
#include <math.h>
#include <fstream>
#include <sstream>

const float AliMagFast::fgkSolR2Max[AliMagFast::kNSolRRanges] =
  {80.f*80.f,250.f*250.f,400.f*400.f,423.f*423.f, 500.f*500.f};

const float AliMagFast::fgkSolZMax = 260.0f;

ClassImp(AliMagFast)

AliMagFast::AliMagFast(const char* inpFName)
{
  // c-tor
  memset(fSolPar,0,sizeof(SolParam_t)*kNSolRRanges*kNSolZRanges*kNQuadrants);
  if (inpFName && !LoadData(inpFName)) {
    AliFatalF("Failed to initialize from %s",inpFName);
  }
}


Bool_t AliMagFast::LoadData(const char* inpFName)
{
  // load field from text file

  std::ifstream in(inpFName,std::ifstream::in);
  if (in.bad()) {
    AliFatalF("Failed to open input file %s",inpFName);
    return kFALSE;
  }
  std::string line;
  int valI,component = -1, nParams = 0, header[4] = {-1,-1,-1,-1}; // iR, iZ, iQuadrant, nVal
  SolParam_t* curParam = 0; //std::nullptr;
  
  while (std::getline(in, line)) {

    if (line.empty() || line[0]=='#') continue; // empy or comment
    std::stringstream ss(line);
    int cnt=0;
    
    if (component<0) {
      while (cnt<4 && (ss>>header[cnt++]));
      if (cnt!=4) {
	AliFatalF("Wrong header: %s",line.c_str());
	return false;
      }
      curParam = &fSolPar[header[0]][header[1]][header[2]];
    }
    else {
      while (cnt<header[3] && (ss>>curParam->mParBxyz[component][cnt++]));
      if (cnt!=header[3]) {
	AliFatalF("Wrong data (npar=%d) for param %d %d %d %d: %s",cnt,header[0],header[1],header[2],header[3],line.c_str());
	return false;	
      }
    }    
    component++;
    if (component>2) {
      component = -1; // next header expected
      nParams++;
    }
  }
  //
  AliInfoF("loaded %d params from %s",nParams,inpFName);
  return kTRUE;
}
  
Bool_t AliMagFast::Field(const double xyz[3], double bxyz[3]) const
{
  float x(xyz[0]),y(xyz[1]),z(xyz[2]);

  // Z segment
  int zSeg = -1;
  if (z<fgkSolZMax) {
    if (zSeg>-fgkSolZMax) zSeg = z<0.f ? 0:1; // solenoid params
    else { // need to check dipole params
      return kFALSE;
    }
  }
  // R segment
  float xx = x*x, yy = y*y, rr = xx+yy;
  int rSeg;

  for (rSeg=0;rSeg<kNSolRRanges;rSeg++) if (rr < fgkSolR2Max[rSeg]) break;
  if (rSeg==kNSolRRanges) { // redefine to max allowed radius
    float scl2Bond = sqrtf(fgkSolR2Max[kNSolRRanges-1]/rr);
    x *= scl2Bond;
    y *= scl2Bond;
  }

  // quadrant
  int quadrant = GetQuadrant(x,y);
  const SolParam_t *par = &fSolPar[rSeg][zSeg][quadrant];
  bxyz[kX] = CalcPol(par->mParBxyz[kX],x,y,z);
  bxyz[kY] = CalcPol(par->mParBxyz[kY],x,y,z);
  bxyz[kZ] = CalcPol(par->mParBxyz[kZ],x,y,z);
  //
  return kTRUE;
}

Bool_t AliMagFast::GetBz(const double xyz[3], double& bz) const
{
  float x(xyz[0]),y(xyz[1]),z(xyz[2]);

  // Z segment
  int zSeg = -1;
  if (z<fgkSolZMax) {
    if (zSeg>-fgkSolZMax) zSeg = z<0.f ? 0:1; // solenoid params
    else { // need to check dipole params
      return kFALSE;
    }
  }
  // R segment
  float xx = x*x, yy = y*y, rr = xx+yy;
  int rSeg;

  for (rSeg=0;rSeg<kNSolRRanges;rSeg++) if (rr < fgkSolR2Max[rSeg]) break;
  if (rSeg==kNSolRRanges) { // redefine to max allowed radius
    float scl2Bond = sqrtf(fgkSolR2Max[kNSolRRanges-1]/rr);
    x *= scl2Bond;
    y *= scl2Bond;
  }

  // quadrant
  int quadrant = GetQuadrant(x,y);
  const SolParam_t *par = &fSolPar[rSeg][zSeg][quadrant];
  bz = CalcPol(par->mParBxyz[kZ],x,y,z);
  //
  return kTRUE;
}
