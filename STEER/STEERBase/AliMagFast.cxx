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
#include <TString.h>
#include <TSystem.h>
#include <math.h>
#include <fstream>
#include <sstream>

const float AliMagFast::fgkSolR2Max[AliMagFast::kNSolRRanges] =
  {80.f*80.f,250.f*250.f,400.f*400.f,423.f*423.f, 500.f*500.f};

const float AliMagFast::fgkSolZMax = 550.0f;

ClassImp(AliMagFast)

AliMagFast::AliMagFast(const char* inpFName) :
fFactorSol(1.f)
{
  // c-tor
  memset(fSolPar,0,sizeof(SolParam_t)*kNSolRRanges*kNSolZRanges*kNQuadrants);
  if (inpFName && !LoadData(inpFName)) {
    AliFatalF("Failed to initialize from %s",inpFName);
  }
}

AliMagFast::AliMagFast(Float_t factor, Int_t nomField, const char* inpFmt) :
fFactorSol(factor)
{
  // c-tor
  if (nomField!=2 && nomField!=5) {
    AliFatalF("No parametrization for nominal field of %d kG",nomField);
  }
  TString pth = Form(inpFmt,nomField);
  memset(fSolPar,0,sizeof(SolParam_t)*kNSolRRanges*kNSolZRanges*kNQuadrants);
  if (!LoadData(pth.Data())) {
    AliFatalF("Failed to initialize from %s",pth.Data());
  }
}

//_______________________________________________________________________
AliMagFast::AliMagFast(const AliMagFast &src):
  fFactorSol(src.fFactorSol)
{
  memcpy(fSolPar,src.fSolPar, kNSolRRanges*kNSolZRanges*kNQuadrants*sizeof(SolParam_t));
}

AliMagFast& AliMagFast::operator=(const AliMagFast& src)
{
  if (this != &src) {
    fFactorSol = src.fFactorSol;
    memcpy(fSolPar,src.fSolPar, kNSolRRanges*kNSolZRanges*kNQuadrants*sizeof(SolParam_t));
  }
  return *this;
}


Bool_t AliMagFast::LoadData(const char* inpFName)
{
  // load field from text file

  std::ifstream in(gSystem->ExpandPathName(inpFName),std::ifstream::in);
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
  // get field
  const float fxyz[3]={float(xyz[0]),float(xyz[1]),float(xyz[2])};
  int zSeg,rSeg,quadrant;
  if (!GetSegment(fxyz,zSeg,rSeg,quadrant)) return kFALSE;
  const SolParam_t *par = &fSolPar[rSeg][zSeg][quadrant];
  bxyz[kX] = CalcPol(par->mParBxyz[kX],fxyz[kX],fxyz[kY],fxyz[kZ])*fFactorSol;
  bxyz[kY] = CalcPol(par->mParBxyz[kY],fxyz[kX],fxyz[kY],fxyz[kZ])*fFactorSol;
  bxyz[kZ] = CalcPol(par->mParBxyz[kZ],fxyz[kX],fxyz[kY],fxyz[kZ])*fFactorSol;
  //
  return kTRUE;
}

Bool_t AliMagFast::GetBz(const double xyz[3], double& bz) const
{
  // get field
  const float fxyz[3]={float(xyz[0]),float(xyz[1]),float(xyz[2])};
  int zSeg,rSeg,quadrant;
  if (!GetSegment(fxyz,zSeg,rSeg,quadrant)) return kFALSE;
  const SolParam_t *par = &fSolPar[rSeg][zSeg][quadrant];
  bz = CalcPol(par->mParBxyz[kZ],fxyz[kX],fxyz[kY],fxyz[kZ])*fFactorSol;
  //
  return kTRUE;
}
  
Bool_t AliMagFast::Field(const float xyz[3], float bxyz[3]) const
{
  // get field
  int zSeg,rSeg,quadrant;
  if (!GetSegment(xyz,zSeg,rSeg,quadrant)) return kFALSE;
  const SolParam_t *par = &fSolPar[rSeg][zSeg][quadrant];
  bxyz[kX] = CalcPol(par->mParBxyz[kX],xyz[kX],xyz[kY],xyz[kZ])*fFactorSol;
  bxyz[kY] = CalcPol(par->mParBxyz[kY],xyz[kX],xyz[kY],xyz[kZ])*fFactorSol;
  bxyz[kZ] = CalcPol(par->mParBxyz[kZ],xyz[kX],xyz[kY],xyz[kZ])*fFactorSol;
  //
  return kTRUE;
}

Bool_t AliMagFast::GetBz(const float xyz[3], float& bz) const
{
  // get field
  int zSeg,rSeg,quadrant;
  if (!GetSegment(xyz,zSeg,rSeg,quadrant)) return kFALSE;
  const SolParam_t *par = &fSolPar[rSeg][zSeg][quadrant];
  bz = CalcPol(par->mParBxyz[kZ],xyz[kX],xyz[kY],xyz[kZ])*fFactorSol;
  //
  return kTRUE;
}

Bool_t AliMagFast::GetSegment(const float xyz[3], int& zSeg,int &rSeg, int &quadrant) const
{
  // get segment of point location
  const float &x = xyz[kX], &y = xyz[kY], &z = xyz[kZ];
  const float zGridSpaceInv = 1/(fgkSolZMax*2/kNSolZRanges);
  zSeg = -1;
  if (z<fgkSolZMax) {
    if (z>-fgkSolZMax) zSeg = (z+fgkSolZMax)*zGridSpaceInv; // solenoid params
    else { // need to check dipole params
      return kFALSE;
    }
  }
  else return kFALSE;
  // R segment
  float xx = x*x, yy = y*y, rr = xx+yy;
  for (rSeg=0;rSeg<kNSolRRanges;rSeg++) if (rr < fgkSolR2Max[rSeg]) break;
  if (rSeg==kNSolRRanges) return kFALSE;
  quadrant = GetQuadrant(x,y);
  return kTRUE;
}
