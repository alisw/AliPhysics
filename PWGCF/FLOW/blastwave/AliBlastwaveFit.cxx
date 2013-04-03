#include<stdio.h>

#include"AliBlastwaveFit.h"

#include"TMath.h"

ClassImp(AliBlastwaveFit);

AliBlastwaveFit::AliBlastwaveFit(const char *name,Double_t mass) :
  TNamed(name,name),
  fMass(mass),
  fFunctionYield(NULL),
  fFunctionV2(NULL),
  fSpectraObj(NULL),
  fV2Obj(NULL),
  fSpectraObjCopy(NULL),
  fXmin(0.2),
  fXmax(3.0)
{
}
//------------------------------------------------------------------------------
AliBlastwaveFit::AliBlastwaveFit() :
  TNamed("BlastwaveFit","BlastwaveFit"),
  fMass(0.0),
  fFunctionYield(NULL),
  fFunctionV2(NULL),
  fSpectraObj(NULL),
  fV2Obj(NULL),
  fSpectraObjCopy(NULL),
  fXmin(0.2),
  fXmax(3.0)
{  
}
//------------------------------------------------------------------------------
AliBlastwaveFit::~AliBlastwaveFit(){
}
