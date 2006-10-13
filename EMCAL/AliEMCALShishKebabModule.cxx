/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//_________________________________________________________________________
// Main class for "twist" geometry of Shish-Kebab case.
// Author: Aleksei Pavlinov(WSU).
// Sep 20004.
// See web page with description of Shish-Kebab geometries:
// http://pdsfweb01.nersc.gov/~pavlinov/ALICE/SHISHKEBAB/RES/shishkebabALICE.html
//_________________________________________________________________________

#include "AliEMCALShishKebabModule.h"
#include "AliEMCALGeometry.h"
#include <TGraph.h>
#include <TMath.h>

ClassImp(AliEMCALShishKebabModule)

  AliEMCALGeometry *AliEMCALShishKebabModule::fgGeometry=0; 
  Double_t AliEMCALShishKebabModule::fga=0.; 
  Double_t AliEMCALShishKebabModule::fgb=0.; 
  Double_t AliEMCALShishKebabModule::fgr=0.; 

//_________________________________________________________________________
AliEMCALShishKebabModule::AliEMCALShishKebabModule() 
  : TNamed(),
    fOK(0),
    fA(0.),
    fB(0.),
    fTheta(0.)
{ 
  // theta in radians ; first object shold be with theta=pi/2.
  if(fgGeometry==0) {
    fTheta = TMath::PiOver2();
    if(GetParameters()) {
      DefineFirstModule();
      DefineName(fTheta);
    }
  } else {
    Warning("AliEMCALShishKebabModule(theta)","You should call this constractor just once !!");
  }
}

//_________________________________________________________________________
AliEMCALShishKebabModule::AliEMCALShishKebabModule(AliEMCALShishKebabModule &leftNeighbor) 
  : TNamed(),
    fOK(0),
    fA(0.),
    fB(0.),
    fTheta(0.)
{ 
  // 22-sep-04
  TObject::SetUniqueID(leftNeighbor.GetUniqueID()+1);
  Init(leftNeighbor.GetA(),leftNeighbor.GetB());
}

//_________________________________________________________________________
AliEMCALShishKebabModule::AliEMCALShishKebabModule(const AliEMCALShishKebabModule& mod) 
  : TNamed(mod.GetName(),mod.GetTitle()),
    fOK(mod.fOK),
    fA(mod.fA),
    fB(mod.fB),
    fTheta(mod.fTheta)
{
  //copy ctor
}

//_________________________________________________________________________
void AliEMCALShishKebabModule::Init(Double_t A, Double_t B)
{ 
  //
  // Initialisation method
  //
  Double_t thetaMin, thetaMax, par[4];
  Int_t npar=0;
  if(A<0){
    //    DefineSecondModuleFirstAssumption();
    thetaMax = TMath::ATan2(fgr, 1.4*fga);
    thetaMin = TMath::ATan2(fgr, 1.6*fga);
    fTheta   = Solve(AliEMCALShishKebabModule::Y2, thetaMin,thetaMax,npar,par);
  } else{
    npar = 4;
    par[0] = fga;
    par[1] = fgr;
    par[2] = A;
    par[3] = B;
    Double_t x = fgr/A;
    thetaMax = TMath::ATan2(fgr,x + 0.5*fga);
    thetaMin = TMath::ATan2(fgr,x + 1.5*fga);
    fTheta   = Solve(AliEMCALShishKebabModule::YALL, thetaMin,thetaMax,npar,par);
  }

  Double_t rOK = fgr / TMath::Sin(fTheta) + (fga/2.)/TMath::Tan(fTheta) + fgb/2.;
  fOK.SetMagPhi(rOK, fTheta); 
  // have to define A and B
  fA = TMath::Tan(fTheta);
  fB = -fga/(2.*TMath::Cos(fTheta));
  DefineName(fTheta);
}

//_________________________________________________________________________
void AliEMCALShishKebabModule::DefineFirstModule()
{
  // Define first module
  fOK.Set(fga/2., fgr + fgb/2.); // position the center of module vs o

  fB = fga/2.;    // z=fB
  fA = -999.;     // fA=infinite in this case
  TObject::SetUniqueID(1); //
}

//_________________________________________________________________________
void AliEMCALShishKebabModule::DefineSecondModuleFirstAssumption()
{ // Keep for testing and checking
  // cos(theta) << 1, theta ~ pi/2.; a/r = 11.4/462.54 = 0.0246465 << 1; 
  // theta=1.53382  from this case; theta=1.533869 from TGraph::Zero 
  Double_t x = (3.*fga)/(2.*fgr);
  fTheta = TMath::ACos(x);
  /*
  Double_t rOK = fgr / TMath::Sin(fTheta) + (fga/2.)/TMath::Tan(fTheta) + fgb/2.;
  fOK.SetMagPhi(rOK, fTheta); 
  // have to define A and B
  fA = TMath::Tan(fTheta);
  fB = -fga/(2.*TMath::Cos(fTheta));
  DefineName(fTheta);
  */
}

//_________________________________________________________________________
Double_t AliEMCALShishKebabModule::Solve(Double_t (*fcn)(Double_t*,Double_t*), 
Double_t xmin, Double_t xmax, Int_t npar, Double_t *par, Double_t eps, Int_t maxIter)
{
  // Find out "zero" using TGraph method
  if(npar); // unused now
  TGraph gr;
  Double_t x,y;
  Int_t k = 0;
  gr.Zero(k, xmin,xmax, eps, x,y, maxIter); // remember initial interval
  while(k!=2) {
    y = fcn(&x, par); 
    gr.Zero(k, xmin,xmax, eps, x,y, maxIter);
  }
  return x;
}

//_________________________________________________________________________
Double_t AliEMCALShishKebabModule::Y2(Double_t *x, Double_t *par)
{ 
  // For position calulation of second module
  if(par);
  Double_t theta = x[0];
  Double_t cos = TMath::Cos(theta);
  Double_t sin = TMath::Sin(theta);
  Double_t y1  = fgr*cos/sin + fga/(2.*sin) - fga*sin;       
  Double_t y2  = fga, y = y1-y2;;
  //  printf(" theta %f Y %12.5e \n", theta, y);
  return y;
}

//_________________________________________________________________________
Double_t AliEMCALShishKebabModule::YALL(Double_t *x, Double_t *par)
{ 
  // For position calulation of 3th, 4th to 30th modules
  Double_t a=par[0], r=par[1], aa=par[2], bb=par[3]; 
  Double_t theta = x[0];
  Double_t cos = TMath::Cos(theta);
  Double_t sin = TMath::Sin(theta);

  Double_t y1  = r + a*cos;       
  Double_t y2  = aa*(r*cos/sin + a/(2.*sin) - a*sin) + bb;
  Double_t y   = y1-y2;
  //  printf(" theta %f Y %12.5e \n", theta, y);
  return y;
}

//_________________________________________________________________________
void AliEMCALShishKebabModule::DefineName(Double_t theta)
{
  // Define name of object
  SetName(Form("%2i(%5.2f)", TObject::GetUniqueID(), theta*TMath::RadToDeg()));
}

//_________________________________________________________________________
Bool_t AliEMCALShishKebabModule::GetParameters()
{
  // Get needing module parameters from EMCAL geometry
  fgGeometry = AliEMCALGeometry::GetInstance();
  if(!fgGeometry) {
    Warning("GetParameters()"," No geometry ");
    return kFALSE; 
  }

  fga = (Double_t)fgGeometry->GetPhiModuleSize();
  fgb = (Double_t)fgGeometry->GetLongModuleSize();
  fgr = (Double_t)(fgGeometry->GetIPDistance() + fgGeometry->GetSteelFrontThickness());
  PrintShish(0);
  return kTRUE;
}

//_________________________________________________________________________
void AliEMCALShishKebabModule::PrintShish(Int_t pri) const
{
  // service method
  if(pri>=0) {
    Info("PrintShish()", " a %7.2f | b %7.2f | r %7.2f ", fga, fgb, fgr);
    printf(" fTheta %f : %5.2f : cos(theta) %f\n", fTheta, GetThetaInDegree(),TMath::Cos(fTheta)); 
    if(pri>0) {
      printf("%i %s | theta %f -> %f\n", GetUniqueID(), GetName(), fTheta, fOK.Phi());
      printf(" A %f B %f \n", fA, fB);

      fOK.Dump();
    }
  }
}

//_________________________________________________________________________
Double_t AliEMCALShishKebabModule::GetThetaInDegree() const 
{
  return fTheta*TMath::RadToDeg();
}
