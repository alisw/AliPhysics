/**************************************************************************
 * Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
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

//*-- Author: Aleksei Pavlinov(WSU)

#include "AliEMCALShishKebabTrd1Module.h"
#include <assert.h>
#include "AliEMCALGeometry.h"

#include "Riostream.h"
#include <TMath.h>

ClassImp(AliEMCALShishKebabTrd1Module)

  AliEMCALGeometry *AliEMCALShishKebabTrd1Module::fgGeometry=0; 
  Double_t AliEMCALShishKebabTrd1Module::fga=0.; 
  Double_t AliEMCALShishKebabTrd1Module::fga2=0.; 
  Double_t AliEMCALShishKebabTrd1Module::fgb=0.; 
  Double_t AliEMCALShishKebabTrd1Module::fgr=0.; 
  Double_t AliEMCALShishKebabTrd1Module::fgangle=0.;   // around one degree 
  Double_t AliEMCALShishKebabTrd1Module::fgtanBetta=0; //

AliEMCALShishKebabTrd1Module::AliEMCALShishKebabTrd1Module(double theta, AliEMCALGeometry *g) : TNamed()
{ // theta in radians ; first object shold be with theta=pi/2.
  cout<< " theta " << theta << " geometry " << g << endl;  
  fTheta = theta;
  if(fgGeometry==0) {
    fgGeometry = g;
    if(GetParameters()) {
      DefineFirstModule();
    }
  } else Warning("AliEMCALShishKebabTrd1Module(theta)","You should call this constractor just once !!");
  DefineName(fTheta);
}

AliEMCALShishKebabTrd1Module::AliEMCALShishKebabTrd1Module(AliEMCALShishKebabTrd1Module &leftNeighbor) : TNamed()
{ // 22-sep-04
  //  printf("** Left Neighbor : %s **\n", leftNeighbor.GetName());
  TObject::SetUniqueID(leftNeighbor.GetUniqueID()+1);
  fTheta  = leftNeighbor.GetTheta() - fgangle; 
  Init(leftNeighbor.GetA(),leftNeighbor.GetB());
}

void AliEMCALShishKebabTrd1Module::Init(double A, double B)
{ // Define parameter module from parameters A,B from previos.
  Double_t yl = (fgb/2)*TMath::Sin(fTheta) + (fga/2)*TMath::Cos(fTheta) + fgr, y = yl;
  Double_t xl = (yl - B) / A;     // y=A*x+B

  //  Double_t xp1 = (fga/2. + fgb/2.*fgtanBetta)/(TMath::Sin(fTheta) + fgtanBetta*TMath::Cos(fTheta));
  //  printf(" xp1 %9.3f \n ", xp1);
  // xp1 == xp => both methods give the same results - 3-feb-05
  Double_t alpha = TMath::Pi()/2. + fgangle/2;
  Double_t xt = (fga+fga2)*TMath::Tan(fTheta)*TMath::Tan(alpha)/(4.*(1.-TMath::Tan(fTheta)*TMath::Tan(alpha)));
  Double_t yt = xt / TMath::Tan(fTheta), xp = TMath::Sqrt(xt*xt + yt*yt);
  Double_t x  = xl + xp;
  fOK.Set(x, y);
  //  printf(" yl %9.3f | xl %9.3f | xp %9.3f \n", yl, xl, xp);

  // have to define A and B; 
  Double_t yCprev = fgr + fga*TMath::Cos(fTheta);
  Double_t xCprev = (yCprev - B) / A;
  Double_t xA     = xCprev + fga*TMath::Sin(fTheta), yA = fgr;

  fThetaA = fTheta - fgangle/2.;
  fA = TMath::Tan(fThetaA); // !!
  fB = yA - fA*xA;

  DefineName(fTheta);
  // Centers of module
  Double_t kk1 = (fga+fga2)/(2.*4.); // kk1=kk2 

  Double_t xk1 = fOK.X() - kk1*TMath::Sin(fTheta);
  Double_t yk1 = fOK.Y() + kk1*TMath::Cos(fTheta) - fgr;
  fOK1.Set(xk1,yk1);

  Double_t xk2 = fOK.X() + kk1*TMath::Sin(fTheta);
  Double_t yk2 = fOK.Y() - kk1*TMath::Cos(fTheta) - fgr;
  fOK2.Set(xk2,yk2);
}

void AliEMCALShishKebabTrd1Module::DefineFirstModule()
{
  fOK.Set(fga2/2., fgr + fgb/2.); // position the center of module vs o

  // parameters of right line : y = A*z + B in system where zero point is IP.
  fThetaA = fTheta - fgangle/2.;
  fA      = TMath::Tan(fThetaA);
  Double_t xA = fga/2. + fga2/2., yA = fgr;
  fB = yA - fA*xA;

  Double_t kk1 = (fga+fga2)/(2.*4.); // kk1=kk2 
  fOK1.Set(fOK.X() - kk1, fOK.Y()-fgr);
  fOK2.Set(fOK.X() + kk1, fOK.Y()-fgr);

  TObject::SetUniqueID(1); //
}

void AliEMCALShishKebabTrd1Module::DefineName(double theta)
{
  char name[100];
  // sprintf(name,"theta_%5.2f",theta*180./TMath::Pi());
  sprintf(name,"%2i(%5.2f)", TObject::GetUniqueID(), theta*180./TMath::Pi());
  SetName(name);
}

Bool_t AliEMCALShishKebabTrd1Module::GetParameters()
{
  if(!fgGeometry) fgGeometry = AliEMCALGeometry::GetInstance();
  TString sn(fgGeometry->GetName()); // 2-Feb-05
  sn.ToUpper();
  if(!fgGeometry) {
    Warning("GetParameters()"," No geometry ");
    return kFALSE; 
  }

  fga        = (Double_t)fgGeometry->GetEtaModuleSize();
  fgb        = (Double_t)fgGeometry->GetLongModuleSize();
  fgangle    = Double_t(fgGeometry->GetTrd1Angle())*TMath::DegToRad();
  fgtanBetta = TMath::Tan(fgangle/2.);
  fgr        = (Double_t)fgGeometry->GetIPDistance();
  if(!sn.Contains("TRD2")) fgr += fgGeometry->GetSteelFrontThickness();
  fga2       = Double_t(fgGeometry->Get2Trd1Dx2());
  PrintShish(0);
  return kTRUE;
}

// service methods
void AliEMCALShishKebabTrd1Module::PrintShish(int pri) const
{
  if(pri>=0) {
    Info("\n PrintShish()", "\n a %7.3f:%7.3f | b %7.2f | r %7.2f \n TRD1 angle %7.6f(%5.2f) | tanBetta %7.6f", 
    fga, fga2, fgb, fgr, fgangle, fgangle*TMath::RadToDeg(), fgtanBetta);
    printf(" fTheta %f : %5.2f : cos(theta) %f\n", 
    fTheta, GetThetaInDegree(),TMath::Cos(fTheta)); 
    if(pri>=1) {
      printf("\n%i |%s| theta %f :  fOK.Phi = %f(%5.2f)\n", 
      GetUniqueID(), GetName(), fTheta, fOK.Phi(),fOK.Phi()*TMath::RadToDeg());
      printf(" A %f B %f | fThetaA %7.6f(%5.2f)\n", fA,fB, fThetaA,fThetaA*TMath::RadToDeg());
      printf(" fOK  : X %8.3f: Y %8.3f \n",  fOK.X(), fOK.Y());
      printf(" fOK1 : X %8.3f: Y %8.3f (in SM, ieta=2)\n", fOK1.X(), fOK1.Y());
      printf(" fOK2 : X %8.3f: Y %8.3f (in SM, ieta=1)\n\n", fOK2.X(), fOK2.Y());
      fOK.Dump();
    }
  }
}
