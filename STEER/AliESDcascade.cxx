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

/* $Id$ */

//-------------------------------------------------------------------------
//               Implementation of the cascade vertex class
//              This is part of the Event Summary Data 
//              which contains the result of the reconstruction
//              and is the main set of classes for analaysis
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliESDcascade.h"

ClassImp(AliESDcascade)

AliESDcascade::AliESDcascade() : TObject() {
  //--------------------------------------------------------------------
  // Default constructor  (Xi-)
  //--------------------------------------------------------------------
  fPdgCode=kXiMinus;
  fEffMass=1.32131;
  fChi2=1.e+33;
  fPos[0]=fPos[1]=fPos[2]=0.;
  fPosCov[0]=fPosCov[1]=fPosCov[2]=fPosCov[3]=fPosCov[4]=fPosCov[5]=0.;
}

inline Double_t det(Double_t a00, Double_t a01, Double_t a10, Double_t a11){
  // determinant 2x2
  return a00*a11 - a01*a10;
}

inline Double_t det (Double_t a00,Double_t a01,Double_t a02,
                     Double_t a10,Double_t a11,Double_t a12,
                     Double_t a20,Double_t a21,Double_t a22) {
  // determinant 3x3
  return 
  a00*det(a11,a12,a21,a22)-a01*det(a10,a12,a20,a22)+a02*det(a10,a11,a20,a21);
}

Double_t AliESDcascade::ChangeMassHypothesis(Double_t &v0q, Int_t code) {
  //--------------------------------------------------------------------
  // This function changes the mass hypothesis for this cascade
  // and returns the "kinematical quality" of this hypothesis
  // together with the "quality" of associated V0 (argument v0q) 
  //--------------------------------------------------------------------
  Double_t nmass=0.13957, pmass=0.93827, ps0=0.101; 
  Double_t bmass=0.13957, mass =1.3213,  ps =0.139;

  fPdgCode=code;

  switch (code) {
  case 213: 
       bmass=0.93827; 
       break;
  case kXiMinus:
       break;
  case kXiPlusBar:
       nmass=0.93827; pmass=0.13957; 
       break;
  case kOmegaMinus: 
       bmass=0.49368; mass=1.67245; ps=0.211;
       break;
  case kOmegaPlusBar: 
       nmass=0.93827; pmass=0.13957; 
       bmass=0.49368; mass=1.67245; ps=0.211;
       break;
  default:
       Info("AliCascadeVertex::ChangeMassHypothesis",
            "Invalide PDG code !  Assuming XiMinus's...");
       fPdgCode=kXiMinus;
    break;
  }

  Double_t pxn=fV0mom[0][0], pyn=fV0mom[0][1], pzn=fV0mom[0][2];
  Double_t pxp=fV0mom[1][0], pyp=fV0mom[1][1], pzp=fV0mom[1][2];
  Double_t px0=pxn+pxp, py0=pyn+pyp, pz0=pzn+pzp;
  Double_t p0=TMath::Sqrt(px0*px0 + py0*py0 + pz0*pz0);

  Double_t e0=TMath::Sqrt(1.11568*1.11568 + p0*p0);
  Double_t beta0=p0/e0;
  Double_t pln=(pxn*px0 + pyn*py0 + pzn*pz0)/p0;
  Double_t plp=(pxp*px0 + pyp*py0 + pzp*pz0)/p0;
  Double_t pt2=pxp*pxp + pyp*pyp + pzp*pzp - plp*plp;

  Double_t a=(plp-pln)/(plp+pln);
  a -= (pmass*pmass-nmass*nmass)/(1.11568*1.11568);
  a = 0.25*beta0*beta0*1.11568*1.11568*a*a + pt2;


  v0q=a - ps0*ps0;


  Double_t pxb=fBachMom[0], pyb=fBachMom[1], pzb=fBachMom[2]; 

  Double_t eb=TMath::Sqrt(bmass*bmass + pxb*pxb + pyb*pyb + pzb*pzb);
  Double_t pxl=px0+pxb, pyl=py0+pyb, pzl=pz0+pzb;
  Double_t pl=TMath::Sqrt(pxl*pxl + pyl*pyl + pzl*pzl);
  
  fEffMass=TMath::Sqrt((e0+eb)*(e0+eb) - pl*pl);

  Double_t beta=pl/(e0+eb);
  Double_t pl0=(px0*pxl + py0*pyl + pz0*pzl)/pl;
  Double_t plb=(pxb*pxl + pyb*pyl + pzb*pzl)/pl;
  pt2=p0*p0 - pl0*pl0;

  a=(pl0-plb)/(pl0+plb);
  a -= (1.11568*1.11568-bmass*bmass)/(mass*mass);
  a = 0.25*beta*beta*mass*mass*a*a + pt2;

  return (a - ps*ps);
}

void 
AliESDcascade::GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
  //--------------------------------------------------------------------
  // This function returns the cascade momentum (global)
  //--------------------------------------------------------------------
  px=fV0mom[0][0]+fV0mom[1][0]+fBachMom[0]; 
  py=fV0mom[0][1]+fV0mom[1][1]+fBachMom[1]; 
  pz=fV0mom[0][2]+fV0mom[1][2]+fBachMom[2]; 
}

void AliESDcascade::GetXYZ(Double_t &x, Double_t &y, Double_t &z) const {
  //--------------------------------------------------------------------
  // This function returns cascade position (global)
  //--------------------------------------------------------------------
  x=fPos[0]; 
  y=fPos[1]; 
  z=fPos[2]; 
}

Double_t AliESDcascade::GetD(Double_t x0, Double_t y0, Double_t z0) const {
  //--------------------------------------------------------------------
  // This function returns the cascade impact parameter
  //--------------------------------------------------------------------

  Double_t x=fPos[0],y=fPos[1],z=fPos[2];
  Double_t px=fV0mom[0][0]+fV0mom[1][0]+fBachMom[0];
  Double_t py=fV0mom[0][1]+fV0mom[1][1]+fBachMom[1];
  Double_t pz=fV0mom[0][2]+fV0mom[1][2]+fBachMom[2];

  Double_t dx=(y0-y)*pz - (z0-z)*py; 
  Double_t dy=(x0-x)*pz - (z0-z)*px;
  Double_t dz=(x0-x)*py - (y0-y)*px;
  Double_t d=TMath::Sqrt((dx*dx+dy*dy+dz*dz)/(px*px+py*py+pz*pz));

  return d;
}

