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

//-------------------------------------------------------------------------
//               Implementation of the V0 vertex class
//
//     Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------
#include <Riostream.h>
#include <TMath.h>
#include <TPDGCode.h>

#include "AliV0vertex.h"
#include "AliITStrackV2.h"

ClassImp(AliV0vertex)

AliV0vertex::AliV0vertex() : TObject() {
  //--------------------------------------------------------------------
  // Default constructor  (K0s)
  //--------------------------------------------------------------------
  fPdgCode=kK0Short;
  fEffMass=0.497672;
  fChi2=1.e+33;
  fPos[0]=fPos[1]=fPos[2]=0.;
  fPosCov[0]=fPosCov[1]=fPosCov[2]=fPosCov[3]=fPosCov[4]=fPosCov[5]=0.;
}

AliV0vertex::AliV0vertex(const AliITStrackV2 &n, const AliITStrackV2 &p) {
  //--------------------------------------------------------------------
  // Main constructor
  //--------------------------------------------------------------------
  fPdgCode=kK0Short;
  fNlab=n.GetLabel(); fPlab=p.GetLabel();

  //Trivial estimation of the vertex parameters
  Double_t pt, phi, x, par[5];
  Double_t alpha, cs, sn;

  n.GetExternalParameters(x,par); alpha=n.GetAlpha();
  pt=1./TMath::Abs(par[4]);
  phi=TMath::ASin(par[2]) + alpha;  
  Double_t px1=pt*TMath::Cos(phi), py1=pt*TMath::Sin(phi), pz1=pt*par[3];
  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  Double_t x1=x*cs - par[0]*sn;
  Double_t y1=x*sn + par[0]*cs;
  Double_t z1=par[1];
  Double_t sx1=sn*sn*n.GetSigmaY2(), sy1=cs*cs*n.GetSigmaY2(); 

  p.GetExternalParameters(x,par); alpha=p.GetAlpha();
  pt=1./TMath::Abs(par[4]);
  phi=TMath::ASin(par[2]) + alpha;  
  Double_t px2=pt*TMath::Cos(phi), py2=pt*TMath::Sin(phi), pz2=pt*par[3];
  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  Double_t x2=x*cs - par[0]*sn;
  Double_t y2=x*sn + par[0]*cs;
  Double_t z2=par[1];
  Double_t sx2=sn*sn*p.GetSigmaY2(), sy2=cs*cs*p.GetSigmaY2(); 
    
  Double_t sz1=n.GetSigmaZ2(), sz2=p.GetSigmaZ2();
  Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
  Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
  Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
  fPos[0]=wx1*x1 + wx2*x2; fPos[1]=wy1*y1 + wy2*y2; fPos[2]=wz1*z1 + wz2*z2;

  //fPos[0]=0.5*(x1+x2); fPos[1]=0.5*(y1+y2); fPos[2]=0.5*(z1+z2);
  fNmom[0]=px1; fNmom[1]=py1; fNmom[2]=pz1; 
  fPmom[0]=px2; fPmom[1]=py2; fPmom[2]=pz2;

  Double_t e1=TMath::Sqrt(0.13957*0.13957 + px1*px1 + py1*py1 + pz1*pz1);
  Double_t e2=TMath::Sqrt(0.13957*0.13957 + px2*px2 + py2*py2 + pz2*pz2);
  fEffMass=TMath::Sqrt((e1+e2)*(e1+e2)-
    (px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));

  fChi2=7.;   
}
/*
Double_t AliV0vertex::ChangeMassHypothesis(Int_t code) {
  //--------------------------------------------------------------------
  // This function changes the mass hypothesis for this V0
  // and returns the "kinematical quality" of this hypothesis 
  //--------------------------------------------------------------------
  Double_t nmass=0.13957, pmass=0.13957, mass=0.49767, des=0;

  fPdgCode=code;

  switch (code) {
  case kLambda0:
    nmass=0.13957; pmass=0.93827; mass=1.1157; des=0.9437-0.1723; break;
  case kLambda0Bar:
    pmass=0.13957; nmass=0.93827; mass=1.1157; des=0.1723-0.9437; break;
  case kK0Short: 
    break;
  default:
    cerr<<"AliV0vertex::ChangeMassHypothesis: ";
    cerr<<"invalide PDG code ! Assuming K0s...\n";
    fPdgCode=kK0Short;
    break;
  }

  Double_t pxn=fNmom[0], pyn=fNmom[1], pzn=fNmom[2]; 
  Double_t pxp=fPmom[0], pyp=fPmom[1], pzp=fPmom[2];

  Double_t en=TMath::Sqrt(nmass*nmass + pxn*pxn + pyn*pyn + pzn*pzn);
  Double_t ep=TMath::Sqrt(pmass*pmass + pxp*pxp + pyp*pyp + pzp*pzp);
  Double_t pxl=pxn+pxp, pyl=pyn+pyp, pzl=pzn+pzp;
  Double_t pl=TMath::Sqrt(pxl*pxl + pyl*pyl + pzl*pzl);

  fEffMass=TMath::Sqrt((en+ep)*(en+ep)-pl*pl);

  Double_t gamma=(en+ep)/mass, betagamma=pl/mass;
  Double_t pln=(pxn*pxl + pyn*pyl + pzn*pzl)/pl;
  Double_t plp=(pxp*pxl + pyp*pyl + pzp*pzl)/pl;
  Double_t plps=gamma*plp - betagamma*ep;

  Double_t diff=2*gamma*plps + betagamma*des;

  return (plp-pln-diff);
  
}
*/
Double_t AliV0vertex::ChangeMassHypothesis(Int_t code) {
  //--------------------------------------------------------------------
  // This function changes the mass hypothesis for this V0
  // and returns the "kinematical quality" of this hypothesis 
  //--------------------------------------------------------------------
  Double_t nmass=0.13957, pmass=0.13957, mass=0.49767, ps=0.206;

  fPdgCode=code;

  switch (code) {
  case kLambda0:
    nmass=0.13957; pmass=0.93827; mass=1.1157; ps=0.101; break;
  case kLambda0Bar:
    pmass=0.13957; nmass=0.93827; mass=1.1157; ps=0.101; break;
  case kK0Short: 
    break;
  default:
    cerr<<"AliV0vertex::ChangeMassHypothesis: ";
    cerr<<"invalide PDG code ! Assuming K0s...\n";
    fPdgCode=kK0Short;
    break;
  }

  Double_t pxn=fNmom[0], pyn=fNmom[1], pzn=fNmom[2]; 
  Double_t pxp=fPmom[0], pyp=fPmom[1], pzp=fPmom[2];

  Double_t en=TMath::Sqrt(nmass*nmass + pxn*pxn + pyn*pyn + pzn*pzn);
  Double_t ep=TMath::Sqrt(pmass*pmass + pxp*pxp + pyp*pyp + pzp*pzp);
  Double_t pxl=pxn+pxp, pyl=pyn+pyp, pzl=pzn+pzp;
  Double_t pl=TMath::Sqrt(pxl*pxl + pyl*pyl + pzl*pzl);

  fEffMass=TMath::Sqrt((en+ep)*(en+ep)-pl*pl);

  Double_t beta=pl/(en+ep);
  Double_t pln=(pxn*pxl + pyn*pyl + pzn*pzl)/pl;
  Double_t plp=(pxp*pxl + pyp*pyl + pzp*pzl)/pl;

  Double_t pt2=pxp*pxp + pyp*pyp + pzp*pzp - plp*plp;

  Double_t a=(plp-pln)/(plp+pln);
  a -= (pmass*pmass-nmass*nmass)/(mass*mass);
  a = 0.25*beta*beta*mass*mass*a*a + pt2;

  return (a - ps*ps);
  
}

void AliV0vertex::GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
  //--------------------------------------------------------------------
  // This function returns V0's momentum (global)
  //--------------------------------------------------------------------
  px=fNmom[0]+fPmom[0]; 
  py=fNmom[1]+fPmom[1]; 
  pz=fNmom[2]+fPmom[2]; 
}

void AliV0vertex::GetXYZ(Double_t &x, Double_t &y, Double_t &z) const {
  //--------------------------------------------------------------------
  // This function returns V0's position (global)
  //--------------------------------------------------------------------
  x=fPos[0]; 
  y=fPos[1]; 
  z=fPos[2]; 
}

Double_t AliV0vertex::GetD(Double_t x0, Double_t y0, Double_t z0) const {
  //--------------------------------------------------------------------
  // This function returns V0's impact parameter
  //--------------------------------------------------------------------
  Double_t x=fPos[0],y=fPos[1],z=fPos[2];
  Double_t px=fNmom[0]+fPmom[0];
  Double_t py=fNmom[1]+fPmom[1];
  Double_t pz=fNmom[2]+fPmom[2];

  Double_t dx=(y0-y)*pz - (z0-z)*py; 
  Double_t dy=(x0-x)*pz - (z0-z)*px;
  Double_t dz=(x0-x)*py - (y0-y)*px;
  Double_t d=TMath::Sqrt((dx*dx+dy*dy+dz*dz)/(px*px+py*py+pz*pz));
  return d;
}
