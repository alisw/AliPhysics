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

#include <TDatabasePDG.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"

ClassImp(AliESDcascade)

AliESDcascade::AliESDcascade() : 
  TObject(),
  fPdgCode(kXiMinus),
  fEffMass(TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass()),
  fChi2(1.e+33),
  fBachIdx(-1)
{
  //--------------------------------------------------------------------
  // Default constructor  (Xi-)
  //--------------------------------------------------------------------
  for (Int_t j=0; j<3; j++) {
    fPos[j]=0.;
    fBachMom[j]=0.;
  }

  for (Int_t i=0; i<2; i++)
    for (Int_t j=0; j<3; j++)
      fV0mom[i][j]=0.;

  fV0idx[0]=fV0idx[1]=-1;

  fPosCov[0]=1e10;
  fPosCov[1]=fPosCov[2]=0.;
  fPosCov[3]=1e10;
  fPosCov[4]=0.;
  fPosCov[5]=1e10;

  fV0momCov[0]=1e10;
  fV0momCov[1]=fV0momCov[2]=0.;
  fV0momCov[3]=1e10;
  fV0momCov[4]=0.;
  fV0momCov[5]=1e10;

  fBachMomCov[0]=1e10;
  fBachMomCov[1]=fBachMomCov[2]=0.;
  fBachMomCov[3]=1e10;
  fBachMomCov[4]=0.;
  fBachMomCov[5]=1e10;
}

AliESDcascade::AliESDcascade(const AliESDv0 &v,
			     const AliExternalTrackParam &t, Int_t i) : 
  TObject(),
  fPdgCode(kXiMinus),
  fEffMass(TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass()),
  fChi2(1.e+33),
  fBachIdx(i)
{
  //--------------------------------------------------------------------
  // Main constructor  (Xi-)
  //--------------------------------------------------------------------

  fV0idx[0]=v.GetNindex(); fV0idx[1]=v.GetPindex();

  Double_t r[3]; t.GetXYZ(r);
  Double_t x1=r[0], y1=r[1], z1=r[2]; // position of the bachelor
  Double_t p[3]; t.GetPxPyPz(p);
  Double_t px1=p[0], py1=p[1], pz1=p[2];// momentum of the bachelor track

  Double_t x2,y2,z2;          // position of the V0 
  v.GetXYZ(x2,y2,z2);    
  Double_t px2,py2,pz2;       // momentum of V0
  v.GetPxPyPz(px2,py2,pz2);

  Double_t a2=((x1-x2)*px2+(y1-y2)*py2+(z1-z2)*pz2)/(px2*px2+py2*py2+pz2*pz2);

  Double_t xm=x2+a2*px2;
  Double_t ym=y2+a2*py2;
  Double_t zm=z2+a2*pz2;

  // position of the cascade decay
  
  fPos[0]=0.5*(x1+xm); fPos[1]=0.5*(y1+ym); fPos[2]=0.5*(z1+zm);
    

  // invariant mass of the cascade (default is Ximinus)
  
  Double_t e1=TMath::Sqrt(0.13957*0.13957 + px1*px1 + py1*py1 + pz1*pz1);
  Double_t e2=TMath::Sqrt(1.11568*1.11568 + px2*px2 + py2*py2 + pz2*pz2);
  
  fEffMass=TMath::Sqrt((e1+e2)*(e1+e2)-
    (px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));


  // momenta of the bachelor and the V0
  
  fBachMom[0]=px1; fBachMom[1]=py1; fBachMom[2]=pz1; 
  v.GetNPxPyPz(px2,py2,pz2);
  fV0mom[0][0]=px2; fV0mom[0][1]=py2; fV0mom[0][2]=pz2;
  v.GetPPxPyPz(px2,py2,pz2);
  fV0mom[1][0]=px2; fV0mom[1][1]=py2; fV0mom[1][2]=pz2;

  //PH Covariance matrices: to be calculated correctly in the future
  fPosCov[0]=1e10;
  fPosCov[1]=fPosCov[2]=0.;
  fPosCov[3]=1e10;
  fPosCov[4]=0.;
  fPosCov[5]=1e10;

  fV0momCov[0]=1e10;
  fV0momCov[1]=fV0momCov[2]=0.;
  fV0momCov[3]=1e10;
  fV0momCov[4]=0.;
  fV0momCov[5]=1e10;

  fBachMomCov[0]=1e10;
  fBachMomCov[1]=fBachMomCov[2]=0.;
  fBachMomCov[3]=1e10;
  fBachMomCov[4]=0.;
  fBachMomCov[5]=1e10;

  fChi2=7.;   

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
       AliError("Invalide PDG code !  Assuming XiMinus's...");
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

