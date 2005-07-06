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

////////////////////////////////////////////////////////////////////////////
//               Implementation of the cascade vertex class               //
//                                                                        //
//  Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr//
////////////////////////////////////////////////////////////////////////////
#include <TMath.h>

#include "AliCascadeVertex.h"
#include "AliITStrackV2.h"
#include "AliV0vertex.h"

ClassImp(AliCascadeVertex)


AliCascadeVertex::AliCascadeVertex(const AliV0vertex &v,const AliITStrackV2 &t) {
  //--------------------------------------------------------------------
  // Main constructor
  //--------------------------------------------------------------------
  fPdgCode=kXiMinus;

  fV0idx[0]=v.GetNindex(); fV0idx[1]=v.GetPindex();
  fBachIdx=t.GetLabel(); 

  //Trivial estimation of the vertex parameters
  Double_t pt, phi, x, par[5];
  Double_t alpha, cs, sn;

  t.GetExternalParameters(x,par); alpha=t.GetAlpha();
  pt=1./TMath::Abs(par[4]);
  phi=TMath::ASin(par[2]) + alpha;  

  // momentum of the bachelor track

  Double_t px1=pt*TMath::Cos(phi), py1=pt*TMath::Sin(phi), pz1=pt*par[3];

  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);

  Double_t x1=x*cs - par[0]*sn; // position of the bachelor at dca (bachelor,V0)
  Double_t y1=x*sn + par[0]*cs;
  Double_t z1=par[1];

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


  fChi2=7.;   

}


