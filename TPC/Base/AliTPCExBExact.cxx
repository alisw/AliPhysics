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

////
// This implementation AliTPCExB is using an "exact" calculation of the ExB
// effect. That is, it uses the drift ODE to calculate the distortion
// without any further assumption.
// Due to the numerical integration of the ODE, there are some numerical
// uncertencies involed.
////

#include "TMath.h"
#include "TTreeStream.h"
#include "AliMagF.h"
#include "AliTPCExBExact.h"

ClassImp(AliTPCExBExact)

const Double_t AliTPCExBExact::fgkEM=1.602176487e-19/9.10938215e-31;
const Double_t AliTPCExBExact::fgkDriftField=-40.e3;

AliTPCExBExact::AliTPCExBExact()
  : fDriftVelocity(0),
    //fkMap(0),
    fkField(0),fkN(0),
    fkNX(0),fkNY(0),fkNZ(0),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.),
    fkNLook(0),fkLook(0) {
  //
  // purely for I/O
  //
}

AliTPCExBExact::AliTPCExBExact(const AliMagF *bField,
			       Double_t driftVelocity,
			       Int_t nx,Int_t ny,Int_t nz,Int_t n)
  : fDriftVelocity(driftVelocity),
    //fkMap(0),
    fkField(bField),fkN(n),
    fkNX(nx),fkNY(ny),fkNZ(nz),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.),
    fkNLook(0),fkLook(0) {
  //
  // The constructor. One has to supply a magnetic field and an (initial)
  // drift velocity. Since some kind of lookuptable is created the
  // number of its meshpoints can be supplied.
  // n sets the number of integration steps to be used when integrating
  // over the full drift length.
  //
  CreateLookupTable();
}

/*
AliTPCExBExact::AliTPCExBExact(const AliFieldMap *bFieldMap,
			       Double_t driftVelocity,Int_t n) 
  : fDriftVelocity(driftVelocity),
    fkMap(bFieldMap),fkField(0),fkN(n), 
    fkNX(0),fkNY(0),fkNZ(0),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.),
    fkNLook(0),fkLook(0) {
  //
  // The constructor. One has to supply a field map and an (initial)
  // drift velocity.
  // n sets the number of integration steps to be used when integrating
  // over the full drift length.
  //

  fkXMin=bFieldMap->Xmin()
    -TMath::Ceil( (bFieldMap->Xmin()+250.0)/bFieldMap->DelX())
    *bFieldMap->DelX();
  fkXMax=bFieldMap->Xmax()
    -TMath::Floor((bFieldMap->Xmax()-250.0)/bFieldMap->DelX())
    *bFieldMap->DelX();
  fkYMin=bFieldMap->Ymin()
    -TMath::Ceil( (bFieldMap->Ymin()+250.0)/bFieldMap->DelY())
    *bFieldMap->DelY();
  fkYMax=bFieldMap->Ymax()
    -TMath::Floor((bFieldMap->Ymax()-250.0)/bFieldMap->DelY())
    *bFieldMap->DelY();
  fkZMax=bFieldMap->Zmax()
    -TMath::Floor((bFieldMap->Zmax()-250.0)/bFieldMap->DelZ())
    *bFieldMap->DelZ();
  fkZMax=TMath::Max(0.,fkZMax); // I really hope that this is unnecessary!

  fkNX=static_cast<Int_t>((fkXMax-fkXMin)/bFieldMap->DelX()+1.1);
  fkNY=static_cast<Int_t>((fkYMax-fkYMin)/bFieldMap->DelY()+1.1);
  fkNZ=static_cast<Int_t>((fkZMax-fkZMin)/bFieldMap->DelZ()+1.1);

  CreateLookupTable();
}
*/

AliTPCExBExact::~AliTPCExBExact() {
  //
  // destruct the poor object.
  //
  delete[] fkLook;
}

void AliTPCExBExact::Correct(const Double_t *position, Double_t *corrected) {
  //
  // correct for the distortion
  //
  Double_t x=(position[0]-fkXMin)/(fkXMax-fkXMin)*(fkNX-1);
  Int_t xi1=static_cast<Int_t>(x);
  xi1=TMath::Max(TMath::Min(xi1,fkNX-2),0);
  Int_t xi2=xi1+1;
  Double_t dx=(x-xi1);
  Double_t dx1=(xi2-x);
  
  Double_t y=(position[1]-fkYMin)/(fkYMax-fkYMin)*(fkNY-1);
  Int_t yi1=static_cast<Int_t>(y);
  yi1=TMath::Max(TMath::Min(yi1,fkNY-2),0);
  Int_t yi2=yi1+1;
  Double_t dy=(y-yi1);
  Double_t dy1=(yi2-y);
  
  Double_t z=position[2]/fkZMax*(fkNZ-1);
  Int_t side;
  if (z>0) {
    side=1;
  }
  else {
    z=-z;
    side=0;
  }
  Int_t zi1=static_cast<Int_t>(z);
  zi1=TMath::Max(TMath::Min(zi1,fkNZ-2),0);
  Int_t zi2=zi1+1;
  Double_t dz=(z-zi1);
  Double_t dz1=(zi2-z);
  
  for (int i=0;i<3;++i)
    corrected[i]
      =fkLook[(((xi1*fkNY+yi1)*fkNZ+zi1)*2+side)*3+i]*dx1*dy1*dz1
      +fkLook[(((xi1*fkNY+yi1)*fkNZ+zi2)*2+side)*3+i]*dx1*dy1*dz
      +fkLook[(((xi1*fkNY+yi2)*fkNZ+zi1)*2+side)*3+i]*dx1*dy *dz1
      +fkLook[(((xi1*fkNY+yi2)*fkNZ+zi2)*2+side)*3+i]*dx1*dy *dz
      +fkLook[(((xi2*fkNY+yi2)*fkNZ+zi1)*2+side)*3+i]*dx *dy *dz1
      +fkLook[(((xi2*fkNY+yi2)*fkNZ+zi2)*2+side)*3+i]*dx *dy *dz
      +fkLook[(((xi2*fkNY+yi1)*fkNZ+zi1)*2+side)*3+i]*dx *dy1*dz1
      +fkLook[(((xi2*fkNY+yi1)*fkNZ+zi2)*2+side)*3+i]*dx *dy1*dz ;
  //    corrected[2]=position[2];
}

/*
void AliTPCExBExact::TestThisBeautifulObject(const AliFieldMap *bFieldMap,
					     const char* fileName) {
  //
  // Have a look at the common part "TestThisBeautifulObjectGeneric".
  //
  fkMap=bFieldMap;
  fkField=0;
  TestThisBeautifulObjectGeneric(fileName);
}
*/

void AliTPCExBExact::TestThisBeautifulObject(const AliMagF *bField,
					     const char* fileName) {
  //
  // Have a look at the common part "TestThisBeautifulObjectGeneric".
  //
  fkField=bField;
  //fkMap=0;
  TestThisBeautifulObjectGeneric(fileName);
}

void AliTPCExBExact::TestThisBeautifulObjectGeneric(const char* fileName) {
  //
  // Well, as the name sais... it tests the object.
  //
  TTreeSRedirector ts(fileName);
  Double_t x[3];
  for (x[0]=-250.;x[0]<=250.;x[0]+=10.)
    for (x[1]=-250.;x[1]<=250.;x[1]+=10.)
      for (x[2]=-250.;x[2]<=250.;x[2]+=10.) {
	Double_t d[3];
	Double_t dnl[3];
	Correct(x,d);
	CalculateDistortion(x,dnl);
	Double_t r=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
	Double_t rd=TMath::Sqrt(d[0]*d[0]+d[1]*d[1]);
	Double_t dr=r-rd;
	Double_t phi=TMath::ATan2(x[0],x[1]);
	Double_t phid=TMath::ATan2(d[0],d[1]);
	Double_t dphi=phi-phid;
	if (dphi<0.) dphi+=TMath::TwoPi();
	if (dphi>TMath::Pi()) dphi=TMath::TwoPi()-dphi;
	Double_t drphi=r*dphi;
	Double_t dx=x[0]-d[0];
	Double_t dy=x[1]-d[1];
	Double_t dz=x[2]-d[2];
	Double_t dnlx=x[0]-dnl[0];
	Double_t dnly=x[1]-dnl[1];
	Double_t dnlz=x[2]-dnl[2];
	Double_t b[3];
	GetB(b,x);
	ts<<"positions"
	  <<"bx="<<b[0]
	  <<"by="<<b[1]
	  <<"bz="<<b[2]
	  <<"x0="<<x[0]
	  <<"x1="<<x[1]
	  <<"x2="<<x[2]
	  <<"dx="<<dx
	  <<"dy="<<dy
	  <<"dz="<<dz
	  <<"dnlx="<<dnlx
	  <<"dnly="<<dnly
	  <<"dnlz="<<dnlz
	  <<"r="<<r
	  <<"phi="<<phi
	  <<"dr="<<dr
	  <<"drphi="<<drphi
	  <<"\n";
      }
}

void AliTPCExBExact::CreateLookupTable() {
  //
  // Helper function to fill the lookup table.
  //
  fkNLook=fkNX*fkNY*fkNZ*2*3;
  fkLook=new Double_t[fkNLook];
  Double_t x[3];
  for (int i=0;i<fkNX;++i) {
    x[0]=fkXMin+(fkXMax-fkXMin)/(fkNX-1)*i;
    for (int j=0;j<fkNY;++j) {
      x[1]=fkYMin+(fkYMax-fkYMin)/(fkNY-1)*j;
      for (int k=0;k<fkNZ;++k) {
	x[2]=1.*fkZMax/(fkNZ-1)*k;
	x[2]=TMath::Max((Double_t)0.0001,x[2]); //ugly
	CalculateDistortion(x,&fkLook[(((i*fkNY+j)*fkNZ+k)*2+1)*3]);
	x[2]=-x[2];
	CalculateDistortion(x,&fkLook[(((i*fkNY+j)*fkNZ+k)*2+0)*3]);
      }
    }
  }
}

void AliTPCExBExact::GetE(Double_t *e,const Double_t *x) const {
  //
  // Helper function returning the E field in SI units (V/m).
  //
  e[0]=0.;
  e[1]=0.;
  e[2]=(x[2]<0.?-1.:1.)*fgkDriftField; // in V/m
}

void AliTPCExBExact::GetB(Double_t *b,const Double_t *x) const {
  //
  // Helper function returning the B field in SI units (T).
  //
  Double_t xm[3];
  // the beautiful m to cm (and the ugly "const_cast") and Double_t 
  // to Float_t read the NRs introduction!:
  for (int i=0;i<3;++i) xm[i]=x[i]*100.;
  Double_t bf[3];
  //if (fkMap!=0)
  //  fkMap->Field(xm,bf);
  //else
  ((AliMagF*)fkField)->Field(xm,bf);
  for (int i=0;i<3;++i) b[i]=bf[i]/10.;
}

void AliTPCExBExact::Motion(const Double_t *x,Double_t,
			    Double_t *dxdt) const {
  //
  // The differential equation of motion of the electrons.
  //
  Double_t tau=fDriftVelocity/fgkDriftField/fgkEM;
  Double_t tau2=tau*tau;
  Double_t e[3];
  Double_t b[3];
  GetE(e,x);
  GetB(b,x);
  Double_t wx=fgkEM*b[0];
  Double_t wy=fgkEM*b[1];
  Double_t wz=fgkEM*b[2];
  Double_t ex=fgkEM*e[0];
  Double_t ey=fgkEM*e[1];
  Double_t ez=fgkEM*e[2];
  Double_t w2=(wx*wx+wy*wy+wz*wz);
  dxdt[0]=(1.+wx*wx*tau2)*ex+(wz*tau+wx*wy*tau2)*ey+(-wy*tau+wx*wz*tau2)*ez;
  dxdt[1]=(-wz*tau+wx*wy*tau2)*ex+(1.+wy*wy*tau2)*ey+(wx*tau+wy*wz*tau2)*ez;
  dxdt[2]=(wy*tau+wx*wz*tau2)*ex+(-wx*tau+wy*wz*tau2)*ey+(1.+wz*wz*tau2)*ez;
  Double_t fac=tau/(1.+w2*tau2);
  dxdt[0]*=fac;
  dxdt[1]*=fac;
  dxdt[2]*=fac;
}

void AliTPCExBExact::CalculateDistortion(const Double_t *x0,
					 Double_t *dist) const {
  //
  // Helper function that calculates one distortion by integration
  // (only used to fill the lookup table).
  //
  Double_t h=0.01*250./fDriftVelocity/fkN;
  Double_t t=0.;
  Double_t xt[3];
  Double_t xo[3];
  for (int i=0;i<3;++i)
    xo[i]=xt[i]=x0[i]*0.01;
  while (TMath::Abs(xt[2])<250.*0.01) {
    for (int i=0;i<3;++i)
      xo[i]=xt[i];
    DGLStep(xt,t,h);
    t+=h;
  }
  if (t!=0.) {
    Double_t p=((xt[2]<0.?-1.:1.)*250.*0.01-xo[2])/(xt[2]-xo[2]);
    dist[0]=(xo[0]+p*(xt[0]-xo[0]))*100.;
    dist[1]=(xo[1]+p*(xt[1]-xo[1]))*100.;
    //    dist[2]=(xo[2]+p*(xt[2]-xo[2]))*100.;
    dist[2]=(x0[2]>0.?-1:1.)*(t-h+p*h)*fDriftVelocity*100.;
    dist[2]+=(x0[2]<0.?-1:1.)*250.;
  }
  else {
    dist[0]=x0[0];
    dist[1]=x0[1];
    dist[2]=x0[2];
  }
  // reverse the distortion, i.e. get the correction
  dist[0]=x0[0]-(dist[0]-x0[0]);
  dist[1]=x0[1]-(dist[1]-x0[1]);
}

void AliTPCExBExact::DGLStep(Double_t *x,Double_t t,Double_t h) const {
  //
  // An elementary integration step.
  // (simple Euler Method)
  //
  Double_t dxdt[3];
  Motion(x,t,dxdt);
  for (int i=0;i<3;++i)
    x[i]+=h*dxdt[i];

  /* suggestions about how to write it this way are welcome!
     void DGLStep(void (*f)(const Double_t *x,Double_t t,Double_t *dxdt),
                   Double_t *x,Double_t t,Double_t h,Int_t n) const;
  */

}
