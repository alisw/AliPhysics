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

//
// This implementation AliTPCExB is using an aproximate calculation of the ExB
// effect. Therefore the drift ODE is Taylor expanded and only the first
// order part is taken.
//
// 
// The ExB correction map is stored in the calib DB
// To test current version:
/*
  char *storage = "local://OCDBres"
  Int_t RunNumber=0;
  AliCDBManager::Instance()->SetDefaultStorage(storage);
  AliCDBManager::Instance()->SetRun(RunNumber) 
  AliTPCExBFirst * exb = AliTPCcalibDB::Instance()->GetExB();

  //  exb->TestExB("exb.root");
  // TFile f("exb.root");
  //positions->Draw("drphi");
*/



#include "TMath.h"
//#include "AliFieldMap.h"
#include "AliMagF.h"
#include "TTreeStream.h"
#include "AliTPCExBFirst.h"

ClassImp(AliTPCExBFirst)

const Double_t AliTPCExBFirst::fgkEM=1.602176487e-19/9.10938215e-31;
const Double_t AliTPCExBFirst::fgkDriftField=-40.e3;

AliTPCExBFirst::AliTPCExBFirst()
  : fDriftVelocity(0),
    fkNX(0),fkNY(0),fkNZ(0),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.),
    fkNMean(0),
    fkMeanBx(0),fkMeanBy(0),fkMeanBz(0.) {
  //
  // purely for I/O
  //
  SetInstance(this);
}

AliTPCExBFirst::AliTPCExBFirst(const AliMagF *bField,
			       Double_t driftVelocity,
			       Int_t nx,Int_t ny,Int_t nz)
  : fDriftVelocity(driftVelocity),
    fkNX(nx),fkNY(ny),fkNZ(nz),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.),
    fkNMean(0),
    fkMeanBx(0),fkMeanBy(0),fkMeanBz(0.) {
  //
  // The constructor. One has to supply a magnetic field and an (initial)
  // drift velocity. Since some kind of lookuptable is created the
  // number of its meshpoints can be supplied.
  //
  //  ConstructCommon(0,bField);
  ConstructCommon(bField);
  SetInstance(this);
}

/*
AliTPCExBFirst::AliTPCExBFirst(const AliFieldMap *bFieldMap,
			       Double_t driftVelocity) 
  : fDriftVelocity(driftVelocity),
    fkNX(0),fkNY(0),fkNZ(0),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.),
    fkNMean(0),
    fkMeanBx(0),fkMeanBy(0),fkMeanBz(0.) {
  //
  // The constructor. One has to supply a field map and an (initial)
  // drift velocity.
  //
  SetInstance(this);
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
  fkZMin=bFieldMap->Zmin()
    -TMath::Ceil( (bFieldMap->Zmin()+250.0)/bFieldMap->DelZ())
    *bFieldMap->DelZ();
  fkZMax=bFieldMap->Zmax()
    -TMath::Floor((bFieldMap->Zmax()-250.0)/bFieldMap->DelZ())
    *bFieldMap->DelZ();

  fkNX=static_cast<Int_t>((fkXMax-fkXMin)/bFieldMap->DelX()+1.1);
  fkNY=static_cast<Int_t>((fkYMax-fkYMin)/bFieldMap->DelY()+1.1);
  fkNZ=static_cast<Int_t>((fkZMax-fkZMin)/bFieldMap->DelZ()+1.1);

  ConstructCommon(bFieldMap,0);
}
*/

AliTPCExBFirst::~AliTPCExBFirst() { 
  //
  // destruct the poor object.
  //
  delete[] fkMeanBx;
  delete[] fkMeanBy;
}

void AliTPCExBFirst::Correct(const Double_t *position,Double_t *corrected) {
  //
  // correct for the distortion
  //
  Double_t bx,by;
  GetMeanFields(position[0],position[1],position[2],&bx,&by);
  if (position[2]>0.) {
    Double_t bxe,bye;
    GetMeanFields(position[0],position[1],250.,&bxe,&bye);
    if (position[2]!=250.) {
      bx=(500.*bxe-(position[2]+250.)*bx)/(250.-position[2]);
      by=(500.*bye-(position[2]+250.)*by)/(250.-position[2]);
    }
    else {
      bx=bxe;
      by=bye;
    }
  }
  
  Double_t mu=fDriftVelocity/fgkDriftField;
  Double_t wt=mu*fkMeanBz;
  
  corrected[0]=mu*(wt*bx-by)/(1.+wt*wt);
  corrected[1]=mu*(wt*by+bx)/(1.+wt*wt);
  
  if (position[2]>0.) {
    corrected[0]*=(250.-position[2]);
    corrected[1]*=(250.-position[2]);
  }
  else {
    corrected[0]*=(-250.-position[2]);
    corrected[1]*=(-250.-position[2]);
  }

  corrected[0]=position[0]-corrected[0];
  corrected[1]=position[1]-corrected[1];
  corrected[2]=position[2];
}

void AliTPCExBFirst::TestThisBeautifulObject(const char* fileName) {
  //
  // well, as the name sais...
  //
  TTreeSRedirector ts(fileName);
  Double_t x[3];
  for (x[0]=-250.;x[0]<=250.;x[0]+=10.)
    for (x[1]=-250.;x[1]<=250.;x[1]+=10.)
      for (x[2]=-250.;x[2]<=250.;x[2]+=10.) {
	Double_t d[3];
	Correct(x,d);
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
	ts<<"positions"
	  <<"x0="<<x[0]
	  <<"x1="<<x[1]
	  <<"x2="<<x[2]
	  <<"dx="<<dx
	  <<"dy="<<dy
	  <<"dz="<<dz
	  <<"r="<<r
	  <<"phi="<<phi
	  <<"dr="<<dr
	  <<"drphi="<<drphi
	  <<"\n";
      }
}


void AliTPCExBFirst::ConstructCommon(//const AliFieldMap *bFieldMap,
				     const AliMagF *bField) {
  //
  // THIS IS PRIVATE! (a helper for the constructor)
  //
  fkNMean=fkNX*fkNY*fkNZ;
  fkMeanBx=new Double_t[fkNMean];
  fkMeanBy=new Double_t[fkNMean];

  Double_t x[3];
  Double_t nBz=0;
  fkMeanBz=0.;
  for (int i=0;i<fkNX;++i) {
    x[0]=fkXMin+i*(fkXMax-fkXMin)/(fkNX-1);
    for (int j=0;j<fkNY;++j) {
      x[1]=fkYMin+j*(fkYMax-fkYMin)/(fkNY-1);
      Double_t r=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
      Double_t bx=0.,by=0.;
      for (int k=0;k<fkNZ;++k) {
	x[2]=fkZMin+k*(fkZMax-fkZMin)/(fkNZ-1);
	Double_t b[3];
	// the x is not const in the Field function...
	Double_t xt[3];
	for (int l=0;l<3;++l) xt[l]=x[l];
	// that happens due to the lack of a sophisticated class design:
	//	if (bFieldMap!=0)
	//	  bFieldMap->Field(xt,b);
	//	else 
	((AliMagF*)bField)->Field(xt,b);
	bx+=b[0]/10.;
	by+=b[1]/10.;
	fkMeanBx[(k*fkNY+j)*fkNX+i]=bx;
	fkMeanBy[(k*fkNY+j)*fkNX+i]=by;
	if (90.<=r&&r<=250.) {
	  fkMeanBz+=b[2]/10.;
	  ++nBz;
	}
      }
    }
  }
  fkMeanBz/=nBz;
}


void AliTPCExBFirst::GetMeanFields(Double_t rx,Double_t ry,Double_t rz,
				   Double_t *Bx,Double_t *By) const {
  //
  // THIS IS PRIVATE! (calculates the mean field utilising a lookup table)
  //
  Double_t x=(fkNX-1)*(rx-fkXMin)/(fkXMax-fkXMin);
  Int_t xi1=static_cast<Int_t>(x);
  xi1=TMath::Max(TMath::Min(xi1,fkNX-2),0);
  Int_t xi2=xi1+1;
  Double_t dx=(x-xi1);
  Double_t dx1=(xi2-x);

  Double_t y=(fkNY-1)*(ry-fkYMin)/(fkYMax-fkYMin);
  Int_t yi1=static_cast<Int_t>(y);
  yi1=TMath::Max(TMath::Min(yi1,fkNY-2),0);
  Int_t yi2=yi1+1;
  Double_t dy=(y-yi1);
  Double_t dy1=(yi2-y);
  
  Double_t z=(fkNZ-1)*(rz-fkZMin)/(fkZMax-fkZMin);
  Int_t zi1=static_cast<Int_t>(z);
  zi1=TMath::Max(TMath::Min(zi1,fkNZ-2),0);
  Int_t zi2=zi1+1;
  Double_t dz=(z-zi1);

  double s0x=fkMeanBx[yi1*fkNX+xi1]*dx1*dy1
	    +fkMeanBx[yi2*fkNX+xi1]*dx1*dy
            +fkMeanBx[yi1*fkNX+xi2]*dx *dy1
            +fkMeanBx[yi2*fkNX+xi2]*dx *dy;
  double s0y=fkMeanBy[yi1*fkNX+xi1]*dx1*dy1
	    +fkMeanBy[yi2*fkNX+xi1]*dx1*dy
            +fkMeanBy[yi1*fkNX+xi2]*dx *dy1
            +fkMeanBy[yi2*fkNX+xi2]*dx *dy;
  Int_t zi0=zi1-1;
  double snmx,snmy;
  if (zi0>=0) {
    snmx=fkMeanBx[(zi0*fkNY+yi1)*fkNX+xi1]*dx1*dy1
        +fkMeanBx[(zi0*fkNY+yi2)*fkNX+xi1]*dx1*dy
        +fkMeanBx[(zi0*fkNY+yi1)*fkNX+xi2]*dx *dy1
        +fkMeanBx[(zi0*fkNY+yi2)*fkNX+xi2]*dx *dy;
    snmy=fkMeanBy[(zi0*fkNY+yi1)*fkNX+xi1]*dx1*dy1
        +fkMeanBy[(zi0*fkNY+yi2)*fkNX+xi1]*dx1*dy
        +fkMeanBy[(zi0*fkNY+yi1)*fkNX+xi2]*dx *dy1
        +fkMeanBy[(zi0*fkNY+yi2)*fkNX+xi2]*dx *dy;
  }
  else
    snmx=snmy=0.;
  double snx=fkMeanBx[(zi1*fkNY+yi1)*fkNX+xi1]*dx1*dy1
 	    +fkMeanBx[(zi1*fkNY+yi2)*fkNX+xi1]*dx1*dy
            +fkMeanBx[(zi1*fkNY+yi1)*fkNX+xi2]*dx *dy1
            +fkMeanBx[(zi1*fkNY+yi2)*fkNX+xi2]*dx *dy;
  double sny=fkMeanBy[(zi1*fkNY+yi1)*fkNX+xi1]*dx1*dy1
 	    +fkMeanBy[(zi1*fkNY+yi2)*fkNX+xi1]*dx1*dy
            +fkMeanBy[(zi1*fkNY+yi1)*fkNX+xi2]*dx *dy1
            +fkMeanBy[(zi1*fkNY+yi2)*fkNX+xi2]*dx *dy;
  double snpx=fkMeanBx[(zi2*fkNY+yi1)*fkNX+xi1]*dx1*dy1
	     +fkMeanBx[(zi2*fkNY+yi2)*fkNX+xi1]*dx1*dy
             +fkMeanBx[(zi2*fkNY+yi1)*fkNX+xi2]*dx *dy1
             +fkMeanBx[(zi2*fkNY+yi2)*fkNX+xi2]*dx *dy;
  double snpy=fkMeanBy[(zi2*fkNY+yi1)*fkNX+xi1]*dx1*dy1
	     +fkMeanBy[(zi2*fkNY+yi2)*fkNX+xi1]*dx1*dy
             +fkMeanBy[(zi2*fkNY+yi1)*fkNX+xi2]*dx *dy1
             +fkMeanBy[(zi2*fkNY+yi2)*fkNX+xi2]*dx *dy;



  *Bx=0.5*(((snpx-2.*snx+snmx)*dz+2.*(snx-snmx))*dz+snx-s0x+snmx);
  *By=0.5*(((snpy-2.*sny+snmy)*dz+2.*(sny-snmy))*dz+sny-s0y+snmy);
  //TODO: make this nice
  if (TMath::Abs(z)>0.001) {
    *Bx/=z;
    *By/=z;
  }
}
