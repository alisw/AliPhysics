#include "TMath.h"
#include "TTreeStream.h"
#include "AliTPCExBFirst.h"

ClassImp(AliTPCExBFirst)

const Double_t AliTPCExBFirst::fgkEM=1.602176487e-19/9.10938215e-31;
const Double_t AliTPCExBFirst::fgkDriftField=40.e3;

AliTPCExBFirst::AliTPCExBFirst(const AliMagF *bField,
			       Double_t driftVelocity,
			       Int_t nx,Int_t ny,Int_t nz)
  : fkNX(nx),fkNY(ny),fkNZ(nz),
    fkXMin(-250.),fkXMax(250.),fkYMin(-250.),fkYMax(250.),
    fkZMin(-250.),fkZMax(250.) {
  //
  // The constructor. One has to supply a magnetic field and an (initial)
  // drift velocity. Since some kind of lookuptable is created the
  // number of its meshpoints can be supplied.
  //
  fDriftVelocity=driftVelocity;
  ConstructCommon(0,bField);
}

AliTPCExBFirst::AliTPCExBFirst(const AliFieldMap *bFieldMap,
			       Double_t driftVelocity) {
  //
  // The constructor. One has to supply a field map and an (initial)
  // drift velocity.
  //
  fDriftVelocity=driftVelocity;

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

AliTPCExBFirst::~AliTPCExBFirst() { 
  //
  // destruct the poor object.
  //
  delete[] fMeanBx;
  delete[] fMeanBy;
}

void AliTPCExBFirst::Correct(const Double_t *position,Double_t *corrected) {
  //
  // correct for the distortion
  //
  Double_t r=TMath::Sqrt(position[0]*position[0]+position[1]*position[1]);
  if (TMath::Abs(position[2])>250.||r<90.||250.<r) {
    for (Int_t i=0;i<3;++i) corrected[i]=position[i];
  }
  else {
    Double_t Bx,By;
    GetMeanFields(position[0],position[1],position[2],&Bx,&By);
    if (position[2]>0.) {
      Double_t Bxe,Bye;
      GetMeanFields(position[0],position[1],250.,&Bxe,&Bye);
      if (position[2]!=250.) {
	Bx=(500.*Bxe-(position[2]+250.)*Bx)/(250.-position[2]);
	By=(500.*Bye-(position[2]+250.)*By)/(250.-position[2]);
      }
      else {
	Bx=Bxe;
	By=Bye;
      }
    }

    Double_t mu=fDriftVelocity/fgkDriftField;
    Double_t wt=mu*fMeanBz;

    corrected[0]=mu*(wt*Bx-By)/(1.+wt*wt);
    corrected[1]=mu*(wt*By+Bx)/(1.+wt*wt);

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

void AliTPCExBFirst::ConstructCommon(const AliFieldMap *bFieldMap,
				     const AliMagF *bField) {
  //
  // THIS IS PRIVATE! (a helper for the constructor)
  //
  fMeanBx=new Double_t[fkNX*fkNY*fkNZ];
  fMeanBy=new Double_t[fkNX*fkNY*fkNZ];

  Double_t x[3];
  Double_t nBz=0;
  fMeanBz=0.;
  for (int i=0;i<fkNX;++i) {
    x[0]=fkXMin+i*(fkXMax-fkXMin)/(fkNX-1);
    for (int j=0;j<fkNY;++j) {
      x[1]=fkYMin+j*(fkYMax-fkYMin)/(fkNY-1);
      Double_t R=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
      Double_t Bx=0.,By=0.;
      for (int k=0;k<fkNZ;++k) {
	x[2]=fkZMin+k*(fkZMax-fkZMin)/(fkNZ-1);
	Float_t B[3];
	// the x is not const in the Field function...
	Float_t xt[3];
	for (int l=0;l<3;++l) xt[l]=x[l];
	// that happens due to the lack of a sophisticated class design:
	if (bFieldMap!=0)
	  bFieldMap->Field(xt,B);
	else 
	  bField->Field(xt,B);
	Bx+=B[0]/10.;
	By+=B[1]/10.;
	/*old
	fMeanBx[(k*fkNY+j)*fkNX+i]=Bx/(k+1);
	fMeanBy[(k*fkNY+j)*fkNX+i]=By/(k+1);
	*/
	fMeanBx[(k*fkNY+j)*fkNX+i]=Bx;
	fMeanBy[(k*fkNY+j)*fkNX+i]=By;
	if (90.<=R&&R<=250.) {
	  fMeanBz+=B[2]/10.;
	  ++nBz;
	}
      }
    }
  }
  fMeanBz/=nBz;
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

  double s0x=fMeanBx[yi1*fkNX+xi1]*dx1*dy1
	    +fMeanBx[yi2*fkNX+xi1]*dx1*dy
            +fMeanBx[yi1*fkNX+xi2]*dx *dy
            +fMeanBx[yi2*fkNX+xi2]*dx *dy1;
  double s0y=fMeanBy[yi1*fkNX+xi1]*dx1*dy1
	    +fMeanBy[yi2*fkNX+xi1]*dx1*dy
            +fMeanBy[yi1*fkNX+xi2]*dx *dy
            +fMeanBy[yi2*fkNX+xi2]*dx *dy1;
  Int_t zi0=zi1-1;
  double snmx,snmy;
  if (zi0>=0) {
    snmx=fMeanBx[(zi0*fkNY+yi1)*fkNX+xi1]*dx1*dy1
        +fMeanBx[(zi0*fkNY+yi2)*fkNX+xi1]*dx1*dy
        +fMeanBx[(zi0*fkNY+yi1)*fkNX+xi2]*dx *dy
        +fMeanBx[(zi0*fkNY+yi2)*fkNX+xi2]*dx *dy1;
    snmy=fMeanBy[(zi0*fkNY+yi1)*fkNX+xi1]*dx1*dy1
        +fMeanBy[(zi0*fkNY+yi2)*fkNX+xi1]*dx1*dy
        +fMeanBy[(zi0*fkNY+yi1)*fkNX+xi2]*dx *dy
        +fMeanBy[(zi0*fkNY+yi2)*fkNX+xi2]*dx *dy1;
  }
  else
    snmx=snmy=0.;
  double snx=fMeanBx[(zi1*fkNY+yi1)*fkNX+xi1]*dx1*dy1
 	    +fMeanBx[(zi1*fkNY+yi2)*fkNX+xi1]*dx1*dy
            +fMeanBx[(zi1*fkNY+yi1)*fkNX+xi2]*dx *dy
            +fMeanBx[(zi1*fkNY+yi2)*fkNX+xi2]*dx *dy1;
  double sny=fMeanBy[(zi1*fkNY+yi1)*fkNX+xi1]*dx1*dy1
 	    +fMeanBy[(zi1*fkNY+yi2)*fkNX+xi1]*dx1*dy
            +fMeanBy[(zi1*fkNY+yi1)*fkNX+xi2]*dx *dy
            +fMeanBy[(zi1*fkNY+yi2)*fkNX+xi2]*dx *dy1;
  double snpx=fMeanBx[(zi2*fkNY+yi1)*fkNX+xi1]*dx1*dy1
	     +fMeanBx[(zi2*fkNY+yi2)*fkNX+xi1]*dx1*dy
             +fMeanBx[(zi2*fkNY+yi1)*fkNX+xi2]*dx *dy
             +fMeanBx[(zi2*fkNY+yi2)*fkNX+xi2]*dx *dy1;
  double snpy=fMeanBy[(zi2*fkNY+yi1)*fkNX+xi1]*dx1*dy1
	     +fMeanBy[(zi2*fkNY+yi2)*fkNX+xi1]*dx1*dy
             +fMeanBy[(zi2*fkNY+yi1)*fkNX+xi2]*dx *dy
             +fMeanBy[(zi2*fkNY+yi2)*fkNX+xi2]*dx *dy1;
  *Bx=0.5*(((snpx-2.*snx+snmx)*dz+2.*(snx-snmx))*dz+snx-s0x+snmx);
  *By=0.5*(((snpy-2.*sny+snmy)*dz+2.*(sny-snmy))*dz+sny-s0y+snmy);
  //TODO: make this nice
  if (TMath::Abs(z)>0.001) {
    *Bx/=z;
    *By/=z;
  }
}
