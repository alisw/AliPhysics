#include "AliMUONSegResV0.h"
#include "TMath.h"
#include "TRandom.h"
#include "TArc.h"
#include "AliMUONchamber.h"
ClassImp(AliMUONsegmentationV0)
    void AliMUONsegmentationV0::Init(AliMUONchamber* Chamber)
{
    fNpx=(Int_t) (Chamber->ROuter()/fDpx+1);
    fNpy=(Int_t) (Chamber->ROuter()/fDpy+1);
    fRmin=Chamber->RInner();
    fRmax=Chamber->ROuter();    
    fCorr=0;
    
}


Float_t AliMUONsegmentationV0::GetAnod(Float_t xhit)
{
    Float_t wire= (xhit>0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
    return fWireD*wire;
}

void AliMUONsegmentationV0::SetPADSIZ(Float_t p1, Float_t p2)
{
    fDpx=p1;
    fDpy=p2;
}
void AliMUONsegmentationV0::
    GetPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
    ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx)-1;
    iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy)-1;
    if (iy >  fNpy) iy= fNpy;
    if (iy < -fNpy) iy=-fNpy;
    if (ix >  fNpx) ix= fNpx;
    if (ix < -fNpx) ix=-fNpx;
}
void AliMUONsegmentationV0::
GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    x = (ix>0) ? Float_t(ix*fDpx)-fDpx/2. : Float_t(ix*fDpx)+fDpx/2.;
    y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)+fDpy/2.;
}

void AliMUONsegmentationV0::
SetHit(Float_t xhit, Float_t yhit)
{
    //
    // Find the wire position (center of charge distribution)
//    Float_t x0a=GetAnod(xhit);
    fxhit=xhit;
    fyhit=yhit;
}

void AliMUONsegmentationV0::
SetPad(Int_t ix, Int_t iy)
{
    GetPadCxy(ix,iy,fx,fy);
}

void AliMUONsegmentationV0::
FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{
    //
    // Find the wire position (center of charge distribution)
    Float_t x0a=GetAnod(xhit);
    fxhit=x0a;
    fyhit=yhit;
    //
    // and take fNsigma*sigma around this center
    Float_t x01=x0a  - dx;
    Float_t x02=x0a  + dx;
    Float_t y01=yhit - dy;
    Float_t y02=yhit + dy;
    //
    // find the pads over which the charge distributes
    GetPadIxy(x01,y01,fixmin,fiymin);
    GetPadIxy(x02,y02,fixmax,fiymax);    
//    printf("\n %f %f %d %d \n",x02,y02,fixmax,fiymax);
//    printf("\n FirstPad called %f %f \n", fDpx, fDpy);
//    printf("\n Hit Position %f %f \n",xhit,yhit);
//    printf("\n Integration limits: %i %i %i %i",fixmin,fixmax,fiymin,fiymax);
//    printf("\n Integration limits: %f %f %f %f \n",x01,x02,y01,y02);
    // 
    // Set current pad to lower left corner
    fix=fixmin;
    fiy=fiymin;
    GetPadCxy(fix,fiy,fx,fy);
}

void AliMUONsegmentationV0::NextPad()
{
  // 
  // Step to next pad in integration region
    if (fix != fixmax) {
	if (fix==-1) fix++;
	fix++;
    } else if (fiy != fiymax) {
	fix=fixmin;
	if (fiy==-1) fiy++;
	fiy++;
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadCxy(fix,fiy,fx,fy);
}

Int_t AliMUONsegmentationV0::MorePads()
//
// Are there more pads in the integration region
{
    if (fix == fixmax && fiy == fiymax) {
	return 0;
    } else {
	return 1;
	
    }
}

void AliMUONsegmentationV0::SigGenInit(Float_t x,Float_t y,Float_t z)
{
//
//  Initialises pad and wire position during stepping
    fxt =x;
    fyt =y;
    GetPadIxy(x,y,fixt,fiyt);
    fiwt= (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}

Int_t AliMUONsegmentationV0::SigGenCond(Float_t x,Float_t y,Float_t z)
{
//
//  Signal will be generated if particle crosses pad boundary or
//  boundary between two wires. 
    Int_t ixt, iyt;
    GetPadIxy(x,y,ixt,iyt);
    Int_t iwt=(x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1;
    if ((ixt != fixt) || (iyt !=fiyt) || (iwt != fiwt)) {
	return 1;
    } else {
	return 0;
    }
}
void AliMUONsegmentationV0::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
//    x1=GetAnod(fxt)-fx-fDpx/2.;
    x1=fxhit-fx-fDpx/2.;
    x2=x1+fDpx;
    y1=fyhit-fy-fDpy/2.;
    y2=y1+fDpy;    
}

void AliMUONsegmentationV0::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])
{
  /*
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
  */
    *Nlist=8;
    Xlist[0]=Xlist[1]=iX;
    Xlist[2]=iX-1;
    Xlist[3]=iX+1;
    Ylist[0]=iY-1;
    Ylist[1]=iY+1;
    Ylist[2]=Ylist[3]=iY;

   // Diagonal elements
    Xlist[4]=iX+1;
    Ylist[4]=iY+1;

    Xlist[5]=iX-1;
    Ylist[5]=iY-1;

    Xlist[6]=iX-1;
    Ylist[6]=iY+1;

    Xlist[7]=iX+1;
    Ylist[7]=iY-1;
}

Float_t AliMUONsegmentationV0::Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y
, Int_t *dummy)
// Returns the square of the distance between 1 pad
// labelled by its Channel numbers and a coordinate
{
  Float_t x,y;
  GetPadCxy(iX,iY,x,y);
  return (x-X)*(x-X) + (y-Y)*(y-Y);
}


void AliMUONsegmentationV0::
FitXY(AliMUONRecCluster* Cluster,TClonesArray* MUONdigits)
    // Default : Centre of gravity method
{
    ;
}

void  AliMUONsegmentationV0::GiveTestPoints(Int_t &n, Float_t *x, Float_t *y)
{
    n=1;
    x[0]=(fRmax+fRmin)/2/TMath::Sqrt(2.);
    y[0]=x[0];
}

void AliMUONsegmentationV0::Draw()
{
    TArc *circle;
    Float_t scale=0.95/fRmax/2.;
    

    circle = new TArc(0.5,0.5,fRmax*scale,0.,360.);
    circle->SetFillColor(2);
    circle->Draw();

    circle = new TArc(0.5,0.5,fRmin*scale,0.,360.);
    circle->SetFillColor(1);
    circle->Draw();
}



//___________________________________________
ClassImp(AliMUONresponseV0)	
Float_t AliMUONresponseV0::IntPH(Float_t eloss)
{
  // Get number of electrons and return charge
     
  Int_t nel;
  nel= Int_t(eloss*1.e9/32.);
  Float_t charge=0;
  if (nel == 0) nel=1;
  for (Int_t i=1;i<=nel;i++) {
    charge -= fChargeSlope*TMath::Log(gRandom->Rndm());    
  }
  return charge;
}
// -------------------------------------------

Float_t AliMUONresponseV0::IntXY(AliMUONsegmentation * segmentation)
{

    const Float_t invpitch = 1/fPitch;
//
//  Integration limits defined by segmentation model
//  
    Float_t xi1, xi2, yi1, yi2;
    segmentation->IntegrationLimits(xi1,xi2,yi1,yi2);
    xi1=xi1*invpitch;
    xi2=xi2*invpitch;
    yi1=yi1*invpitch;
    yi2=yi2*invpitch;
//
// The Mathieson function 
    Double_t ux1=fSqrtKx3*TMath::TanH(fKx2*xi1);
    Double_t ux2=fSqrtKx3*TMath::TanH(fKx2*xi2);

    Double_t uy1=fSqrtKy3*TMath::TanH(fKy2*yi1);
    Double_t uy2=fSqrtKy3*TMath::TanH(fKy2*yi2);

    
    return Float_t(4.*fKx4*(TMath::ATan(ux2)-TMath::ATan(ux1))*
		      fKy4*(TMath::ATan(uy2)-TMath::ATan(uy1)));
}









