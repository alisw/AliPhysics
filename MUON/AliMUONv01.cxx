/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliMUONv01.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"

//___________________________________________
ClassImp(AliMUONsegmentationV01)

AliMUONsegmentationV01::AliMUONsegmentationV01() 
{
    fNsec=4;
    fRSec.Set(fNsec);    
    fNDiv.Set(fNsec);      
    fDpxD.Set(fNsec);      
}
    
void   AliMUONsegmentationV01::SetSegRadii(Float_t  r[4])
{
    for (Int_t i=0; i<4; i++) {
	fRSec[i]=r[i];
	printf("\n R %d %f",i,fRSec[i]);
	
    }
}

void AliMUONsegmentationV01::SetPadDivision(Int_t ndiv[4])
{
//
// Defines the pad size perp. to the anode wire (y) for different sectors. 
//
    for (Int_t i=0; i<4; i++) {
	fNDiv[i]=ndiv[i];
	printf("\n Ndiv %d %d",i,fNDiv[i]);
    }
    ndiv[0]=ndiv[1];
}


void AliMUONsegmentationV01::Init(AliMUONchamber*)
{
//
//  Fill the arrays fCx (x-contour) and fNpx (ix-contour) for each sector
//  These arrays help in converting from real to pad co-ordinates and
//  vice versa
//
    Int_t isec;
    printf("\n Initialise segmentation v01");
    fNpy=Int_t(fRSec[fNsec-1]/fDpy)+1;

    fDpxD[fNsec-1]=fDpx;
    if (fNsec > 1) {
	for (Int_t i=fNsec-2; i>=0; i--){
	    fDpxD[i]=fDpxD[fNsec-1]/fNDiv[i];
	    printf("\n dx %d %f",i,fDpxD[i]);
	}
    }
    fWireD=fDpxD[1]/3;
//
// fill the arrays defining the pad segmentation boundaries
    Float_t ry;
    Int_t   dnx;
//
//  loop over sections
    for(isec=0; isec<fNsec; isec++) {
//  
//  loop over pads along the aode wires
	for (Int_t iy=1; iy<=fNpy; iy++) {
//
	    Float_t x=iy*fDpy-fDpy/2;
	    if (x > fRSec[isec]) {
		fNpx[isec][iy]=0;
		 fCx[isec][iy]=0;
	    } else {
		ry=TMath::Sqrt(fRSec[isec]*fRSec[isec]-x*x);
		if (isec > 0) {
		    dnx= Int_t((ry-fCx[isec-1][iy])/fDpxD[isec]);
		    if (TMath::Odd(dnx)) dnx--;
		    fNpx[isec][iy]=fNpx[isec-1][iy]+dnx;
		    fCx[isec][iy]=fCx[isec-1][iy]+dnx*fDpxD[isec];
		} else {
		    dnx=Int_t(ry/fDpxD[isec]);
		    fNpx[isec][iy]=dnx;
		    if (TMath::Odd(dnx)) dnx--;
		    fCx[isec][iy]=dnx*fDpxD[isec];
		}
	    }
	} // y-pad loop
    } // sector loop
 
    //   
    // for debugging only 
    for (Int_t iy=0; iy<fNpy; iy++) {
	printf("\n iy %d",iy);
	for(isec=0; isec<fNsec; isec++) {
	    printf("  %d", 
		   fNpx[isec][iy]);
	}
	printf("\n iy %d",iy);
	for(isec=0; isec<fNsec; isec++) {
	    printf("  %f", 
		   fCx[isec][iy] );
	}
	printf("\n");
    }
    Float_t x,y;
    Int_t jx,jy;
    //   
    // for debugging only 
    for (Int_t ix=1; ix<100; ix++) {
	for (Int_t iy=1; iy<100; iy++) {
	    GetPadCxy(ix,iy,x,y);
	    GetPadIxy(x,y,jx,jy);
	    if ((ix != jx) && (jx!=-1)) {
		printf("\n %d %f %d %f", ix,x,iy,y);
		printf("\n %d %f %d %f", jx,x,jy,y);
		printf("\n \n **********");
	    }
	}
    }
}

Int_t AliMUONsegmentationV01::Sector(Int_t ix, Int_t iy)
{
    Int_t absix=TMath::Abs(ix);
    Int_t absiy=TMath::Abs(iy);
    Int_t isec=0;
    for (Int_t i=0; i<fNsec; i++) {
	if (absix<=fNpx[i][absiy]){
	    isec=i;
	    break;
	}
    }
    return isec;
}


Float_t AliMUONsegmentationV01::GetAnod(Float_t xhit)
{
//
// Finds anod-wire position closest to xhit
    Float_t wire= (xhit<0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
//    printf("getanod: %f %f %f ",xhit, wire, fWireD*wire);
    return fWireD*wire;
}

    void AliMUONsegmentationV01::
    GetPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
    iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy)-1;
    if (iy >  fNpy) iy= fNpy;
    if (iy < -fNpy) iy=-fNpy;
//
//  Find sector isec
    Int_t isec=-1;
    Float_t absx=TMath::Abs(x);
    Int_t absiy=TMath::Abs(iy);
    for (Int_t i=0; isec < fNsec; i++) {
	if (absx <= fCx[i][absiy]) {
	    isec=i;
	    break;
	}
    }
    if (isec>0) {
	ix= Int_t((absx-fCx[isec-1][absiy])/fDpxD[isec])
	    +fNpx[isec-1][absiy]+1;
    } else if (isec == 0) {
	ix= Int_t(absx/fDpxD[isec])+1;
    } else {
	ix=fNpx[fNsec-1][absiy]+1;	
    }
//    printf("\n %d %d",isec,absiy);
    
    ix = (x>0) ? ix:-ix;
}

void AliMUONsegmentationV01::
GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)+fDpy/2.;
//
//  Find sector isec
    Int_t isec=Sector(ix,iy);
//
    Int_t absix=TMath::Abs(ix);
    Int_t absiy=TMath::Abs(iy);
    if (isec) {
	x=fCx[isec-1][absiy]+(absix-fNpx[isec-1][absiy])*fDpxD[isec];
	x=(ix>0) ?  x-fDpxD[isec]/2 : -x+fDpxD[isec]/2;
    } else {
	x=y=0;
    }
}


void AliMUONsegmentationV01::
GetSuperPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
    ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx)-1;
    iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy)-1;
}

void AliMUONsegmentationV01::
GetSuperPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y) 
{
    x = (ix>0) ? Float_t(ix*fDpx)-fDpx/2. : Float_t(ix*fDpx)+fDpx/2.;
    y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)+fDpy/2.;
}

void AliMUONsegmentationV01::SetPADSIZ(Float_t p1, Float_t p2)
{
    fDpx=p1;
    fDpy=p2;
}

void AliMUONsegmentationV01::
FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{
    //
    // Find the wire position (center of charge distribution)
    Float_t x0a=GetAnod(xhit);
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
    fxmin=x01;
    fxmax=x02;
    //
    // upper and lower limits should be checked
    /*
    printf("\n FirstPad called");
    printf("\n Hit Position %f %f",xhit,yhit);
    printf("\n Closest wire %f", x0a);
    printf("\n Integration limits: %i %i %i %i",fixmin,fixmax,fiymin,fiymax);
    printf("\n Integration limits: %f %f %f %f \n",x01,x02,y01,y02);
    */
    // 
    // Set current pad to lower left corner
    if (fixmax < fixmin) fixmax=fixmin;
    if (fiymax < fiymin) fiymax=fiymin;    
    fix=fixmin;
    fiy=fiymin;
    GetPadCxy(fix,fiy,fx,fy);
    
}

void AliMUONsegmentationV01::NextPad()
{
  // 
  // Step to next pad in integration region
    Float_t xc,yc;
    Int_t   iyc;
    
//  step from left to right    
    if (fx < fxmax && fx != 0) {
	fix++;
//  step up 
    } else if (fiy != fiymax) {
	fiy++;
//      get y-position of next row (yc), xc not used here 	
	GetPadCxy(fix,fiy,xc,yc);
//      get x-pad coordiante for 1 pad in row (fix)
	GetPadIxy(fxmin,yc,fix,iyc);
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadCxy(fix,fiy,fx,fy);
    fSector=Sector(fix,fiy);
//    printf("\n this pad %f %f %d %d",fx,fy,fix,fiy);
    
}

Int_t AliMUONsegmentationV01::MorePads()
//
// Are there more pads in the integration region
{
    if ((fx >= fxmax  && fiy >= fiymax) || fy==0) {
	return 0;
    } else {
	return 1;
    }
}

void AliMUONsegmentationV01::SigGenInit(Float_t x,Float_t y,Float_t)
{
//
//  Initialises pad and wire position during stepping
    fxt=x;
    fyt=y;
    GetPadIxy(x,y,fixt,fiyt);
    fiwt= (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}

Int_t AliMUONsegmentationV01::SigGenCond(Float_t x,Float_t y,Float_t)
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

void AliMUONsegmentationV01::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
    x1=fxt-fx-fDpxD[fSector]/2.;
    x2=x1+fDpxD[fSector];
    y1=fyt-fy-fDpy/2.;
    y2=y1+fDpy;    
}

void AliMUONsegmentationV01::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])
{
    const Float_t epsilon=fDpy/1000;
    
    Float_t x,y;
    Int_t   ixx, iyy, isec1;
//
    Int_t   isec0=Sector(iX,iY);
    Int_t i=0;
//    
//  step right
    Xlist[i]=iX+1;
    Ylist[i++]=iY;
//
//  step left    
    Xlist[i]=iX-1;
    Ylist[i++]=iY;
//
//  step up 
    GetPadCxy(iX,iY,x,y);
    GetPadIxy(x+epsilon,y+fDpy,ixx,iyy);
    Xlist[i]=ixx;
    Ylist[i++]=iY+1;
    isec1=Sector(ixx,iyy);
    if (isec1==isec0) {
//
//  no sector boundary crossing
    Xlist[i]=ixx+1;
    Ylist[i++]=iY+1;
	
    Xlist[i]=ixx-1;
    Ylist[i++]=iY+1;
    } else if (isec1 < isec0) {
// finer segmentation
	Xlist[i]=ixx+1;
	Ylist[i++]=iY+1;
	
	Xlist[i]=ixx-1;
	Ylist[i++]=iY+1;
	
	Xlist[i]=ixx-2;
	Ylist[i++]=iY+1;
    } else {
// coarser segmenation	
	if (TMath::Odd(iX-fNpx[isec1-1][iY+1])) {
	    Xlist[i]=ixx-1;
	    Ylist[i++]=iY+1;
	} else {
	    Xlist[i]=ixx+1;
	    Ylist[i++]=iY+1;
	}
    }
//
// step down 
    GetPadCxy(iX,iY,x,y);
    GetPadIxy(x+epsilon,y-fDpy,ixx,iyy);
    Xlist[i]=ixx;
    Ylist[i++]=iY-1;
    isec1=Sector(ixx,iyy);
    if (isec1==isec0) {
//
//  no sector boundary crossing
    Xlist[i]=ixx+1;
    Ylist[i++]=iY-1;
	
    Xlist[i]=ixx-1;
    Ylist[i++]=iY-1;
    } else if (isec1 < isec0) {
// finer segmentation
	Xlist[i]=ixx+1;
	Ylist[i++]=iY-1;
	
	Xlist[i]=ixx-1;
	Ylist[i++]=iY-1;
	
	Xlist[i]=ixx-2;
	Ylist[i++]=iY-1;
    } else {
// coarser segmenation	
	if (TMath::Odd(iX-fNpx[isec1-1][iY-1])) {
	    Xlist[i]=ixx-1;
	    Ylist[i++]=iY-1;
	} else {
	    Xlist[i]=ixx+1;
	    Ylist[i++]=iY-1;
	}
    }
    *Nlist=i;
}

//___________________________________________
void AliMUONsegmentationV01::
FitXY(AliMUONRecCluster* Cluster,TClonesArray* MUONdigits)
    // Default : Centre of gravity method
{
    Float_t x=0;
    Float_t y=0;
    Float_t q=0;
    Float_t xToAdd;
    Float_t yToAdd;
    
    if (gAlice->TreeD()->GetReadEvent() != Cluster->GetCathod()+1)
	// next line warns if in the future cathod 1 is not event 2 !
	printf("ClusterFillXY : not reading the right cathod !\n");
    for(Int_t clusDigit=Cluster->FirstDigitIndex();
	clusDigit!=Cluster->InvalidDigitIndex();
	clusDigit=Cluster->NextDigitIndex()) {
	AliMUONdigit* pDigit=(AliMUONdigit*)MUONdigits->UncheckedAt(clusDigit);
	GetPadCxy(pDigit->fPadX,pDigit->fPadY,xToAdd,yToAdd);
	x+= xToAdd*pDigit->fSignal;
	y+= yToAdd*pDigit->fSignal;
	q+= (Float_t) pDigit->fSignal;
    }
    Cluster->fX=x/q;
    Cluster->fY=y/q;
}








