/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TBox.h> 
#include <TRandom.h> 

#include "AliMUONSegResV01.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"

//___________________________________________
ClassImp(AliMUONsegmentationV01)

Float_t AliMUONsegmentationV01::Dpx(Int_t isec)
{
   return fDpxD[isec];
}

Float_t AliMUONsegmentationV01::Dpy(Int_t )
{
   return fDpy;
}

AliMUONsegmentationV01::AliMUONsegmentationV01() 
{
    fNsec=4;
    fRSec.Set(fNsec);    
    fNDiv.Set(fNsec);      
    fDpxD.Set(fNsec);      
    fRSec[0]=fRSec[1]=fRSec[2]=fRSec[3]=0;     
    fNDiv[0]=fNDiv[1]=fNDiv[2]=fNDiv[3]=0;     
    fDpxD[0]=fDpxD[1]=fDpxD[2]=fDpxD[3]=0;     
    fCorr = new TObjArray(3);
    (*fCorr)[0]=0;
    (*fCorr)[1]=0;
    (*fCorr)[2]=0;
} 
    
void   AliMUONsegmentationV01::SetSegRadii(Float_t  r[4])
{
    for (Int_t i=0; i<4; i++) {
	fRSec[i]=r[i];
	//printf("\n R %d %f \n",i,fRSec[i]);
	
    }
}


void AliMUONsegmentationV01::SetPadDivision(Int_t ndiv[4])
{
//
// Defines the pad size perp. to the anode wire (y) for different sectors. 
//
    for (Int_t i=0; i<4; i++) {
	fNDiv[i]=ndiv[i];
	//printf("\n Ndiv %d %d \n",i,fNDiv[i]);
    }
    ndiv[0]=ndiv[1];
}


void AliMUONsegmentationV01::Init(AliMUONchamber* )
{
//
//  Fill the arrays fCx (x-contour) and fNpxS (ix-contour) for each sector
//  These arrays help in converting from real to pad co-ordinates and
//  vice versa
//
    Int_t isec;
    //printf("\n Initialise segmentation v01 -- test !!!!!!!!!!!!!! \n");
    fNpy=Int_t(fRSec[fNsec-1]/fDpy)+1;

    fDpxD[fNsec-1]=fDpx;
    if (fNsec > 1) {
	for (Int_t i=fNsec-2; i>=0; i--){
	    fDpxD[i]=fDpxD[fNsec-1]/fNDiv[i];
	    //printf("\n test ---dx %d %f \n",i,fDpxD[i]);
	}
    }
//
// fill the arrays defining the pad segmentation boundaries
    Float_t ry;
    Int_t   dnx;
    Int_t   add;
//
//  loop over sections
    for(isec=0; isec<fNsec; isec++) {
//  
//  loop over pads along the aode wires
	for (Int_t iy=1; iy<=fNpy; iy++) {
//
	    Float_t x=iy*fDpy-fDpy/2;
	    if (x > fRSec[isec]) {
		fNpxS[isec][iy]=0;
		 fCx[isec][iy]=0;
	    } else {
		ry=TMath::Sqrt(fRSec[isec]*fRSec[isec]-x*x);
		if (isec > 1) {
		    dnx= Int_t((ry-fCx[isec-1][iy])/fDpxD[isec]);
		    if (isec < fNsec-1) {
			if (TMath::Odd((Long_t)dnx)) dnx++;    		
		    }
                    fNpxS[isec][iy]=fNpxS[isec-1][iy]+dnx;
          	    fCx[isec][iy]=fCx[isec-1][iy]+dnx*fDpxD[isec];
		} else if (isec == 1) {
		    dnx= Int_t((ry-fCx[isec-1][iy])/fDpxD[isec]);
		    fNpxS[isec][iy]=fNpxS[isec-1][iy]+dnx;
                    add=4 - (fNpxS[isec][iy])%4;
                    if (add < 4) fNpxS[isec][iy]+=add; 
		    dnx=fNpxS[isec][iy]-fNpxS[isec-1][iy];
		    fCx[isec][iy]=fCx[isec-1][iy]+dnx*fDpxD[isec];
		} else {
		    dnx=Int_t(ry/fDpxD[isec]);
                    fNpxS[isec][iy]=dnx;
		    fCx[isec][iy]=dnx*fDpxD[isec];
		}
	    }
	} // y-pad loop
    } // sector loop
    //   
    // for debugging only 
   
    //printf("segmentationv01 - I was here ! \n");
}

Int_t AliMUONsegmentationV01::Sector(Int_t ix, Int_t iy)
{
    Int_t absix=TMath::Abs(ix);
    Int_t absiy=TMath::Abs(iy);
    Int_t isec=0;
    for (Int_t i=0; i<fNsec; i++) {
	if (absix<=fNpxS[i][absiy]){
	    isec=i;
	    break;
	}
    }
    return isec;
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
    for (Int_t i=0; i < fNsec; i++) {
	if (absx <= fCx[i][absiy]) {
	    isec=i;
	    break;
	}
    }
    if (isec>0) {
	ix= Int_t((absx-fCx[isec-1][absiy])/fDpxD[isec])
	    +fNpxS[isec-1][absiy]+1;
    } else if (isec == 0) {
	ix= Int_t(absx/fDpxD[isec])+1;
    } else {
	ix=fNpxS[fNsec-1][absiy]+1;	
    }
//    printf("\n something %d %d \n",isec,absiy);
    
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
    Int_t isec=AliMUONsegmentationV01::Sector(ix,iy);
//
    Int_t absix=TMath::Abs(ix);
    Int_t absiy=TMath::Abs(iy);
    if (isec) {
	x=fCx[isec-1][absiy]+(absix-fNpxS[isec-1][absiy])*fDpxD[isec];
	x=(ix>0) ?  x-fDpxD[isec]/2 : -x+fDpxD[isec]/2;
    } else {
	x=y=0;
    }
}

void AliMUONsegmentationV01::
SetPad(Int_t ix, Int_t iy)
{
    GetPadCxy(ix,iy,fx,fy);
    fSector=Sector(ix,iy);
}


void AliMUONsegmentationV01::
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
    fxmin=x01;
    fxmax=x02;
    fymin=y01;
    fymax=y02;
    
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
	if (fix==-1) fix++;
	fix++;
//  step up 
    } else if (fiy != fiymax) {
	if (fiy==-1) fiy++;
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
    if (MorePads() && 
	(fSector ==-1 || fSector==0 || 
	 TMath::Abs(fx)<1.5 || TMath::Abs(fy)<1.5)) 
	NextPad();
    
//    printf("\n this pad %f %f %d %d \n",fx,fy,fix,fiy);
    
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

void AliMUONsegmentationV01::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
    x1=fxhit-fx-Dpx(fSector)/2.;
    x2=x1+Dpx(fSector);
    y1=fyhit-fy-Dpy(fSector)/2.;
    y2=y1+Dpy(fSector);    
}

void AliMUONsegmentationV01::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])
{
    const Float_t epsilon=fDpy/1000;
    
    Float_t x,y;
    Int_t   ixx, iyy, isec1;
//
    Int_t   isec0=AliMUONsegmentationV01::Sector(iX,iY);
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
    AliMUONsegmentationV01::GetPadCxy(iX,iY,x,y);
    AliMUONsegmentationV01::GetPadIxy(x+epsilon,y+fDpy,ixx,iyy);
    Xlist[i]=ixx;
    Ylist[i++]=iY+1;
    isec1=AliMUONsegmentationV01::Sector(ixx,iyy);
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

	if (TMath::Odd(iX-fNpxS[isec1-1][iY+1])) {
	    Xlist[i]=ixx-1;
	    Ylist[i++]=iY+1;
	} else {
	    Xlist[i]=ixx+1;
	    Ylist[i++]=iY+1;
	}
    }
//
// step down 
    AliMUONsegmentationV01::GetPadCxy(iX,iY,x,y);
    AliMUONsegmentationV01::GetPadIxy(x+epsilon,y-fDpy,ixx,iyy);
    Xlist[i]=ixx;
    Ylist[i++]=iY-1;
    isec1=AliMUONsegmentationV01::Sector(ixx,iyy);
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
// coarser segmentation	

	if (TMath::Odd(iX-fNpxS[isec1-1][iY-1])) {
	    Xlist[i]=ixx-1;
	    Ylist[i++]=iY-1;
	} else {
	    Xlist[i]=ixx+1;
	    Ylist[i++]=iY-1;
	}

    }
    *Nlist=i;
}

void AliMUONsegmentationV01::GiveTestPoints(Int_t &n, Float_t *x, Float_t *y)
{
    n=3;
    x[0]=(fRSec[0]+fRSec[1])/2/TMath::Sqrt(2.);
    y[0]=x[0];
    x[1]=(fRSec[1]+fRSec[2])/2/TMath::Sqrt(2.);
    y[1]=x[1];
    x[2]=(fRSec[2]+fRSec[3])/2/TMath::Sqrt(2.);
    y[2]=x[2];
}

void AliMUONsegmentationV01::Draw(Option_t *)
{
    TBox *box;
    
    Float_t dx=0.95/fCx[3][1]/2;
    Float_t dy=0.95/(Float_t(Npy()))/2;
    Float_t x0,y0,x1,y1;
    Float_t xc=0.5;
    Float_t yc=0.5;
    
    for (Int_t iy=1; iy<Npy(); iy++)
    {
	for (Int_t isec=0; isec<4; isec++) {
	    if (isec==0) {
		x0=0;
		x1=fCx[isec][iy]*dx;
	    } else {
		x0=fCx[isec-1][iy]*dx;
		x1=fCx[isec][iy]*dx;
	    }
	    y0=Float_t(iy-1)*dy;
	    y1=y0+dy;
	    box=new TBox(x0+xc,y0+yc,x1+xc,y1+yc);
	    box->SetFillColor(isec+1);
	    box->Draw();

	    box=new TBox(-x1+xc,y0+yc,-x0+xc,y1+yc);
	    box->SetFillColor(isec+1);
	    box->Draw();

	    box=new TBox(x0+xc,-y1+yc,x1+xc,-y0+yc);
	    box->SetFillColor(isec+1);
	    box->Draw();

	    box=new TBox(-x1+xc,-y1+yc,-x0+xc,-y0+yc);
	    box->SetFillColor(isec+1);
	    box->Draw();
	}
    }
}
void AliMUONsegmentationV01::SetCorrFunc(Int_t isec, TF1* func)
{
    (*fCorr)[isec]=func;
    
}

TF1* AliMUONsegmentationV01::CorrFunc(Int_t isec)
{ 
    return (TF1*) (*fCorr)[isec];
}

