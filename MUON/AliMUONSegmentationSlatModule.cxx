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

/*
$Log$
*/

/////////////////////////////////////////////////////
//  Segmentation classes for slat modules          //
//  to be used with AluMUONSegmentationSlat        //
/////////////////////////////////////////////////////


#include "AliMUONSegmentationSlatModule.h"
#include <TMath.h>
#include <iostream.h>

#include "AliMUONSegmentationV01.h"

//___________________________________________
ClassImp(AliMUONSegmentationSlatModule)

AliMUONSegmentationSlatModule::AliMUONSegmentationSlatModule() 
{
// Default constructor
    fNsec=4;
    fNDiv = new TArrayI(fNsec);      
    fDpxD = new TArrayF(fNsec);      
    (*fNDiv)[0]=(*fNDiv)[1]=(*fNDiv)[2]=(*fNDiv)[3]=0;     
    (*fDpxD)[0]=(*fDpxD)[1]=(*fDpxD)[2]=(*fDpxD)[3]=0;     
}

void AliMUONSegmentationSlatModule::SetPcbBoards(Int_t n[4])
{
//
// Set Pcb Board segmentation zones
    for (Int_t i=0; i<4; i++) fPcbBoards[i]=n[i];
}


void AliMUONSegmentationSlatModule::SetPadDivision(Int_t ndiv[4])
{
//
// Defines the pad size perp. to the anode wire (y) for different sectors. 
// Pad sizes are defined as integral fractions ndiv of a basis pad size
// fDpx
// 
    for (Int_t i=0; i<4; i++) {
	(*fNDiv)[i]=ndiv[i];
    }
    ndiv[0]=ndiv[1];
}

Float_t AliMUONSegmentationSlatModule::Dpx(Int_t isec) const
{
// Return x-strip width
    return (*fDpxD)[isec];
} 


Float_t AliMUONSegmentationSlatModule::Dpy(Int_t isec) const
{
// Return y-strip width

    return fDpy;
}


void AliMUONSegmentationSlatModule::
GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy) 
{
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
    iy = Int_t(y/fDpy)+1;
    if (iy >  fNpy) iy= fNpy;
//
//  Find sector isec
    
    Int_t isec=-1;
    for (Int_t i=fNsec-1; i > 0; i--) {
	if (x >= fCx[i-1]) {
	    isec=i;
	    break;
	}
    }

    if (isec>0) {
	ix= Int_t((x-fCx[isec-1])/(*fDpxD)[isec])
	    +fNpxS[isec-1]+1;
    } else if (isec == 0) {
	ix= Int_t(x/(*fDpxD)[isec])+1;
    } else {
	ix=0;
	iy=0;
    }
}

void AliMUONSegmentationSlatModule::
GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y) 
{
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    y = Float_t(iy*fDpy)-fDpy/2.;
//
//  Find sector isec
    Int_t isec=AliMUONSegmentationSlatModule::Sector(ix,iy);
    if (isec == -1) printf("\n PadC %d %d %d  \n ", isec, ix, iy);
//
    if (isec>0) {
	x = fCx[isec-1]+(ix-fNpxS[isec-1])*(*fDpxD)[isec];
	x = x-(*fDpxD)[isec]/2;
    } else {
	x=y=0;
    }
}

void AliMUONSegmentationSlatModule::
SetPad(Int_t ix, Int_t iy)
{
    //
    // Sets virtual pad coordinates, needed for evaluating pad response 
    // outside the tracking program 
    GetPadC(ix,iy,fX,fY);
    fSector=Sector(ix,iy);
}

void AliMUONSegmentationSlatModule::
SetHit(Float_t x, Float_t y)
{
    fXhit = x;
    fYhit = y;
    
    if (x < 0) fXhit = 0;
    if (y < 0) fYhit = 0;
    
    if (x >= fCx[fNsec-1]) fXhit = fCx[fNsec-1];
    if (y >= fDyPCB)       fYhit = fDyPCB;

    
}


void AliMUONSegmentationSlatModule::
FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{
// Initialises iteration over pads for charge distribution algorithm
//
    //
    // Find the wire position (center of charge distribution)
    Float_t x0a=GetAnod(xhit);
    fXhit=x0a;
    fYhit=yhit;
    //
    // and take fNsigma*sigma around this center
    Float_t x01=x0a  - dx;
    Float_t x02=x0a  + dx;
    Float_t y01=yhit - dy;
    Float_t y02=yhit + dy;
    
    if (x01 < 0) x01 = 0;
    if (y01 < 0) y01 = 0;
    
    Int_t isec=-1;
    for (Int_t i=fNsec-1; i > 0; i--) {
	if (x02 >= fCx[i-1]) {
	    isec=i;
	    break;
	}
    }
   
    if (x02 >= fCx[fNsec-1]) x02 = fCx[fNsec-1];
    if (y02 >= fDyPCB) y02 = fDyPCB;
    //
    // find the pads over which the charge distributes
    GetPadI(x01,y01,fIxmin,fIymin);
    GetPadI(x02,y02,fIxmax,fIymax);
    if (fIxmax > fNpx) fIxmax=fNpx;
    if (fIymax > fNpyS[isec]) fIymax = fNpyS[isec];    
    fXmin=x01;
    fXmax=x02;
    fYmin=y01;
    fYmax=y02;
    
    // 
    // Set current pad to lower left corner
    if (fIxmax < fIxmin) fIxmax=fIxmin;
    if (fIymax < fIymin) fIymax=fIymin;    
    fIx=fIxmin;
    fIy=fIymin;
    
    GetPadC(fIx,fIy,fX,fY);
    fSector=Sector(fIx,fIy);
//    printf("\n \n First Pad: %d %d %f %f %d %d %d %f" , 
//	   fIxmin, fIxmax, fXmin, fXmax, fNpx, fId, isec, Dpy(isec));    
//    printf("\n \n First Pad: %d %d %f %f %d %d %d %f",
//	   fIymin, fIymax, fYmin, fYmax, fNpy, fId, isec, Dpy(isec));
}

void AliMUONSegmentationSlatModule::NextPad()
{
// Stepper for the iteration over pads
//
// Step to next pad in the integration region
//  step from left to right    
    if (fIx != fIxmax) {
	fIx++;
	GetPadC(fIx,fIy,fX,fY);
	fSector=Sector(fIx,fIy);
//  step up 
    } else if (fIy != fIymax) {
	fIx=fIxmin;
	fIy++;
	GetPadC(fIx,fIy,fX,fY);
	fSector=Sector(fIx,fIy);

    } else {
	fIx=-1;
	fIy=-1;
    }
//    printf("\n Next Pad %d %d %f %f %d %d %d %d %d ", 
}


Int_t AliMUONSegmentationSlatModule::MorePads()
// Stopping condition for the iterator over pads
//
// Are there more pads in the integration region
{
    
    return  (fIx != -1  || fIy != -1);
}


Int_t AliMUONSegmentationSlatModule::Sector(Int_t ix, Int_t iy) 
{
//
// Determine segmentation zone from pad coordinates
//
    Int_t isec=-1;
    for (Int_t i=0; i < fNsec; i++) {
	if (ix <= fNpxS[i]) {
	    isec=i;
	    break;
	}
    }
    if (isec == -1) printf("\n Sector: Attention isec ! %d %d %d %d \n",
			   fId, ix, iy,fNpxS[3]);

    return isec;

}

void AliMUONSegmentationSlatModule::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2) 
{
//  Returns integration limits for current pad
//

    x1=fXhit-fX-Dpx(fSector)/2.;
    x2=x1+Dpx(fSector);
    y1=fYhit-fY-Dpy(fSector)/2.;
    y2=y1+Dpy(fSector);    
//    printf("\n Integration Limits %f %f %f %f %d %f", x1, x2, y1, y2, fSector, Dpx(fSector));

}

void AliMUONSegmentationSlatModule::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]) 
{
// Returns list of next neighbours for given Pad (iX, iY)
//
//
    Int_t i=0;
//    
//  step right
    if (iX+1 <= fNpx) {
	Xlist[i]=iX+1;
	Ylist[i++]=iY;
    }
//
//  step left    
    if (iX-1 > 0) {
	Xlist[i]=iX-1;
	Ylist[i++]=iY;
    }

//    
//  step up
    if (iY+1 <= fNpy) {
	Xlist[i]=iX;
	Ylist[i++]=iY+1;
    }
//
//  step down    
    if (iY-1 > 0) {
	Xlist[i]=iX;
	Ylist[i++]=iY-1;
    }

    *Nlist=i;
}


void AliMUONSegmentationSlatModule::Init(Int_t chamber)
{
    printf("\n Initialise segmentation SlatModule \n");
//
//  Fill the arrays fCx (x-contour) and fNpxS (ix-contour) for each sector
//  These arrays help in converting from real to pad co-ordinates and
//  vice versa
//   
//  Segmentation is defined by rectangular modules approximating
//  concentric circles as shown below
//
//  PCB module size in cm
    fDxPCB=40;
    fDyPCB=40;
//
// number of pad rows per PCB
//    
    Int_t nPyPCB=Int_t(fDyPCB/fDpy);
//
// maximum number of pad rows    
    fNpy=nPyPCB;
//
//  Calculate padsize along x
    (*fDpxD)[fNsec-1]=fDpx;
    if (fNsec > 1) {
	for (Int_t i=fNsec-2; i>=0; i--){
	    (*fDpxD)[i]=(*fDpxD)[fNsec-1]/(*fNDiv)[i];
	    printf("\n test ---dx %d %f \n",i,(*fDpxD)[i]);
	}
    }
//
// fill the arrays defining the pad segmentation boundaries
//
//  
//  Loop over sectors (isec=0 is the dead space surounding the beam pipe)
    for (Int_t isec=0; isec<4; isec++) {
	if (isec==0) {
	    fNpxS[0] = 0;
	    fNpyS[0] = 0;
	    fCx[0]   = 0;
	} else {
	    fNpxS[isec]=fNpxS[isec-1] + fPcbBoards[isec]*Int_t(fDxPCB/(*fDpxD)[isec]);
	    fNpyS[isec]=fNpy;
	    fCx[isec]=fCx[isec-1] + fPcbBoards[isec]*fDxPCB;
	}
    } // sectors
// maximum number of pad rows    
    fNpy=nPyPCB;
    fNpx=fNpxS[3];
}












