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


#include "AliMUONSegmentationSlatModuleN.h"
#include <TMath.h>
#include <iostream.h>

#include "AliMUONSegmentationV01.h"

//___________________________________________
ClassImp(AliMUONSegmentationSlatModuleN)

AliMUONSegmentationSlatModuleN::AliMUONSegmentationSlatModuleN() 
{
// Default constructor
}


Float_t AliMUONSegmentationSlatModuleN::Dpx(Int_t isec) const
{
//
// Returns x-pad size for given sector isec
    return fDpx;
}

Float_t AliMUONSegmentationSlatModuleN::Dpy(Int_t isec) const
{
//
// Returns y-pad size for given sector isec
    return (*fDpxD)[isec];
} 


void AliMUONSegmentationSlatModuleN::
GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
    ix = Int_t(x/fDpx)+1;
    if (ix >  fNpx) ix= fNpx;
//
//  Find sector isec
    Int_t isec=-1;
    for (Int_t i = fNsec-1; i > 0; i--) {
	if (x >= fCx[i-1]) {
	    isec=i;
	    break;
	}
    }
//
//
    if (isec == -1) {
	ix = 0;
	iy = 0;
    } else {
	iy = Int_t(y/(*fDpxD)[isec])+1;
    }
}

void AliMUONSegmentationSlatModuleN::
GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    x = Float_t(ix*fDpx)-fDpx/2.;
//
//  Find sector isec
    Int_t isec=Sector(ix,iy);
    if (isec == -1) printf("\n  gtpadc2 Warning isec !\n");  
    y = iy*(*fDpxD)[isec]-(*fDpxD)[isec]/2.;
}


void AliMUONSegmentationSlatModuleN::NextPad()
{
// Stepper for the iteration over pads
//
    Float_t xc,yc;
    Int_t   ixc;
 //  step up    
    if ((fY + Dpy(fSector)) < fYmax) {
	fIy++;
	GetPadC(fIx,fIy,fX,fY);
//  step right 
    } else if (fIx != fIxmax) {
	fIx++;
//      get y-position of next row (yc), xc not used here 	
	GetPadC(fIx,fIy,xc,yc);
//      get x-pad coordiante for 1 pad in row (fIx)
	GetPadI(xc,fYmin,ixc,fIy);
	GetPadC(fIx,fIy,fX,fY);
	fSector=Sector(fIx,fIy);
    } else {
	fIx=fIy=-1;
    }

    if (fIy > fNpyS[fSector])  printf("\n this pad %f %f %d %d \n",fX,fY,fIx,fIy);
    GetPadC(fIx, fIy, xc, yc);
//    printf("\n Next Pad (n)%d %d %f %f %d", fIx,fIy,fX,fY,fSector);
}

Int_t AliMUONSegmentationSlatModuleN::MorePads()
// Stopping condition for the iterator over pads
//
//
// Are there more pads in the integration region
{
    if ((fY >= fYmax  && fIx >= fIxmax) || fIy == -1) {
	return 0;
    } else {
	return 1;
    }
}

void AliMUONSegmentationSlatModuleN::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])
{

    
// Returns list of next neighbours for given Pad (iX, iY)
//
//
    Int_t i=0;
    Float_t x,y;
    Int_t ix,iy,isec1,isec2;
    Float_t kEpsil= 0.001;
    
//    
//  step up
    Int_t isec=Sector(iX, iY);
    
    if (iY+1 <= fNpyS[isec]) {
	Xlist[i]=iX;
	Ylist[i++]=iY+1;
    }
//
//  step down    
    if (iY-1 > 0) {
	Xlist[i]=iX;
	Ylist[i++]=iY-1;
    }
//
//    
//  step right
    
    if (iX+1 <= fNpx) {

	
	GetPadC(iX, iY, x, y);
	GetPadI(x+fDpx, y, ix, iy);
	Xlist[i]=iX+1;
	Ylist[i++]=iy;
    }
//
//  step left    
    if (iX-1 > 0) {
	isec1=Sector(iX,   iY);
	isec2=Sector(iX-1, iY);
	if (isec1==isec2) {
	    Xlist[i]=iX-1;
	    Ylist[i++]=iY;
	} else {
	    GetPadC(iX, iY, x, y);
	    GetPadI(x-fDpx, y+kEpsil, ix, iy);
	    if (ix != -1) {
		Xlist[i]=iX-1;
		Ylist[i++]=iy;
	    }
	    GetPadI(x-fDpx, y-kEpsil, ix, iy);
	    if (ix != -1) {
		Xlist[i]=iX-1;
		Ylist[i++]=iy;
	    }
	}
    }
    *Nlist=i;
}


void AliMUONSegmentationSlatModuleN::Init(Int_t chamber)
{
    printf("\n Initialise segmentation SlatModuleN \n");
//
//  Fill the arrays fCx (x-contour) for each sector
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
    fNpxPCB = Int_t(fDxPCB/fDpx) ;
//
//  Calculate padsize along y
    (*fDpxD)[fNsec-1]=fDpy;
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
	    fCx[0]   = 0;
	    fNpxS[0] = 0;
	    fNpyS[0] = 0;
	} else {
	    fNpxS[isec] = fNpxS[isec-1] + fPcbBoards[isec]*fNpxPCB;
	    fNpyS[isec] = Int_t(fDyPCB/fDpy)*(*fNDiv)[isec];

	    printf("\n %d %d ",isec, fNpxS[isec]);
	    
	    fCx[isec] = fCx[isec-1] + fPcbBoards[isec]*fDxPCB;
	    fNpx += fPcbBoards[isec] * fNpxPCB;
	}
    } // sectors

    fNpx=fNpxS[3];
    fNpy=Int_t(fDyPCB/fDpy)*(*fNDiv)[1];
}












