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
Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/06/09 21:37:56  morsch
AliMUONSegmentationV02 code  from  AliMUONSegResV02.cxx

*/



/////////////////////////////////////////////////////
//  Segmentation and Response classes version 02   //
/////////////////////////////////////////////////////


#include "AliMUONSegmentationV02.h"
#include "iostream.h"

//___________________________________________
ClassImp(AliMUONSegmentationV02)

void AliMUONSegmentationV02::SetPadSize(Float_t p1, Float_t p2)
{
//  Sets the padsize 
//
    fDpy=p1;
    fDpx=p2;
}

Int_t AliMUONSegmentationV02::Npx()
// Returns maximum number if pads in x
{return AliMUONSegmentationV01::Npy();}

Int_t AliMUONSegmentationV02::Npy()
// Returns maximum number if pads in y
{return AliMUONSegmentationV01::Npx();}


Float_t AliMUONSegmentationV02::Dpx(Int_t isec)
// Returns pad-size in x
{return fDpy;}

Float_t AliMUONSegmentationV02::Dpy(Int_t isec)
// Returns pad-size in y
{return fDpxD[isec];}
Int_t AliMUONSegmentationV02::Sector(Int_t ix, Int_t iy)
// Returns sector number for given pad position
//
{return AliMUONSegmentationV01::Sector(iy, ix);}

void AliMUONSegmentationV02::

GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
{
AliMUONSegmentationV01::GetPadI(y, x, iy, ix);
// printf("\n x,y,ix,iy %f %f %d %d", x,y,ix,iy);
}

void AliMUONSegmentationV02::
GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
{
    AliMUONSegmentationV01::GetPadC(iy, ix, y, x);
}
void AliMUONSegmentationV02::SetPad(Int_t ix,Int_t iy)
{
    //
    // Sets virtual pad coordinates, needed for evaluating pad response 
    // outside the tracking program 
    GetPadC(ix,iy,fx,fy);
    fSector=Sector(ix,iy);    
}



void AliMUONSegmentationV02::NextPad()
{
// Stepper for the iteration over pads
//
  // 
  // Step to next pad in integration region
    Float_t xc,yc;
    Int_t   ixc;
    
//  step up    
    if (fy < fymax && fy != 0) {
	if (fiy==-1) fiy++;
	fiy++;
//  step right 
    } else if (fix != fixmax) {
	if (fix==-1) fix++;
	fix++;
//      get y-position of next row (yc), xc not used here 	
	GetPadC(fix,fiy,xc,yc);
//      get x-pad coordiante for 1 pad in row (fix)
	GetPadI(xc,fymin,ixc,fiy);
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadC(fix,fiy,fx,fy);
    fSector=Sector(fix,fiy);
    if (MorePads() && 
	(fSector ==-1 || fSector==0 )) 
	NextPad();
//    printf("\n this pad %f %f %d %d \n",fx,fy,fix,fiy);
    
}

Int_t AliMUONSegmentationV02::MorePads()
// Stopping condition for the iterator over pads
//
//
// Are there more pads in the integration region
{
    if ((fy >= fymax  && fix >= fixmax) || fx==0) {
	return 0;
    } else {
	return 1;
    }
}

void AliMUONSegmentationV02::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])
{
// Returns list of next neighbours for given Pad (iX, iY)
//
    Int_t n;
    AliMUONSegmentationV01::Neighbours(iY, iX, &n, Ylist, Xlist);
    *Nlist=n;
}




