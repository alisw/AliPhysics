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
//  Segmentation and Response classes version 02   //
/////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include "AliMUONSegResV02.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"

//___________________________________________
ClassImp(AliMUONsegmentationV02)

void AliMUONsegmentationV02::SetPADSIZ(Float_t p1, Float_t p2)
{
    fDpy=p1;
    fDpx=p2;
}

Int_t AliMUONsegmentationV02::Npx()
{return AliMUONsegmentationV01::Npy();}

Int_t AliMUONsegmentationV02::Npy()
{return AliMUONsegmentationV01::Npx();}


Float_t AliMUONsegmentationV02::Dpx(Int_t)
{return fDpy;}

Float_t AliMUONsegmentationV02::Dpy(Int_t isec)
{return fDpxD[isec];}

Int_t AliMUONsegmentationV02::Sector(Int_t ix, Int_t iy)
{return AliMUONsegmentationV01::Sector(iy, ix);}

void AliMUONsegmentationV02::
GetPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
AliMUONsegmentationV01::GetPadIxy(y, x, iy, ix);
// printf("\n x,y,ix,iy %f %f %d %d", x,y,ix,iy);
}

void AliMUONsegmentationV02::
GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
    AliMUONsegmentationV01::GetPadCxy(iy, ix, y, x);
//    printf("\n ix,iy,x,y %d %d %f %f ", ix,iy,x,y);
}
void    AliMUONsegmentationV02::SetPad(Int_t ix,Int_t iy)
{
    GetPadCxy(ix,iy,fx,fy);
    fSector=Sector(ix,iy);    
}



void AliMUONsegmentationV02::NextPad()
{
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
	GetPadCxy(fix,fiy,xc,yc);
//      get x-pad coordiante for 1 pad in row (fix)
	GetPadIxy(xc,fymin,ixc,fiy);
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

Int_t AliMUONsegmentationV02::MorePads()
//
// Are there more pads in the integration region
{
    if ((fy >= fymax  && fix >= fixmax) || fx==0) {
	return 0;
    } else {
	return 1;
    }
}

void AliMUONsegmentationV02::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])
{
    Int_t N;
    AliMUONsegmentationV01::Neighbours(iY, iX, &N, Ylist, Xlist);
    *Nlist=N;
}


