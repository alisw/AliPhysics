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


#include "AliRICHSegmentationV0.h"


ClassImp(AliRICHSegmentationV0)

//________________________________________________________________________________
AliRICHSegmentationV0::AliRICHSegmentationV0()
{ 

// Default constructor for AliRICHSegmentationV0 (no dead zones)

  fNpx=144;
  fNpy=160;
  fDeadZone=0;
  fSector=-1;
  fCorr = 0;
}

void AliRICHSegmentationV0::Init(Int_t id)
{

// Initialisation of chambers
  id++;//to remove warning  
  Float_t csi_length = fNpy*fDpy + fDeadZone;
  Float_t csi_width = fNpx*fDpx + 2*fDeadZone;
  
  fPadPlane_Width = (csi_width - 2*fDeadZone)/3;
  fPadPlane_Length = (csi_length - fDeadZone)/2;
}


Float_t AliRICHSegmentationV0::GetAnod(Float_t xhit) const
{

// Get anod wire closer to hit

    Float_t wire= (xhit>0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
    return fWireD*wire;
}

void AliRICHSegmentationV0::SetPadSize(Float_t p1, Float_t p2)
{

// Set the pad size

    fDpx=p1;
    fDpy=p2;
}
void AliRICHSegmentationV0::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
  ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx);
  iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy);
  if (iy >  fNpy) iy= fNpy;
  if (iy < -fNpy) iy=-fNpy;
  if (ix >  fNpx) ix= fNpx;
  if (ix < -fNpx) ix=-fNpx;
}
void AliRICHSegmentationV0::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
  x = (ix>0) ? Float_t(ix*fDpx)-fDpx/2. : Float_t(ix*fDpx)-fDpx/2.;
  y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)-fDpy/2.;
}

void AliRICHSegmentationV0::SetHit(Float_t xhit, Float_t yhit)
{// Find the wire position (center of charge distribution)
    fXhit=xhit;
    fYhit=yhit;
}

void AliRICHSegmentationV0::SetPad(Int_t ix, Int_t iy)
{// Move to pad ix, iy

    GetPadC(ix,iy,fX,fY);
}



void AliRICHSegmentationV0::FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{

    // Find the wire position (center of charge distribution)
    Float_t x0a=GetAnod(xhit);
    fXhit=x0a;
    fYhit=yhit;
    // and take fNsigma*sigma around this center
    Float_t x01=x0a  - dx;
    Float_t x02=x0a  + dx;
    Float_t y01=yhit - dy;
    Float_t y02=yhit + dy;
    // find the pads over which the charge distributes
    GetPadI(x01,y01,fIxmin,fIymin);
    GetPadI(x02,y02,fIxmax,fIymax);    
    // Set current pad to lower left corner
    fIx=fIxmin;
    fIy=fIymin;
    GetPadC(fIx,fIy,fX,fY);    
}

void AliRICHSegmentationV0::NextPad()
{
    // Step to next pad in integration region
    if (fIx <= fIxmax) {
	fIx++;
    } else if (fIy <= fIymax) {
	fIx=fIxmin;
	fIy++;
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadC(fIx,fIy,fX,fY);
}

Int_t AliRICHSegmentationV0::MorePads()
{// Are there more pads in the integration region   
  if (fIx >= fIxmax && fIy >= fIymax) {
    return 0;
  } else {
    return 1;
  }
}
//__________________________________________________________________________________________________
void AliRICHSegmentationV0::SigGenInit(Float_t x,Float_t y,Float_t)
{//  Initialises pad and wire position during stepping
    fXt =x;
    fYt =y;
    GetPadI(x,y,fIxt,fIyt);
    fIwt= (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}

Int_t AliRICHSegmentationV0::SigGenCond(Float_t x,Float_t y,Float_t)
{//  Signal will be generated if particle crosses pad boundary or
//  boundary between two wires. 
    Int_t ixt, iyt;
    GetPadI(x,y,ixt,iyt);
    Int_t iwt=(x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1;
    
    if ((ixt != fIxt) || (iyt !=fIyt) || (iwt != fIwt)) {
      return 1;
    } else {
      return 0;
    }
}
void AliRICHSegmentationV0::IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{// Calculate the integration limits
  x1=fXhit-fX-fDpx/2.;
  x2=x1+fDpx;
  y1=fYhit-fY-fDpy/2.;
  y2=y1+fDpy;
}

void AliRICHSegmentationV0::Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])
{
//Is used for the cluster finder, include diagonal elements
    
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
}

Float_t AliRICHSegmentationV0::Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y,Int_t *dummy)
{// Returns the square of the distance between 1 pad labelled by its Channel numbers and a coordinate
  dummy=0;
  Float_t x,y;
  GetPadC(iX,iY,x,y);
  return (x-X)*(x-X) + (y-Y)*(y-Y);
}

void  AliRICHSegmentationV0::GiveTestPoints(Int_t &n, Float_t *x, Float_t *y) const
{

// Test

    n=1;
    x[0]=0.;
    y[0]=x[0];
}

void AliRICHSegmentationV0::Draw(const char* opt) const
{
  Info("Draw","%s",opt);
    
}
