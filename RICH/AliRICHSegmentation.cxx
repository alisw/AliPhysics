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

/* $Id$ */

#include "AliRICHSegmentation.h"


ClassImp(AliRICHSegmentation)

//________________________________________________________________________________
AliRICHSegmentation::AliRICHSegmentation()
{ 
   fNpx=144;      // Number of pads in x direction
   fNpy=160;      // Number of pads in y direction
   fSector=-1;
   fDeadZone=3.0; // spacer between photocathod planes cm
   fDpx=0.84;     // Pad x size cm
   fDpy=0.8;      // Pad y size cm
   fWireD=0.84/2; // cm set by SetDAnod	 
}//AliRICHSegmentation::ctor()

void AliRICHSegmentation::Init()
{
  Float_t csi_length = fNpy*fDpy + fDeadZone;
  Float_t csi_width = fNpx*fDpx + 2*fDeadZone;

  fPadPlane_Width = (csi_width - 2*fDeadZone)/3;
  fPadPlane_Length = (csi_length - fDeadZone)/2;
}//void AliRICHSegmentation::Init()

// calculate sector from x-y coordinates

Int_t AliRICHSegmentation::Sector(Float_t x, Float_t y)
{
// Calculate to which sector does the hit (x,y) belong
  
  fSector=-1;
  
  //Parametrized definition

  if (y<-fDeadZone/2)
    {
      if (x> fPadPlane_Width/2 +fDeadZone)
	{
	  if ( x<fPadPlane_Width/2 +fDeadZone + fPadPlane_Width)
	    fSector=0;
	}
      if (x< fPadPlane_Width/2)
	{
	  if (x> -( fPadPlane_Width/2))
	    fSector=2;
	}
      if (x< -( fPadPlane_Width/2 +fDeadZone))
	{
	  if (x> -( fPadPlane_Width/2 +fDeadZone +  fPadPlane_Width))
	    fSector=4;
	}
    }
  else if (y>fDeadZone/2)
    {
      if (x> fPadPlane_Width/2 +fDeadZone)
	{
	  if (x< fPadPlane_Width/2 +fDeadZone +  fPadPlane_Width)
	    fSector=1;
	}
      if (x< fPadPlane_Width/2)
	{
	  if (x> -( fPadPlane_Width/2))
	    fSector=3;
	}
      if (x< -( fPadPlane_Width/2 +fDeadZone))
	{
	  if (x> -( fPadPlane_Width/2 +fDeadZone +  fPadPlane_Width))
	    fSector=5;
	}
    }
  
  return fSector;
}//Int_t AliRICHSegmentation::Sector(Float_t x, Float_t y)


void AliRICHSegmentation::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  real coordinates (x,y) -> (ix,iy) pad numbers
//
// Please check origin of pad numbering !!!

  Int_t sector=Sector(x,y);

  //printf("Sector: %d\n",sector);

  if (sector==0)
    {
      //ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx);
      //iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy);
      ix = Int_t ((x-fDeadZone)/fDpx);
      iy = Int_t ((y+fDeadZone/2)/fDpy)-1;
    }
  if (sector==1)
    {
      ix = Int_t ((x-fDeadZone)/fDpx);
      iy = Int_t ((y-fDeadZone/2)/fDpy);
    }
  if (sector==2)
    {
      ix = (x>=0)? ix = Int_t (x/fDpx) : ix = Int_t (x/fDpx)-1;
      iy = Int_t ((y+fDeadZone/2)/fDpy)-1;
    }
  if (sector==3)
    {
      ix = (x>=0)? ix = Int_t (x/fDpx) : ix = Int_t (x/fDpx)-1;
      iy = Int_t ((y-fDeadZone/2)/fDpy);
    }
  if (sector==4)
    {
      ix = Int_t ((x+fDeadZone)/fDpx)-1;
      iy = Int_t ((y+fDeadZone/2)/fDpy)-1;
    }
  if (sector==5)
    {
      ix = Int_t ((x+fDeadZone)/fDpx)-1;
      iy = Int_t ((y-fDeadZone/2)/fDpy);
    }

  
  //ix = Int_t (x/fDpx);
  //iy = Int_t (y/fDpy);

  //ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx);
  //iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy);

  if (sector==-1)
    {
      ix = fIxmax;
      iy = fIymax;
    }

  if (iy >  fNpy) iy= fNpy;
  if (iy < -fNpy) iy=-fNpy;
  if (ix >  fNpx) ix= fNpx;
  if (ix < -fNpx) ix=-fNpx;
}//void AliRICHSegmentation::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)


void AliRICHSegmentation::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//   
//    pad numbers (ix,iy)->(x,y) real coordinates
//

  Int_t sector=-1;


  Float_t padplane_width = fNpx/3;
  
  if (iy<0)
    {
      if (ix < fNpx/2)
	{
	  if (ix >= padplane_width/2)
	    sector=0;
	}
      if (ix< padplane_width/2)
	{
	  if (ix >= -(padplane_width/2))
	    sector=2;
	}
      if (ix >= -(fNpx/2))
	{
	  if (ix < -(padplane_width/2))
	    sector=4;
	}
    }
  if (iy>=0)
    {
      if (ix < fNpx/2)
	{
	  if (ix >= padplane_width/2)
	    sector=1;
	}
      if (ix< padplane_width/2)
	{
	  if (ix >= -(padplane_width/2))
	    sector=3;
	}
      if (ix >= -(fNpx/2))
	{
	  if (ix < -(padplane_width/2))
	    sector=5;
	}
    }

  if (sector==0)
    {
      x = Float_t(ix)*fDpx+fDpx/2+fDeadZone;
      y = Float_t(iy)*fDpy+fDpy/2-fDeadZone/2;
    }
  if (sector==1)
    {
      x = Float_t(ix)*fDpx+fDpx/2+fDeadZone;
      y = Float_t(iy)*fDpy+fDpy/2+fDeadZone/2;
    }
  if (sector==2)
    {
      x = (ix>=0) ? x = Float_t(ix)*fDpx+fDpx/2 : x = Float_t(ix)*fDpx+fDpx/2;
      y = Float_t(iy)*fDpy+fDpy/2-fDeadZone/2;
    }
  if (sector==3)
    {
      x = (ix>=0) ? x = Float_t(ix)*fDpx+fDpx/2 : x = Float_t(ix)*fDpx+fDpx/2;
      y = Float_t(iy)*fDpy+fDpy/2+fDeadZone/2;
    }
  if (sector==4)
    {
      x = Float_t(ix)*fDpx+fDpx/2-fDeadZone;
      y = Float_t(iy)*fDpy+fDpy/2-fDeadZone/2;
    }
  if (sector==5)
    {
      x = Float_t(ix)*fDpx+fDpx/2-fDeadZone;
      y = Float_t(iy)*fDpy+fDpy/2+fDeadZone/2;
    }
   
}//void AliRICHSegmentation::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)


void AliRICHSegmentation::IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{

// Calculates integration limits

  x1=fXhit-fX-fDpx/2.;
  x2=x1+fDpx;
  y1=fYhit-fY-fDpy/2.;
  y2=y1+fDpy;
}//void AliRICHSegmentation::IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)


Int_t AliRICHSegmentation::SigGenCond(Float_t x,Float_t y,Float_t)
{
//
//  Signal will be generated if particle crosses pad boundary or
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


void AliRICHSegmentation::SigGenInit(Float_t x,Float_t y,Float_t)
{
//
//  Initialises pad and wire position during stepping
    fXt =x;
    fYt =y;
    GetPadI(x,y,fIxt,fIyt);
    fIwt= (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}


Float_t AliRICHSegmentation::GetAnod(Float_t xhit) const
{

// Get anod wire closer to hit

    Float_t wire= (xhit>0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
    return fWireD*wire;
}


void AliRICHSegmentation::SetHit(Float_t xhit, Float_t yhit)
{
//
// Find the wire position (center of charge distribution)
//    Float_t x0a=GetAnod(xhit);
    fXhit=xhit;
    fYhit=yhit;
}

void AliRICHSegmentation::
SetPad(Int_t ix, Int_t iy)
{

// Move to pad ix, iy

    GetPadC(ix,iy,fX,fY);
}



void AliRICHSegmentation::
FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{

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
    //
    // find the pads over which the charge distributes
    GetPadI(x01,y01,fIxmin,fIymin);
    GetPadI(x02,y02,fIxmax,fIymax);    
    // 
    // Set current pad to lower left corner
    fIx=fIxmin;
    fIy=fIymin;
    GetPadC(fIx,fIy,fX,fY);
    
    //if (fSector==2)
      //printf("fIx: %d, fIy: %d fX: %f, fY: %f\n",fIx,fIy,fX,fY);
}

void AliRICHSegmentation::NextPad()
{
    //printf("\n Next Pad \n");
    
    // 
    // Step to next pad in integration region
    if (fIx <= fIxmax) {
//	if (fIx==-1) fIx++;
	fIx++;
    } else if (fIy <= fIymax) {
//	if (fIy==-1) fIy++;
	fIx=fIxmin;
	fIy++;
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadC(fIx,fIy,fX,fY);
}

Int_t AliRICHSegmentation::MorePads()
{
//
// Are there more pads in the integration region

    //printf("\n More  Pads ? \n");
  
  
  if (fIx >= fIxmax && fIy >= fIymax) {
    //printf("There are no more pads\n\n\n\n\n");
    return 0;
  } else {
    //printf("There are more pads\n\n");
    return 1;
  }
}



void AliRICHSegmentation::Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])
{
//Is used for the cluster finder, include diagonal elements    
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
}//void AliRICHSegmentation::Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])

Float_t AliRICHSegmentation::Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *dummy)
{
// Returns the square of the distance between 1 pad
// labelled by its Channel numbers and a coordinate

  Float_t x,y;
  GetPadC(iX,iY,x,y);
  return (x-X)*(x-X) + (y-Y)*(y-Y);
}//Float_t AliRICHSegmentation::Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *dummy)
