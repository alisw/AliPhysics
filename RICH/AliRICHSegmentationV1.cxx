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

#include "AliRICHSegmentationV1.h"

ClassImp(AliRICHSegmentationV1)
//__________________________________________________________________________________________________
AliRICHSegmentationV1::AliRICHSegmentationV1()
{// Default constructor for AliRICHSegmantionV1 (with dead zones)

   fNpx=144;      // number of pads along X direction 
   fNpy=160;      // number of pads along Y direction 
   fDeadZone=3.0; // space between CsI photocathods in cm
   fDpx=0.84;     // pad width in cm
   fDpy=0.80;     // pad heights in cm
   fWireD=0.84/2;	 
   fSector=-1;
   Init(0);       // ??? remove 0
}
//__________________________________________________________________________________________________
void AliRICHSegmentationV1::Init(Int_t /*id*/)
{//Recalculates all the values after some of them have been changed
    
   Float_t csi_length = fNpy*fDpy + fDeadZone;
   Float_t csi_width = fNpx*fDpx + 2*fDeadZone;

   fPadPlane_Width = (csi_width - 2*fDeadZone)/3;
   fPadPlane_Length = (csi_length - fDeadZone)/2;
}
//__________________________________________________________________________________________________
Int_t AliRICHSegmentationV1::Sector(Float_t x, Float_t y)
{// Calculate in which sector is the hit
  
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
}//Sector()
//__________________________________________________________________________________________________
void AliRICHSegmentationV1::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{

  ix=9999; //PH Fake values which should not be returned
  iy=9999;

  Int_t sector=Sector(x,y);


  if (sector==0)
    {
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


  if (sector==-1)
    {
      ix = fIxmax;
      iy = fIymax;
    }

  if (iy >  fNpy) iy= fNpy;
  if (iy < -fNpy) iy=-fNpy;
  if (ix >  fNpx) ix= fNpx;
  if (ix < -fNpx) ix=-fNpx;
}//GetPadI()
//__________________________________________________________________________________________________
void AliRICHSegmentationV1::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
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
      //x = (ix>0) ? Float_t(ix*fDpx)-fDpx/2. : Float_t(ix*fDpx)-fDpx/2.;
      //y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)-fDpy/2.;
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
}//GetPadC
//__________________________________________________________________________________________________
void AliRICHSegmentationV1::IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
  x1=fXhit-fX-fDpx/2.;
  x2=x1+fDpx;
  y1=fYhit-fY-fDpy/2.;
  y2=y1+fDpy;
}
