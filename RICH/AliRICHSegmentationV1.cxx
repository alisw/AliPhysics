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
  Revision 1.6  2001/01/24 21:00:29  jbarbosa
  Redefinition of sectors and pad coordinates/real coordinates transformations.

  Revision 1.5  2001/01/22 21:37:39  jbarbosa
  Added parametrised definiton sectors

  Revision 1.4  2001/01/22 21:35:39  jbarbosa
  Added deadzone size to data members

  Revision 1.3  2000/10/03 21:44:09  morsch
  Use AliSegmentation and AliHit abstract base classes.

  Revision 1.2  2000/10/02 15:48:55  jbarbosa
  Fixed coding conventions.

  Revision 1.1  2000/06/12 15:34:28  jbarbosa
  Cleaned up version.

*/

#include "AliRICHSegmentationV1.h"


//--------------------------------------------------------
ClassImp(AliRICHSegmentationV1)

//________________________________________________________________________________
AliRICHSegmentationV1::AliRICHSegmentationV1()
{ 

// Default constructor for AliRICHSegmantionV1 (with dead zones)

  fNpx=144;
  fNpy=160;
  //fNpx=80;
  //fNpy=48;
  fDeadZone=2.6;
  fSector=-1;
}

void AliRICHSegmentationV1::Init(Int_t id)
{

// Initialisation of chambers

  //printf("*            Initialising SegmentationV1 (dead zones) in chamber %d              *\n",id+1);

  // parametrised definition

  Float_t csi_length = fNpy*fDpy + fDeadZone;
  Float_t csi_width = fNpx*fDpx + 2*fDeadZone;

  fPadPlane_Width = (csi_width - 2*fDeadZone)/3;
  fPadPlane_Length = (csi_length - fDeadZone)/2;

}

// calculate sector from x-y coordinates

Int_t AliRICHSegmentationV1::Sector(Float_t x, Float_t y)
{

// Calculate in which sector is the hit
  
  fSector=-1;
  
  //old numerical definition

  /*if (x<-fDeadZone/2)
    {
      if (y>22.75)
	{
	  if (y<63.1)
	    fSector=0;
	}
      if (y<20.15)
	{
	  if (y>(-20.15))
	    fSector=2;
	}
      if (y<(-22.75))
	{
	  if (y>(-63.1))
	    fSector=4;
	}
    }
  else if (x>fDeadZone/2)
    {
      if (y>22.75)
	{
	  if (y<63.1)
	    fSector=1;
	}
      if (y<20.15)
	{
	  if (y>(-20.15))
	    fSector=3;
	}
      if (y<(-22.75))
	{
	  if (y>(-63.1))
	    fSector=5;
	}
    }*/

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

  
  //if (fSector==2)
    //printf("x:%f, y:%f, sector:%d\n",x,y,fSector);

  return fSector;
}


void AliRICHSegmentationV1::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
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
}

void AliRICHSegmentationV1::
GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
//

  //Int_t sector=Sector(ix*.8,iy*.84);

  Int_t sector=-1;

 // old numerical definition

 /*if (ix<=0)
    {
      if (iy<=72)
	{
	  if (iy>24)
	    sector=0;
	}
      if (iy<=24)
	{
	  if (iy>-24)
	    sector=2;
	}
      if (iy<=-24)
	{
	  if (iy>-72)
	    sector=4;
	}
    }
  if (ix>0)
    {
      if (iy<=72)
	{
	  if (iy>24)
	    sector=1;
	}
      if (iy<=24)
	{
	  if (iy>-24)
	    sector=3;
	}
      if (iy<=-24)
	{
	  if (iy>-72)
	    sector=5;
	}
    }*/
  
  //  parametrised definition

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
 
  
  //if (sector==2)
    //printf("fSector:%d x:%f y:%f\n",fSector,x,y);
  
}

void AliRICHSegmentationV1::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{

// Calculates integration limits

/*
  x1=fXt-fX-fDpx/2.;
  x2=x1+fDpx;
  y1=fYt-fY-fDpy/2.;
  y2=y1+fDpy;    
*/
  //Int_t sector=Sector(fX,fY);

  //printf("Sector:%d\n",sector);

  x1=fXhit-fX-fDpx/2.;
  x2=x1+fDpx;
  y1=fYhit-fY-fDpy/2.;
  y2=y1+fDpy;
}






