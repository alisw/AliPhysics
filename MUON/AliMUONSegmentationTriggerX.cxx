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
Revision 1.1.2.1  2000/06/09 21:51:03  morsch
Code from AliMUONSegResTrigger.cxx

*/


/*
Old Log:
Revision 1.1.2.4  2000/04/26 12:33:25  morsch
Minor changes in some methods (CP)

Revision 1.1.2.3  2000/03/20 18:14:16  morsch
Missing sector added.

Revision 1.1.2.2  2000/02/20 07:50:49  morsch
Bugs in Dpx, Dpy and ISector methods corrected (P.C.)

Revision 1.1.2.1  2000/02/17 14:33:49  morsch
Draft version from P. Crochet

*/

#include "AliMUONSegmentationTriggerX.h"
#include "TMath.h"
#include "TRandom.h"
#include "TArc.h"
#include "AliMUONChamber.h"
#include <iostream.h> 
ClassImp(AliMUONSegmentationTriggerX)
    
void AliMUONSegmentationTriggerX::Init(AliMUONChamber* Chamber)
{
  cout << "Initialize Trigger Chamber Geometry X " << "\n";    
  AliMUONSegmentationTrigger::Init(Chamber);

// calculate x & y position of X strips
  for (Int_t imodule=0; imodule<fgNmodule; imodule++) {
    Float_t width=StripSizeX(fgNum[imodule]);     
    for (Int_t istrip=0; istrip<fgNstripx[imodule]; istrip++){    
      fXofxsmin[imodule][istrip] = fgXcmin[imodule]*fZscale;
      fXofxsmax[imodule][istrip] = fgXcmax[imodule]*fZscale;
      
      fYofxsmin[imodule][istrip] = (fgYcmin[imodule]+width*(istrip))*fZscale;
      fYofxsmax[imodule][istrip] = (fgYcmin[imodule]+width*(istrip+1))*fZscale;
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::GetPadIxy(Float_t x,Float_t y,Int_t &ix,Int_t &iy){
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
//  x,y = real coordinates; ix = module number , iy = strip number
  ix = 0;    
  iy = 0;
  for (Int_t imodule=0; imodule<fgNmodule; imodule++) {
    for (Int_t istrip=0; istrip<fgNstripx[imodule]; istrip++){
      if (x>fXofxsmin[imodule][istrip]&&x<fXofxsmax[imodule][istrip]&&
      	  y>fYofxsmin[imodule][istrip]&&y<fYofxsmax[imodule][istrip]){
	  ix = fgNum[imodule];
	  iy = istrip;
      }
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y){
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//  ix = module number , iy = strip number;  x,y = center of strip
  x = 0.;    
  y = 0.; 
  
  for (Int_t imodule=0; imodule<fgNmodule; imodule++) {
    if (fgNum[imodule]==ix){
      x=fXofxsmin[imodule][iy]+(fXofxsmax[imodule][iy]-fXofxsmin[imodule][iy])/2.;
      y=fYofxsmin[imodule][iy]+(fYofxsmax[imodule][iy]-fYofxsmin[imodule][iy])/2.;
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::SetPadSize(Float_t p1, Float_t p2)
{
//  Sets the padsize 
//  
  fDpx=p1;
  fDpy=p2;
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[2], Int_t Ylist[2]){
// Returns list of next neighbours for given Pad (ix, iy)  

  Int_t absiX=TMath::Abs(iX); 
  *Nlist = 0;
    
  if (absiX!=0) {                           
    Int_t numModule=ModuleNumber(absiX);
    
    if (iY<fgNstripx[numModule]-1) {          // strip up in same module 
      *Nlist=1;
      Xlist[0]=absiX;
      Ylist[0]=iY+1;
    } 
    
    if (iY>0) {                               // strip down in same module
      *Nlist=*Nlist+1;
      Xlist[*Nlist-1]=absiX;
      Ylist[*Nlist-1]=iY-1;
    } 
    
    if (iX<0) {                               // left side of chamber 
      for (Int_t i=0; i<*Nlist; i++) {Xlist[i]=-Xlist[i];}
    }
  }
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::SetPad(Int_t ix, Int_t iy)
{
  // Sets virtual pad coordinates, needed for evaluating pad response 
  // outside the tracking program
  GetPadCxy(ix,iy,fx,fy);
  GetPadIxy(fx,fy,fix,fiy);
  fSector=Sector(ix,iy);
}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::ISector()
{ return fSector;}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::Ix()
{ return fix;}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::Iy()
{ return fiy;}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTriggerX::Dpx(Int_t isec)
{ // returns x size of x strips for sector isec
    
  if (isec==1) {
    return 17.0*fZscale;
  } else if (isec==2) {
    return 34.0*fZscale;
  } else if (isec==3) {
    return 34.0*fZscale;
  } else if (isec==4) {
    return 34.0*fZscale;
  } else if (isec==5) {
    return 34.0*fZscale;
  } else if (isec==6) {
    return 68.0*fZscale;
  } else {
    return 0.;
  }
}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTriggerX::Dpy(Int_t isec)
{ // returns y size of x strips for sector isec
    // cout << " In AliMUONSegmentationTriggerX::Dpx" << "\n";  

  if (isec==1) {
    return 1.0625*fZscale;
  } else if (isec==2) {
    return 1.0625*fZscale;
  } else if (isec==3) {
    return 1.0625*fZscale;
  } else if (isec==4) {
    return 2.125*fZscale;
  } else if (isec==5) {
    return 4.25*fZscale;
  } else if (isec==6) {
    return 4.25*fZscale;
  } else {
    return 0.;	
  }   
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::SetHit(Float_t xhit, Float_t yhit)
{ AliMUONSegmentationTrigger::SetHit(xhit,yhit);}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::Sector(Int_t ix, Int_t iy)
{
// Returns sector number for given module
// 
  Int_t absix=TMath::Abs(ix);
  Int_t iwidth=Int_t(StripSizeX(absix));
  
  if (absix==52) {
    return 1;
  } else if (absix==41||absix==61) {
    return 2;
  } else if (iwidth==1) {
    return 3;
  } else if (iwidth==2) {
    return 4;
  } else if ((absix>=11&&absix<17)||(absix>=91&&absix<97)) {
    return 5;
  } else if (iwidth==4) {
    return 6;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::
IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& x3, Float_t& width) 
{ // returns quantities needed to evaluate neighbour strip response
//  cout << " In AliMUONSegmentationTriggerX::IntegrationLimits" << "\n";
  Int_t ix,iy;
  Float_t xstrip,ystrip;
  GetPadIxy(fxhit,fyhit,ix,iy);  
  GetPadCxy(ix,iy,xstrip,ystrip);  
  x1=fyhit;        // hit y position
  x2=ystrip;       // y coordinate of the main strip
  x3=fy;           // current strip real y coordinate  
  width=StripSizeX(ix);   // width of the main strip 
}








