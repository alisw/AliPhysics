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

#include <TMath.h>

#include "AliMUONSegmentationTriggerY.h"
#include "AliMUONTriggerConstants.h"
#include "AliMUON.h"
#include "AliRun.h"

ClassImp(AliMUONSegmentationTriggerY)

//------------------------------------------------------------------
AliMUONSegmentationTriggerY::AliMUONSegmentationTriggerY()
  : AliMUONSegmentationTrigger()
{
// Constructor
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerY::Init(Int_t chamber)
{
// intialize Y segmentation 
  AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
  if(pMUON->GetDebug()>1) printf("%s: Initialize Trigger Chamber Geometry Y\n",ClassName());
  AliMUONSegmentationTrigger::Init(chamber);  

// calculate x & y position of Y strips
  Int_t nModule=AliMUONTriggerConstants::Nmodule();  
  for (Int_t imodule=0; imodule<nModule; imodule++) { 
      Int_t moduleId=AliMUONTriggerConstants::ModuleId(imodule);      
      Int_t nStrip=AliMUONTriggerConstants::NstripY(imodule);
      for (Int_t istrip=0; istrip<nStrip; istrip++){
	  Float_t width=StripSizeY(moduleId,istrip);

	  if (imodule<63) {
	      if (moduleId-10*Int_t(moduleId/10.)==7&&istrip>7) {
		  fXofysmin[imodule][istrip]=
		      ( AliMUONTriggerConstants::XcMin(imodule)+
			(width*2.)*8 + width*(istrip-8) )*fZscale;
		  fXofysmax[imodule][istrip]=
		      ( AliMUONTriggerConstants::XcMin(imodule)+
			(width*2.)*8 + width*(istrip-7) )*fZscale;
	      } else {	      
		  fXofysmin[imodule][istrip]=
		      (AliMUONTriggerConstants::XcMin(imodule)
		       +width*(istrip))*fZscale;
		  fXofysmax[imodule][istrip]=
		      (AliMUONTriggerConstants::XcMin(imodule)
		       +width*(istrip+1))*fZscale;
		  
	      }
	      
	  } else {	
	      fXofysmin[imodule][istrip]=-1.*fXofysmax[imodule-63][istrip];
	      fXofysmax[imodule][istrip]=-1.*fXofysmin[imodule-63][istrip];
	  }      
	  fYofysmin[imodule][istrip] = fYcmin[imodule]*fZscale;
	  fYofysmax[imodule][istrip] = fYcmax[imodule]*fZscale;
      }
  }  
/*
	  if (TMath::Abs(AliMUONTriggerConstants::ModuleId(imodule))==11) {
	      printf("module Id istrip fXofxsmin fXofxsmax fYofxsmin fYofxsmax %d %d %f %f %f %f \n",
		     AliMUONTriggerConstants::ModuleId(imodule),
		     istrip,
		     fXofysmin[imodule][istrip],
		     fXofysmax[imodule][istrip],
		     fYofysmin[imodule][istrip],
		     fYofysmax[imodule][istrip]);
	  }
*/

}


//------------------------------------------------------------------
void AliMUONSegmentationTriggerY::GetPadI(Float_t x,Float_t y,Int_t &ix,Int_t &iy) 
{
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
//  x,y = real coordinates; ix = module number , iy = strip number

  ix = 0;    
  iy = 0;
  Int_t nModule=AliMUONTriggerConstants::Nmodule();
  for (Int_t imodule=0; imodule<nModule; imodule++) {
      Int_t nStrip=AliMUONTriggerConstants::NstripY(imodule);      
    for (Int_t istrip=0; istrip<nStrip; istrip++){
      if (x>fXofysmin[imodule][istrip]&&x<fXofysmax[imodule][istrip]&&
	  y>fYofysmin[imodule][istrip]&&y<fYofysmax[imodule][istrip]){
	ix = AliMUONTriggerConstants::ModuleId(imodule);
	iy = istrip;
      }
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerY::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//  ix = module number , iy = strip number;  x,y = center of strip
  x = 0.;    
  y = 0.;
  Int_t nModule=AliMUONTriggerConstants::Nmodule();
  for (Int_t imodule=0; imodule<nModule; imodule++) {
    if (AliMUONTriggerConstants::ModuleId(imodule)==ix){
      x=fXofysmin[imodule][iy]+(fXofysmax[imodule][iy]-fXofysmin[imodule][iy])/2.;
      y=fYofysmin[imodule][iy]+(fYofysmax[imodule][iy]-fYofysmin[imodule][iy])/2.;
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerY::SetPadSize(Float_t p1, Float_t p2)
{
//  Sets the padsize 
//  
  fDpx=p1;
  fDpy=p2;
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerY::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]){
// Returns list of 10 next neighbours for given Y strip (ix, iy)  
// neighbour number 9 8 7 6 5 (Y strip (ix, iy)) 0 1 2 3 4 in the list
//                  \_______/                    \_______/
//                    left                         right  
// Note : should not be used to return a list of neighbours larger than 16 !

  Int_t absiX = TMath::Abs(iX); 
  Int_t numModule = ModuleNumber(absiX);   // module number Id.
  Int_t nStrip = AliMUONTriggerConstants::NstripY(numModule); //numb of strips
  Int_t iCandidateLeft, iCandidateRight;
  Int_t iNewCandidateRight=0; 
  Int_t iNewCandidateLeft=0;
// first strip number on the right of the left module  
  if ( (absiX-(Int_t(absiX/10))*10)!=1 && absiX!=52 ) 
    iNewCandidateLeft = 
      AliMUONTriggerConstants::NstripY(ModuleNumber(absiX-1))-1;
  Int_t j;
  
  *Nlist = 10;
  for (Int_t i=0; i<10; i++) Xlist[i]=Ylist[i]=0;

  if (iY < nStrip) {

    for (Int_t i=0; i<5; i++) {
      j = i + 5;
      iCandidateRight = iY + (i + 1);
      iCandidateLeft  = iY - (i + 1);
      if (iCandidateRight < nStrip) { // strip in same module  
	Xlist[i] = absiX;
	Ylist[i] = iCandidateRight;  
      } else if ((absiX+1)%10!=8) {   // need to scan the module on the right
	Xlist[i] = absiX+1;
	Ylist[i] = iNewCandidateRight;  
	iNewCandidateRight++;
      }
      
      if (iCandidateLeft >=0 ) { // strip in same module
	Xlist[j] = absiX;
	Ylist[j] = iCandidateLeft;  
      } else if ( iNewCandidateLeft !=0) {
	Xlist[j] = absiX-1;
	Ylist[j] = iNewCandidateLeft;  
	iNewCandidateLeft--;
      }
    }
    
    if (iX<0) {                                  // left side of chamber 
      for (Int_t i=0; i<10; i++) { 
	if (Xlist[i]!=0) Xlist[i]=-Xlist[i]; 
      }
    }
    
  } // iY < nStrip    
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerY::SetPad(Int_t ix, Int_t iy)
{
  // Sets virtual pad coordinates, needed for evaluating pad response 
  // outside the tracking program 
  GetPadC(ix,iy,fX,fY);
  GetPadI(fX,fY,fIx,fIy);
  fSector=Sector(ix,iy);    
}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerY::ISector() 
{ return fSector;}

//------------------------------------------------------------------   

Int_t AliMUONSegmentationTriggerY::Ix()
{ return fIx;}

//------------------------------------------------------------------   

Int_t AliMUONSegmentationTriggerY::Iy()
{ return fIy;}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTriggerY::Dpx(Int_t isec) const
{ 
// returns x size of y strips for sector isec
  if (isec==1) {
    return AliMUONTriggerConstants::StripWidth(1)*fZscale;
  } else if (isec==2) {
    return AliMUONTriggerConstants::StripWidth(1)*fZscale;
  } else if (isec==3) {
    return AliMUONTriggerConstants::StripWidth(1)*fZscale;
  } else if (isec==4) {
    return AliMUONTriggerConstants::StripWidth(2)*fZscale;
  } else {
    return 0.;	
  }       
}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTriggerY::Dpy(Int_t isec) const
{ 
// returns y size of y strips for sector isec
  if (isec==1) {
    return AliMUONTriggerConstants::StripLength(3)*fZscale;
  } else if (isec==2) {
    return AliMUONTriggerConstants::StripLength(2)*fZscale;
  } else if (isec==3) {
    return AliMUONTriggerConstants::StripLength(3)*fZscale;
  } else if (isec==4) {
    return AliMUONTriggerConstants::StripLength(3)*fZscale;
  } else if (isec==5) {
    return AliMUONTriggerConstants::StripLength(3)*fZscale;
  } else {
    return 0.;
  }
}
//------------------------------------------------------------------   
void AliMUONSegmentationTriggerY::GetPadI(Float_t x, Float_t y, Float_t /*z*/, Int_t &ix, Int_t &iy) 
{
  GetPadI(x, y, ix, iy);
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerY::SetHit(Float_t xhit, Float_t yhit)
{ 
// set hits during diintegration
  AliMUONSegmentationTrigger::SetHit(xhit,yhit);
}
//------------------------------------------------------------------   
void AliMUONSegmentationTriggerY::SetHit(Float_t xhit, Float_t yhit, Float_t /*zhit*/)
{
  SetHit(xhit, yhit);
}
//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerY::Sector(Int_t ix, Int_t iy)
{
// Returns sector number for given module
// 
  Int_t absix=TMath::Abs(ix);
  Int_t iwidth=Int_t(StripSizeY(absix,iy));

  if (absix==52) {
    return 1;
  } else if (absix==41||absix==61) {
    return 2;
  } else if (iwidth==2) {
    return 3;
  } else if (iwidth==4) {
    return 4;
  } else {
    return 0;
  }      

}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerY::
IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& x3, Float_t& x4) 
{ 
// returns quantities needed to evaluate neighbour strip response
  Int_t ix,iy;
  Float_t xstrip,ystrip;
  GetPadI(fXhit,fYhit,ix,iy);  
  GetPadC(ix,iy,xstrip,ystrip);  
  x1=fXhit;        // hit x position
  x2=xstrip;       // x coordinate of the main strip
  x3=fX;           // current strip real x coordinate  
  //  width=StripSizeY(ix);   // width of the main strip 

  // find the position of the 2 borders of the current strip
  Float_t xmin = fXofysmin[ModuleNumber(fIx)][fIy];
  Float_t xmax = fXofysmax[ModuleNumber(fIx)][fIy];

  // dist. between the hit and the closest border of the current strip
  x4 = (TMath::Abs(xmax-x1) > TMath::Abs(xmin-x1)) ? 
    TMath::Abs(xmin-x1):TMath::Abs(xmax-x1);    

}




