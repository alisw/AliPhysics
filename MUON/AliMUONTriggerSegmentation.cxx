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

//*********************************************************
// Segmentation classe for trigger chambers.
// The present version is a bit messy because, although the
// x,y positions are local, most of the methods use global 
// strip numbering. The conversion is made via GetPadLoc2Glo.
// To be improved in the future.
//*********************************************************

#include <TArrayI.h>
#include <TArrayF.h>
#include "AliMUONTriggerSegmentation.h"
#include "AliLog.h"

//___________________________________________
ClassImp(AliMUONTriggerSegmentation)

AliMUONTriggerSegmentation::AliMUONTriggerSegmentation() 
  : AliMUONVGeometryDESegmentation(),
    fBending(0),
    fId(0),
    fNsec(7),
    fNDiv(0),
    fDpxD(0),
    fDpyD(0), 
    fDpx(0),
    fDpy(0),
    fNpx(999999),
    fNpy(999999),
    fWireD(0.25),
    fXhit(0.),
    fYhit(0.),
    fIx(0),
    fIy(0),
    fX(0.),
    fY(0.),
    fIxmin(0),
    fIxmax(0),
    fIymin(0),
    fIymax(0),
// add to St345SlatSegmentation
    fLineNumber(0),
    fRpcHalfXsize(0),
    fRpcHalfYsize(0)
{
// add to St345SlatSegmentation
  for (Int_t i=0; i<7; i++) {
      fNstrip[i]=0;
      fStripYsize[i]=0.;   
      fStripXsize[i]=0.;  
      fModuleXmin[i]=0.;
      fModuleXmax[i]=0.;  
      fModuleYmin[i]=0.;  
  }
}


//___________________________________________
AliMUONTriggerSegmentation::AliMUONTriggerSegmentation(Bool_t bending) 
  : AliMUONVGeometryDESegmentation(),
    fBending(bending),
    fId(0),
    fDpx(0),
    fDpy(0),
    fNpx(999999),
    fNpy(999999),
    fWireD(0.25),
    fXhit(0.),
    fYhit(0.),
    fIx(0),
    fIy(0),
    fX(0.),
    fY(0.),
    fIxmin(0),
    fIxmax(0),
    fIymin(0),
    fIymax(0),
// add to St345SlatSegmentation
    fLineNumber(0),
    fRpcHalfXsize(0),
    fRpcHalfYsize(0)
{
  // Non default constructor
  fNsec = 7;  // 4 sector densities at most per slat 
  fNDiv = new TArrayI(fNsec);      
  fDpxD = new TArrayF(fNsec);      
  fDpyD = new TArrayF(fNsec);      
  (*fNDiv)[0]=(*fNDiv)[1]=(*fNDiv)[2]=(*fNDiv)[3]=0;     
  (*fDpxD)[0]=(*fDpxD)[1]=(*fDpxD)[2]=(*fDpxD)[3]=0;       
  (*fDpyD)[0]=(*fDpyD)[1]=(*fDpyD)[2]=(*fDpyD)[3]=0;       
// add to St345SlatSegmentation
  for (Int_t i=0; i<7; i++) {
      fNstrip[i]=0;
      fStripYsize[i]=0.;   
      fStripXsize[i]=0.;  
      fModuleXmin[i]=0.;
      fModuleXmax[i]=0.;  
      fModuleYmin[i]=0.;  
  }
}

//----------------------------------------------------------------------
AliMUONTriggerSegmentation::AliMUONTriggerSegmentation(const AliMUONTriggerSegmentation& rhs) : AliMUONVGeometryDESegmentation(rhs)
{
  AliFatal("Not implemented.");
}
//----------------------------------------------------------------------
AliMUONTriggerSegmentation::~AliMUONTriggerSegmentation() 
{
  // Destructor
  if (fNDiv) delete fNDiv;
  if (fDpxD) delete fDpxD;
  if (fDpyD) delete fDpyD;
}
//----------------------------------------------------------------------
AliMUONTriggerSegmentation& AliMUONTriggerSegmentation::operator=(const AliMUONTriggerSegmentation& rhs)
{
// Protected assignement operator
  if (this == &rhs) return *this;
  AliFatal("Not implemented.");
  return *this;  
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerSegmentation::ModuleColNum(Int_t ixGlo)
{
// returns column number (from 0 to 6) in which the (global) module 
// ixGlo is sitting (could return 7 if ixGlow=isec)
    return TMath::Abs(ixGlo)-Int_t(TMath::Abs(ixGlo)/10)*10-1;
}
//----------------------------------------------------------------------
Bool_t AliMUONTriggerSegmentation::HasPad(Int_t ix, Int_t iy)
{
// check if steping outside the limits (iy=1,2... iyGlo=0,1...)
    Bool_t hasPad = true;    
    Int_t ixGlo = 0;
    Int_t iyGlo = 0; 
    GetPadLoc2Glo(ix, iy, ixGlo, iyGlo);
    if (iyGlo>=fNstrip[ModuleColNum(ixGlo)]) hasPad = false;
    return hasPad;    
}
//____________________________________________________________________________
Float_t AliMUONTriggerSegmentation::Dpx(Int_t isec) const
{
// return x-strip width in sector isec
    Float_t size = (isec<8) ? fStripXsize[isec-1] : fStripXsize[isec-2]/2.;
    return size;
}
//____________________________________________________________________________
Float_t AliMUONTriggerSegmentation::Dpy(Int_t  isec) const
{
// return y-strip width in sector isec
    Float_t size = (isec<8) ? fStripYsize[isec-1] : fStripYsize[isec-2];
    return size;
}
//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadLoc2Glo(Int_t ixLoc, Int_t iyLoc, 
					       Int_t &ixGlo, Int_t &iyGlo)
{    
// converts ixLoc & iyLoc into ixGlo & iyGLo (module,strip number)
    ixGlo = 0; // see AliMUONTriggerConstants::fgkModuleI
    iyGlo = 0; // from 0 to (fNtrip-1) in module   
    if (fBending) { 
	ixGlo = 10*fLineNumber + ixLoc;
	iyGlo = iyLoc - 1;
    } else if (!fBending) {	
	Int_t iCountStrip = 0;	
	for (Int_t iModule=0; iModule<fNsec; iModule++) {		
	    for (Int_t iStrip=0; iStrip<fNstrip[iModule]; iStrip++) {
		if ((ixLoc-1)==iCountStrip) {
		    ixGlo = 10*fLineNumber + iModule + 1;
		    iyGlo = iStrip;
		}
		iCountStrip++;
	    }
	}
    }
//    printf(" in GetPadLoc2Glo ixLoc iyLoc ixGlo iyGlo %i %i %i %i \n",ixLoc,iyLoc,ixGlo,iyGlo);
    
}
//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y) 
{
    Int_t ixGlo = 0;
    Int_t iyGlo = 0;    
    GetPadLoc2Glo(ix, iy, ixGlo, iyGlo);
//    ix = ixGlo;
//    iy = iyGlo;
    
// Returns real coordinates (x,y) for given pad coordinates (ix,iy)
    x = 0.;
    y = 0.;
    Int_t iModule = ModuleColNum(ixGlo); // find column number (0-6)
    if (fBending) {
	if (iModule==0) {
	    x =  fStripXsize[iModule]/ 2.;
	} else {	
	x = fModuleXmax[iModule-1] + fStripXsize[iModule]/2.;
	}
	y =  fModuleYmin[iModule] + 
	    iyGlo*fStripYsize[iModule] + fStripYsize[iModule]/2.;
    } else if (!fBending) {
	if (ModuleColNum(ixGlo)==6 && iyGlo>7) {
	    x = fModuleXmin[iModule] + 8*fStripXsize[iModule] +
		(iyGlo-8)*fStripXsize[iModule]/2. + fStripXsize[iModule]/4.;
	} else {	    
	    x = fModuleXmin[iModule] + 
		iyGlo*fStripXsize[iModule] + fStripXsize[iModule]/2.;
	}
	y =  fModuleYmin[iModule] + fStripYsize[iModule] / 2.;
    }    
    x = x - fRpcHalfXsize;
    y = y - fRpcHalfYsize;

//    printf(" in GetPadC iModule ixloc iyloc ix iy x y %i %i %i %i %i %f %f \n",iModule,ixLoc,iyLoc,ix,iy,x,y);
}

//_____________________________________________________________________________
void AliMUONTriggerSegmentation::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy) 
{
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)

    x = x + fRpcHalfXsize;
    y = y + fRpcHalfYsize;
// find module number    
    Int_t modNum=0;    
    for (Int_t iModule=0; iModule<fNsec; iModule++) { // modules
	if ( x > fModuleXmin[iModule] && x < fModuleXmax[iModule] ) {
	    ix = 10*fLineNumber + iModule;
	    modNum = iModule;	    
	}
    }

// find strip number 
    Float_t yMin = 0.;    
    Float_t yMax = fModuleYmin[modNum];
    Float_t xMin = 0.;
    Float_t xMax = fModuleXmin[modNum];
    if (ix!=0) {
	for (Int_t iStrip=0; iStrip<fNstrip[modNum]; iStrip++) {
	    if (fBending) {
		yMin = yMax;
		yMax = yMin + fStripYsize[modNum];
		if (y > yMin && y < yMax) iy = iStrip;
	    } else if (!fBending) {
		xMin = xMax;
		xMax = xMin + fStripXsize[modNum];
		if (x > xMin && x < xMax) iy = iStrip;
	    }
	} // loop on strips
    } // if ix!=0
//    printf("in GetPadI ix iy x y %i %i %f %f \n",ix,iy,x,y);
}
//-------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadI(Float_t x, Float_t y , Float_t /*z*/, Int_t &ix, Int_t &iy)
{
  GetPadI(x, y, ix, iy);
}
//____________________________________________________________________________
void AliMUONTriggerSegmentation::SetPadSize(Float_t p1, Float_t p2)
{
//  Sets the padsize 
    fDpx=p1;
    fDpy=p2;
}
//-------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetLineNumber(Int_t iLineNumber){
    fLineNumber = iLineNumber;    
}
//-------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetPad(Int_t ix, Int_t iy)
{
  //
  // Sets virtual pad coordinates, needed for evaluating pad response 
  // outside the tracking program 
  GetPadC(ix,iy,fX,fY);
  fSector=Sector(ix,iy);
}
//---------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetHit(Float_t x, Float_t y)
{
  // Set current hit 
  //
  fXhit = x;
  fYhit = y;
}
//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetHit(Float_t xhit, Float_t yhit, Float_t /*zhit*/)
{
  SetHit(xhit, yhit);
}

//--------------------------------------------------------------------------
Int_t AliMUONTriggerSegmentation::Sector(Int_t ix, Int_t iy) 
{
// determine segmentation zone from pad coordinates (from 1 to 8)
    Int_t ixGlo = 0;
    Int_t iyGlo = 0;    
    GetPadLoc2Glo(ix, iy, ixGlo, iyGlo);    
    if (!fBending && ModuleColNum(ixGlo)==6 && iyGlo>7) {
	return 8; // sector 8: diff. strip width within same module
    } else {
	return ModuleColNum(ixGlo)+1;
    }
}

//-----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2) 
{
/*
// need to return (only) x4 = dist. betwwen the hit and the closest border of
// the current strip
    Int_t ix,iy;
    Float_t xstrip,ystrip;
    GetPadI(fXhit,fYhit,ix,iy);  
    GetPadC(ix,iy,xstrip,ystrip);  
    x1=fYhit;        // hit y position
    x2=ystrip;       // y coordinate of the main strip
    x3=fY;           // current strip real y coordinate  
    
    // find the position of the 2 borders of the current strip
    Float_t ymin = ystrip - fStripXsize[ModuleColNum(ixGlo)

fYofxsmin[ModuleNumber(fIx)][fIy];
    Float_t ymax = fYofxsmax[ModuleNumber(fIx)][fIy];
    
    // dist. between the hit and the closest border of the current strip
    x4 = (TMath::Abs(ymax-x1) > TMath::Abs(ymin-x1)) ? 
	TMath::Abs(ymin-x1):TMath::Abs(ymax-x1);
*/
  //  Returns integration limits for current pad
  //
  x1=fXhit-fX-Dpx(fSector)/2.;
  x2=x1+Dpx(fSector);
  y1=fYhit-fY-Dpy(fSector)/2.;
  y2=y1+Dpy(fSector);    
  //    printf("\n Integration Limits %f %f %f %f %d %f", x1, x2, y1, y2, fSector, Dpx(fSector));
}


//-----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]) 
{
  // Returns list of next neighbours for given Pad (iX, iY)
  Int_t i=0;
  //  step right
  if (iX+1 <= fNpx) {
    Xlist[i]=iX+1;
    Ylist[i++]=iY;
  }
  //  step left    
  if (iX-1 > 0) {
    Xlist[i]=iX-1;
    Ylist[i++]=iY;
  } 
  Int_t sector = Sector(iX,iY);
  //  step up
  if (iY+1 <= fNpyS[sector]) {
    Xlist[i]=iX;
    Ylist[i++]=iY+1;
  }
  //  step down    
  if (iY-1 > 0) {
    Xlist[i]=iX;
    Ylist[i++]=iY-1;
  }
  *Nlist=i;
}

//--------------------------------------------------------------------------
void AliMUONTriggerSegmentation::Init(Int_t detectionElementId,
				      Int_t nStrip[7],
				      Float_t stripYsize[7],
				      Float_t stripXsize[7],
				      Float_t offset)
{
//    printf(" fBending: %d \n",fBending);
    
    Int_t nStripMax = 0;
    if (fBending) nStripMax = nStrip[0];
    
    for (Int_t i=0; i<7; i++) {
	fNstrip[i]=nStrip[i];
	fStripYsize[i]=stripYsize[i];
	fStripXsize[i]=stripXsize[i];
	fModuleYmin[0]=0.;
    }
// take care of offset in Y in chamber 5, first module
    fModuleYmin[0] = offset;
    
    Float_t tmp = 0.;  
    Int_t npad = 0;  // number of pad in x and y
    for (Int_t iModule=0; iModule<fNsec; iModule++) { // modules	
	fModuleXmin[iModule] = tmp;      
	npad = npad + fNstrip[iModule];
	if (fBending) {
	    fModuleXmax[iModule] = 
		fModuleXmin[iModule] + fStripXsize[iModule];
	} else if (!fBending) {
	    if (iModule<6) {
		fModuleXmax[iModule] = 
		    fModuleXmin[iModule] + 
		    fStripXsize[iModule]*fNstrip[iModule];
	    } else if (iModule==6) { 
		fModuleXmax[iModule] = 
		    fModuleXmin[iModule] + 
		    (fStripXsize[iModule]*fNstrip[iModule]/2) +
		    (fStripXsize[iModule]/2.*fNstrip[iModule]/2);
	    }	  
	}
	tmp = fModuleXmax[iModule];      
	
// calculate nStripMax in x & y
	if (fBending) {
	    if (fNstrip[iModule] > nStripMax) nStripMax = fNstrip[iModule];
	} else if (!fBending) {
	    for (Int_t iStrip=0; iStrip<fNstrip[iModule]; iStrip++) nStripMax++;
	}
    } // loop on modules

// associate nStripMax
    fNpx = (fBending) ? fNsec : nStripMax;
    fNpy = (fBending) ? nStripMax : 1;
    
// calculate half size in x & y (to shift local coordinate ref. system)
  fRpcHalfXsize = 0;
  fRpcHalfYsize = 0;  
  if (fBending) {
      for (Int_t iModule=0; iModule<fNsec; iModule++)  
	  fRpcHalfXsize = fRpcHalfXsize + fStripXsize[iModule];      
      fRpcHalfYsize = fNstrip[1] * fStripYsize[1];
  } else if (!fBending) {
      fRpcHalfXsize = fModuleXmax[6];
      fRpcHalfYsize = fStripYsize[1];
  }
  fRpcHalfXsize = fRpcHalfXsize / 2.;
  fRpcHalfYsize = fRpcHalfYsize / 2.;  

/*
  printf(" fNpx fNpy fRpcHalfXsize fRpcHalfYsize = %i %i %f %f \n",
	 fNpx,fNpy,fRpcHalfXsize,fRpcHalfYsize);

  for (Int_t iModule=0; iModule<fNsec; iModule++) {
      printf(" iModule fModuleXmin fModuleXmax fModuleYmin fStripXsize fStripYsize %i %f %f %f %f %f\n",
	     iModule,fModuleXmin[iModule],fModuleXmax[iModule],
	     fModuleYmin[iModule],
	     fStripXsize[iModule],fStripYsize[iModule]);
	     }

  for (Int_t iModule=0; iModule<fNsec; iModule++) {
      printf(" iModule fNstrip fStripXsize fStripYsize %i %i %f %f \n",
	     iModule,fNstrip[iModule],
	     fStripXsize[iModule],fStripYsize[iModule]);
  }
*/

  

//  printf("npad = %i",npad);  

  fId = detectionElementId;
}







