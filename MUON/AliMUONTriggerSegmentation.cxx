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

/// \class AliMUONTriggerSegmentation
/// \brief Segmentation classe for trigger chambers.
/// In the present version the method use global strip coordinates except
/// HasPad. The conversion is made via GetPadLoc2Glo.
/// To be improved in the future.

#include "AliMUONTriggerSegmentation.h"
#include "AliLog.h"
#include "Riostream.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerSegmentation)
/// \endcond

//___________________________________________
AliMUONTriggerSegmentation::AliMUONTriggerSegmentation() 
  : AliMUONVGeometryDESegmentation(),
    fBending(0),
    fId(0),
    fNsec(7),
    fNpx(999999),
    fNpy(999999),
    fSector(0),
    fXhit(0.),
    fYhit(0.),
    fIx(0),
    fIy(0),
    fX(0.),
    fY(0.),
// add to St345SlatSegmentation
    fLineNumber(0),
    fRpcHalfXsize(0),
    fRpcHalfYsize(0)
{
/// Default constructor

// add to St345SlatSegmentation
  for (Int_t i=0; i<7; i++) {
      fNstrip[i]=0;
      fStripYsize[i]=0.;   
      fStripXsize[i]=0.;  
      fModuleXmin[i]=0.;
      fModuleXmax[i]=0.;  
      fModuleYmin[i]=0.;  
  }

  AliDebug(1, Form("default (empty) ctor this = %p", this));
}


//___________________________________________
AliMUONTriggerSegmentation::AliMUONTriggerSegmentation(Bool_t bending) 
  : AliMUONVGeometryDESegmentation(),
    fBending(bending),
    fId(0),
    fNsec(7),
    fNpx(999999),
    fNpy(999999),
    fSector(0),
    fXhit(0.),
    fYhit(0.),
    fIx(0),
    fIy(0),
    fX(0.),
    fY(0.),
// add to St345SlatSegmentation
    fLineNumber(0),
    fRpcHalfXsize(0),
    fRpcHalfYsize(0)
{
/// Standard constructor

// add to St345SlatSegmentation
  for (Int_t i=0; i<7; i++) {
      fNstrip[i]=0;
      fStripYsize[i]=0.;   
      fStripXsize[i]=0.;  
      fModuleXmin[i]=0.;
      fModuleXmax[i]=0.;  
      fModuleYmin[i]=0.;  
  }

  AliDebug(1, Form("ctor this = %p", this) ); 
}

//----------------------------------------------------------------------
AliMUONTriggerSegmentation::~AliMUONTriggerSegmentation() 
{
/// Destructor

  AliDebug(1, Form("dtor this = %p", this) ); 
}
//----------------------------------------------------------------------
Int_t AliMUONTriggerSegmentation::ModuleColNum(Int_t ix)
{
/// returns column number (from 0 to 6) in which the (global) module 
/// ix is sitting (could return 7 if ix=isec)

    return TMath::Abs(ix)-Int_t(TMath::Abs(ix)/10)*10-1;
}
//----------------------------------------------------------------------
Bool_t AliMUONTriggerSegmentation::HasPad(Int_t ix, Int_t iy)
{
/// check if steping outside the limits (iy=1,2... iy=0,1...)

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
/// return x-strip width in sector isec

    Float_t size = (isec<8) ? fStripXsize[isec-1] : fStripXsize[isec-2]/2.;
    return size;
}
//____________________________________________________________________________
Float_t AliMUONTriggerSegmentation::Dpy(Int_t  isec) const
{
/// return y-strip width in sector isec

    Float_t size = (isec<8) ? fStripYsize[isec-1] : fStripYsize[isec-2];
    return size;
}
//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadLoc2Glo(Int_t ixLoc, Int_t iyLoc, 
					       Int_t &ixGlo, Int_t &iyGlo)
{    
/// converts ixLoc & iyLoc into ixGlo & iyGLo (module,strip number)

    ixGlo = 0; // 
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
//    printf(" in GetPadLoc2Glo fbending ixLoc iyLoc ixGlo iyGlo %i %i %i %i \n",fBending,ixLoc,iyLoc,ixGlo,iyGlo);
}

//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadGlo2Loc(Int_t ixGlo, Int_t iyGlo, 
					       Int_t &ixLoc, Int_t &iyLoc)
{    
/// converts ixGlo & iyGlo into ixLoc & iyLoc 

    ixLoc = 0; 
    iyLoc = 0; 
    if (fBending) { 
	ixLoc = ModuleColNum(ixGlo) + 1;
	iyLoc = iyGlo + 1;
    } else if (!fBending) {	
	Int_t iCountStrip = 1;	
	for (Int_t iModule=0; iModule<fNsec; iModule++) {		
	    for (Int_t iStrip=0; iStrip<fNstrip[iModule]; iStrip++) {
		if ((iModule==ModuleColNum(ixGlo))&&(iStrip==iyGlo)) {
		    iyLoc = 1;		    
		    ixLoc = iCountStrip;
		}		
		iCountStrip++;
	    }
	}
    }    
//    printf(" in GetPadGlo2Loc fBending ixGlo iyGlo ixLoc iyLoc %i %i %i %i %i \n",fBending,ixGlo,iyGlo,ixLoc,iyLoc);
}

//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y) 
{
/// Returns local real coordinates (x,y) for local pad coordinates (ix,iy)

    x = 0.;
    y = 0.;
    Int_t iModule = ModuleColNum(ix); // find column number (0-6)
    if (fBending) {
	if (iModule==0) {
	    x =  fStripXsize[iModule]/ 2.;
	} else {	
	x = fModuleXmax[iModule-1] + fStripXsize[iModule]/2.;
	}
	y =  fModuleYmin[iModule] + 
	    iy*fStripYsize[iModule] + fStripYsize[iModule]/2.;
    } else if (!fBending) {
	if (ModuleColNum(ix)==6 && iy>7) {
	    x = fModuleXmin[iModule] + 8*fStripXsize[iModule] +
		(iy-8)*fStripXsize[iModule]/2. + fStripXsize[iModule]/4.;
	} else {	    
	    x = fModuleXmin[iModule] + 
		iy*fStripXsize[iModule] + fStripXsize[iModule]/2.;
	}
	y =  fModuleYmin[iModule] + fStripYsize[iModule] / 2.;
    }    
    x = x - fRpcHalfXsize;
    y = y - fRpcHalfYsize;

//    printf(" in GetPadC fBending ix iy x y %i %i %i %f %f \n",fBending,ix,iy,x,y);
}

//_____________________________________________________________________________
void AliMUONTriggerSegmentation::GetPadI(Float_t x, Float_t y, Int_t &ix, Int_t &iy) 
{
///  Returns global pad coordinates (ix,iy) for local real coordinates (x,y)

    ix = -1;
    iy = -1;

    x = x + fRpcHalfXsize;
    y = y + fRpcHalfYsize;
// find module number    
    Int_t modNum=0;    
    for (Int_t iModule=0; iModule<fNsec; iModule++) { // modules
	if ( x > fModuleXmin[iModule] && x < fModuleXmax[iModule] ) {
	    ix = 10*fLineNumber + iModule + 1;
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
		if (modNum==6 && iStrip>7) {
		    xMax = xMin + fStripXsize[modNum]/2.;
		} else {		    
		    xMax = xMin + fStripXsize[modNum];
		}		
		if (x > xMin && x < xMax) iy = iStrip;
	    } //
	} // loop on strips
    } // if ix!=0
//    printf("in GetPadI fBending x y ix iy %i %f %f %i %i \n",fBending,x,y,ix,iy);
}
//-------------------------------------------------------------------------
void AliMUONTriggerSegmentation::GetPadI(Float_t x, Float_t y , Float_t /*z*/, Int_t &ix, Int_t &iy)
{
///  Returns global pad coordinates (ix,iy) for local real coordinates (x,y)

  GetPadI(x, y, ix, iy);
}

//-------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetLineNumber(Int_t iLineNumber){
/// Set line number

    fLineNumber = iLineNumber;    
}
//-------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetPad(Int_t ix, Int_t iy)
{
/// Sets virtual pad coordinates, needed for evaluating pad response 
/// outside the tracking program 

  GetPadC(ix,iy,fX,fY);
  fIx = ix; // used in IntegrationLimits
  fIy = iy;    
  fSector=Sector(ix,iy);
}
//---------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetHit(Float_t x, Float_t y)
{
/// Set current hit 

  fXhit = x;
  fYhit = y;
}
//----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::SetHit(Float_t xhit, Float_t yhit, Float_t /*zhit*/)
{
/// Set current hit 

  SetHit(xhit, yhit);
}

//--------------------------------------------------------------------------
Int_t AliMUONTriggerSegmentation::Sector(Int_t ix, Int_t iy) 
{
/// determine segmentation zone from pad coordinates (from 1 to 8)

    if (!fBending && ModuleColNum(ix)==6 && iy>7) {
	return 8; // sector 8: diff. strip width within same module
    } else {
	return ModuleColNum(ix)+1;
    }
}

//-----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::IntegrationLimits(Float_t& x1,Float_t& x2,
                                                   Float_t& x3, Float_t& x4) 
{
/// need to return (only) x4 = dist. betwwen the hit and the closest border of
/// the current strip

    Int_t ix,iy;
    Float_t xstrip,ystrip;
    GetPadI(fXhit,fYhit,ix,iy);
    GetPadC(ix,iy,xstrip,ystrip);
    AliDebug(1,Form("fXhit,Yhit=%e,%e xstrip,ystrip=%e,%e\n",
                    fXhit,fYhit,xstrip,ystrip));
    x1= (fBending) ? fYhit : fXhit;  // hit y (bending) / x (!bending) position
    x2= (fBending) ? ystrip : xstrip; // y or x coord. of the main strip
    x3= (fBending) ? fY : fX;          // current strip real y or x coord.
    Int_t modNum = ModuleColNum(fIx);

    // find the position of the 2 borders of the current strip
    Float_t min = 0.;
    Float_t max = 0.;    
    if (fBending) {
	min = x3 - fStripYsize[modNum]/2.;
	max = x3 + fStripYsize[modNum]/2.;	
    } else {
	if (modNum==6 && fIy>7) { // equivalent to fSector == 8
 	    min = x3 - fStripXsize[modNum]/4.;
	    max = x3 + fStripXsize[modNum]/4.;	
	} else 
	    min = x3 - fStripXsize[modNum]/2.;
	    max = x3 + fStripXsize[modNum]/2.;		
    }    
    // dist. between the hit and the closest border of the current strip
    x4 = (TMath::Abs(max-x1) > TMath::Abs(min-x1)) ? 
      TMath::Abs(min-x1):TMath::Abs(max-x1);
  
    AliDebug(1,Form("Bending %d x1=%e x2=%e x3=%e x4=%e xmin,max=%e,%e\n",
                    fBending,x1,x2,x3,x4,min,max));
}


//-----------------------------------------------------------------------------
void AliMUONTriggerSegmentation::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]) 
{
/// <pre>
///-----------------BENDING-----------------------------------------
/// Returns list of 10 next neighbours for given X strip (ix, iy)  
/// neighbour number 4 in the list -                     
/// neighbour number 3 in the list  |                    
/// neighbour number 2 in the list  |_ Upper part             
/// neighbour number 1 in the list  |            
/// neighbour number 0 in the list -           
///      X strip (ix, iy) 
/// neighbour number 5 in the list -       
/// neighbour number 6 in the list  | _ Lower part
/// neighbour number 7 in the list  |
/// neighbour number 8 in the list  | 
/// neighbour number 9 in the list -
///
///-----------------NON-BENDING-------------------------------------
/// Returns list of 10 next neighbours for given Y strip (ix, iy)  
/// neighbour number 9 8 7 6 5 (Y strip (ix, iy)) 0 1 2 3 4 in the list
///                 \\_______/                    \\_______/
///                    left                         right

    Int_t absiX = TMath::Abs(iX); 
    Int_t modNum = ModuleColNum(absiX); // from 0 to 6
    Int_t nStrip = fNstrip[modNum];    
    
    if (fBending) {
	Int_t iCandidateUp, iCandidateDo;
	Int_t j;
	
	*Nlist = 10;
	for (Int_t i=0; i<10; i++) Xlist[i]=Ylist[i]=0;
	
	if (iY < nStrip) {	    
	    for (Int_t i=0; i<5; i++) {
		j = i + 5;
		iCandidateUp = iY + (i + 1);
		iCandidateDo = iY - (i + 1);
		if (iCandidateUp < nStrip) { 
		    Xlist[i] = iX;
		    Ylist[i] = iCandidateUp;  
		}
		if (iCandidateDo >= 0) { 
		    Xlist[j] = iX;
		    Ylist[j] = iCandidateDo;  
		}
	    }	    
	} // iY < nStrip	

    } else { // non-bending
	 
      Int_t iCandidateLeft, iCandidateRight;
      Int_t iNewCandidateRight=0; 
      Int_t iNewCandidateLeft=0;
// first strip number on the right of the left module  
      if ( modNum!=0 && absiX!=52 ) 
	  iNewCandidateLeft = fNstrip[modNum-1]-1;
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
	      } else if (modNum!=6) {   // need to scan the module on the right
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
    } // non-bending

//    for (Int_t i=0; i<10; i++) {
//	printf("AliMUONTriggerSegmentation LOC fBending i ix iy = %i %i %i %i \n",fBending,i,Xlist[i],Ylist[i]);
//   }
}

//--------------------------------------------------------------------------
void AliMUONTriggerSegmentation::Init(Int_t detectionElementId,
				      Int_t nStrip[7],
				      Float_t stripYsize[7],
				      Float_t stripXsize[7],
				      Float_t offset)
{
/// Initialize

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
//   fNpx = (fBending) ? fNsec : nStripMax;
//   fNpy = (fBending) ? nStripMax : 1;
    fNpx = 124; // tot num of modules (like with old segmentation)
    fNpy = 64; // max number of y strips within one module

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

  fId = detectionElementId;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::Print(Option_t*) const
{
/// Printing

  cout << "fId=" << fId << " fBending=" << fBending << " fNsec=" 
  << fNsec << " Nx,Ny=" << fNpx << "," << fNpy 
  << " LineNumber=" << fLineNumber 
  << " fRpcHalfSize(X,Y)=" << fRpcHalfXsize << "," << fRpcHalfYsize
  << endl;
  
  for (Int_t iModule=0; iModule<fNsec; iModule++) 
  { 
    cout << "Module " << iModule 
    << " xmin,xmax=" << fModuleXmin[iModule] 
    << "," << fModuleXmax[iModule] 
    << " ymin=" << fModuleYmin[iModule]
    << " StripSize(X,Y)=(" << fStripXsize[iModule] << ","
    << fStripYsize[iModule] << ")"
    << endl;
  }                    
}





