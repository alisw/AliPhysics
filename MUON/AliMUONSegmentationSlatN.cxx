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
Revision 1.1  2000/10/06 09:00:47  morsch
Segmentation class for chambers built out of slats.

*/

#include "AliMUONSegmentationSlatN.h"
#include "AliMUONSegmentationSlatModuleN.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include <TMath.h>
#include <iostream.h>

//___________________________________________
ClassImp(AliMUONSegmentationSlatN);



AliMUONSegmentationSlatN::AliMUONSegmentationSlatN()
{
// Default constructor
}


Float_t AliMUONSegmentationSlatN::Dpx(Int_t isec) const
{
//
// Returns y-pad size for given sector isec
   return fDpx;
}

Float_t AliMUONSegmentationSlatN::Dpy(Int_t isec) const
{
//
// Returns x-pad size for given sector isec
// isec = 100*islat+iregion
//
    Int_t islat, iregion;
    islat    = isec/100;
    iregion  = isec%100;
    return Slat(islat)->Dpy(iregion);
}



void AliMUONSegmentationSlatN::GlobalToLocal(
    Int_t ix, Int_t iy, Int_t &islat, Int_t &ixlocal, Int_t &iylocal)
{
//
// Perform local to global transformation for pad coordinates
//
    Int_t iytemp = TMath::Abs(iy); 
    Int_t index  = 0;
    
    iylocal = iytemp;
    ix=TMath::Abs(ix);
    
//
// Find slat number (index) and iylocal  
    for (Int_t i=0; i<fNSlats; i++) {
	if (ix <= Slat(i)->Npx()) {
	    Int_t isec=Slat(i)->Sector(ix,1);
	    iytemp-=Slat(i)->Npy()*(*fNDiv)[isec]/(*fNDiv)[1];
	}
	if (iytemp <= 0) break;
	iylocal = iytemp;
	index=i+1;
    }
    ixlocal=ix;
    islat=index;
// Done !
}

void AliMUONSegmentationSlatN::LocalToGlobal(
    Int_t islat, Int_t ixlocal, Int_t iylocal, Int_t &ix, Int_t &iy)
{
// Local to global transformation for pad coordinates
    
    Int_t i;
    iy=iylocal;
    
//
// Find iy global by adding iy offset from slats below
    for (i=0; i<islat; i++) {
	if (ixlocal <= Slat(i)->Npx()) {
	    Int_t isec=Slat(i)->Sector(ixlocal,1);
	    iy+=Slat(i)->Npy()*(*fNDiv)[isec]/(*fNDiv)[1];
	}
    }
//
// Perform symmetry transformation
    ix=ixlocal*fSym[0];
    iy=iy*fSym[1];
}


void AliMUONSegmentationSlatN::
GetPadI(Float_t x, Float_t y, Float_t z, Int_t &ix, Int_t &iy) 
{
// Returns pad coordinates for given set of space coordinates

    Int_t islat, i;
    Float_t xlocal, ylocal;
// Transform to local coordinates    
    AliMUONSegmentationSlat::GlobalToLocal(x,y,z,islat,xlocal,ylocal);
    Slat(islat)->GetPadI(xlocal, ylocal, ix, iy);
// add to local iy offfset from slats below  
    for (i=0; i<islat; i++) {
	if (ix <= Slat(i)->Npx()) {
	    Int_t isec=Slat(i)->Sector(ix,1);
	    iy+=Slat(i)->Npy()*(*fNDiv)[isec]/(*fNDiv)[1];
	}
    }
// Determine sign depending on quadrant
    ix=ix*Int_t(TMath::Sign((Float_t)1.,x));
    iy=iy*Int_t(TMath::Sign((Float_t)1.,y));    

}

AliMUONSegmentationSlatModule* AliMUONSegmentationSlatN::
CreateSlatModule()
{
    // Factory method for slat module
    return new AliMUONSegmentationSlatModuleN();
}





