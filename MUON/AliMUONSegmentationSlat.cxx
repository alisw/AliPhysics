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

#include "AliMUONSegmentationSlat.h"
#include "AliMUONSegmentationSlatModule.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "AliRun.h"
#include <TMath.h>
#include <iostream.h>

//___________________________________________
ClassImp(AliMUONSegmentationSlat)

AliMUONSegmentationSlat::AliMUONSegmentationSlat() 
{
// Default constructor
    fSlats=0;            
    fNDiv = new TArrayI(4);   
}

void AliMUONSegmentationSlat::SetPadSize(Float_t p1, Float_t p2)
{
//  Sets the pad (strip) size 
//  
    fDpx=p1;
    fDpy=p2;
}

Float_t AliMUONSegmentationSlat::GetAnod(Float_t xhit) const
{
// Returns for a hit position xhit the position of the nearest anode wire    
    Float_t wire= (xhit>0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
    return fWireD*wire;
}

Float_t AliMUONSegmentationSlat::Dpx(Int_t isec) const
{
//
// Returns x-pad size for given sector isec
// isec = 100*islat+iregion
//
    Int_t islat, iregion;
    islat    = isec/100;
    iregion  = isec%100;
    return Slat(islat)->Dpx(iregion);
}

Float_t AliMUONSegmentationSlat::Dpy(Int_t isec) const
{
//
// Returns y-pad (strip)  size for given sector isec
   return fDpy;
}

void AliMUONSegmentationSlat::SetPadDivision(Int_t ndiv[4])
{
//
// Defines the pad size perp. to the anode wire (y) for different sectors. 
// Pad sizes are defined as integral fractions ndiv of a basis pad size
// fDpx
// 
    for (Int_t i=0; i<4; i++) {
	(*fNDiv)[i]=ndiv[i];
    }
}

void AliMUONSegmentationSlat::GlobalToLocal(
    Float_t x, Float_t y, Float_t z, Int_t &islat, Float_t &xlocal, Float_t &ylocal)
{
//
// Perform local to global transformation for space coordinates
//
    Float_t zlocal;
    Int_t i;
    Int_t index=-1;
// Transform According to slat plane z-position: negative side is shifted down 
//                                                 positive side is shifted up
// by half the overlap
    zlocal = z-fChamber->Z();
    Float_t ys = y-TMath::Sign(fShift,zlocal);

//  Set the signs for the symmetry transformation and transform to first quadrant
    SetSymmetry(x,ys);
    Float_t yabs=TMath::Abs(ys);
    Float_t xabs=TMath::Abs(x);

    Int_t ifirst = (zlocal*ys < Float_t(0))? 0:1;
//
// Find slat number                      
    for (i=ifirst; i<fNSlats; i+=2) {
	index=i;
	if ((yabs >= fYPosition[i]) && (yabs < fYPosition[i]+fSlatY)) break;
    }
    
//
// Transform to local coordinate system

    
    ylocal = yabs-fYPosition[index];
    xlocal = xabs-fXPosition[index];
    islat  = index;
    if (i >= fNSlats) {islat = -1; x=-1; y = -1;}
}

void AliMUONSegmentationSlat::GlobalToLocal(
    Int_t ix, Int_t iy, Int_t &islat, Int_t &ixlocal, Int_t &iylocal)
{
//
// Perform global to local transformation for pad coordinates
//
    Int_t iytemp = TMath::Abs(iy);
    Int_t index  = 0;
    
    iylocal = iytemp;

//
// Find slat number (index) and iylocal  
    for (Int_t i=0; i<fNSlats; i++) {
	iytemp-=Slat(i)->Npy();
	
	
	if (iytemp <= 0) break;
	iylocal = iytemp;
	index=i+1;
    }

    ixlocal=TMath::Abs(ix);
    islat=index;
    
// Done !
}

void AliMUONSegmentationSlat::
LocalToGlobal(Int_t islat, Float_t  xlocal, Float_t  ylocal, Float_t  &x, Float_t  &y, Float_t &z)
{
// Transform from local to global space coordinates
//
// upper plane (y>0) even slat number is shifted down
// upper plane (y>0)  odd slat number is shifted up 
// lower plane (y<0) even slat number is shifted up
// lower plane (y<0)  odd slat number is shifted down
//

    x = (xlocal+fXPosition[islat])*fSym[0];
    if ((TMath::Even(islat) && fSym[1]>0) || (TMath::Odd(islat)&&fSym[1]<0)) {
	y=(ylocal+fYPosition[islat])*fSym[1]-fShift;
	z=-fDz;
    } else {
	y=(ylocal+fYPosition[islat])*fSym[1]+fShift;
	z=fDz;
    }

    z+=fChamber->Z();

}


void AliMUONSegmentationSlat::LocalToGlobal(
    Int_t islat, Int_t ixlocal, Int_t iylocal, Int_t &ix, Int_t &iy)
{
// Transform from local to global pad coordinates
//
    Int_t i;
    iy=iylocal;
    
//
// Find slat number (index) and iylocal  
    for (i=0; i<islat; i++) iy+=Slat(islat)->Npy();

    ix=ixlocal*fSym[0];
    iy=iy*fSym[1];
}


void AliMUONSegmentationSlat::SetSymmetry(Int_t   ix,   Int_t iy)
{
// Set set signs for symmetry transformation
    fSym[0]=TMath::Sign(1,ix);
    fSym[1]=TMath::Sign(1,iy);
    
}

void AliMUONSegmentationSlat::SetSymmetry(Float_t  x, Float_t  y)
{
// Set set signs for symmetry transformation
    fSym[0]=Int_t (TMath::Sign(1.,x));
    fSym[1]=Int_t (TMath::Sign(1.,y));
}

void AliMUONSegmentationSlat::
GetPadI(Float_t x, Float_t y, Float_t z, Int_t &ix, Int_t &iy)
{
// Returns pad coordinates for given set of space coordinates

    Int_t islat, i;
    Float_t xlocal, ylocal;
    
    GlobalToLocal(x,y,z,islat,xlocal,ylocal);
    if (islat == -1) {
	ix=0; iy=0; return;
    }
    
    Slat(islat)->GetPadI(xlocal, ylocal, ix, iy);
    for (i=0; i<islat; i++) iy+=Slat(islat)->Npy();

    ix=ix*Int_t(TMath::Sign(1.,x));    
// Transform y 
    iy=iy*Int_t(TMath::Sign(1.,y));   
}

void AliMUONSegmentationSlat::
GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z)
{
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    Int_t islat, ixlocal, iylocal;
//
// Delegation of transforamtion to slat
    GlobalToLocal(ix,iy,islat,ixlocal,iylocal);
    Slat(islat)->GetPadC(ixlocal, iylocal, x, y);
// Slat offset
    x+=fXPosition[islat];
    y+=fYPosition[islat];    

// Symmetry transformation of quadrants    
    x=x*TMath::Sign(1,ix);
    y=y*TMath::Sign(1,iy);    

// Shift of slat planes
    if ((TMath::Even(islat)&&iy>0) || (TMath::Odd(islat)&&iy<0)) {
	y-=fShift;
	z=-fDz+fChamber->Z();
    } else {
	y+=fShift;
	z=fDz+fChamber->Z();
    }
}

Int_t AliMUONSegmentationSlat::ISector()
{
// Returns current sector during tracking
    Int_t iregion;
    
    iregion =  fCurrentSlat->ISector();
    return 100*fSlatIndex+iregion;
}

Int_t AliMUONSegmentationSlat::Sector(Int_t ix, Int_t iy)
{
    Int_t ixlocal, iylocal, iregion, islat;

    GlobalToLocal(ix,iy,islat,ixlocal,iylocal);
    
    iregion =  Slat(islat)->Sector(ixlocal, iylocal);
    return 100*islat+iregion;
}


void AliMUONSegmentationSlat::SetPad(Int_t ix, Int_t iy)
{
    //
    // Sets virtual pad coordinates, needed for evaluating pad response 
    // outside the tracking program
    Int_t islat, ixlocal, iylocal;

    SetSymmetry(ix,iy);
    
    GlobalToLocal(ix,iy,islat,ixlocal,iylocal);
    fSlatIndex=islat;
    fCurrentSlat=Slat(islat);
    fCurrentSlat->SetPad(ixlocal, iylocal);
}

void  AliMUONSegmentationSlat::SetHit(Float_t xhit, Float_t yhit, Float_t zhit)
{   //
    // Sets current hit coordinates

    Float_t xlocal, ylocal;
    Int_t islat;

    

    GlobalToLocal(xhit,yhit,zhit,islat,xlocal,ylocal);
    fSlatIndex=islat;
    if (islat < 0) printf("\n SetHit: %d", islat);
    
    fCurrentSlat=Slat(islat);
    fCurrentSlat->SetHit(xlocal, ylocal);
}


void AliMUONSegmentationSlat::
FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, Float_t dx, Float_t dy)
{
// Initialises iteration over pads for charge distribution algorithm
//



    Int_t islat;
    Float_t xlocal, ylocal;
    GlobalToLocal(xhit, yhit, zhit, islat, xlocal, ylocal);
    fSlatIndex=islat;
    fCurrentSlat=Slat(islat);
    fCurrentSlat->FirstPad(xlocal, ylocal, dx, dy);

}


void AliMUONSegmentationSlat::NextPad()
{
// Stepper for the iteration over pads
//
    fCurrentSlat->NextPad();
}


Int_t AliMUONSegmentationSlat::MorePads()
// Stopping condition for the iterator over pads
//
// Are there more pads in the integration region
{ 
    return fCurrentSlat->MorePads();
}

void AliMUONSegmentationSlat::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
//  Returns integration limits for current pad
//
    
    fCurrentSlat->IntegrationLimits(x1, x2, y1, y2);

}

void AliMUONSegmentationSlat::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])
{
// Returns list of neighbours of pad with coordinates iX, iY

    Int_t i, xListLocal[10], yListLocal[10], iXlocal, iYlocal, islat;
    
    SetSymmetry(iX,iY);

    GlobalToLocal(iX, iY, islat, iXlocal, iYlocal);
 
    Slat(islat)->Neighbours(iXlocal, iYlocal, Nlist, xListLocal, yListLocal);
    
    for (i=0; i<*Nlist; i++) LocalToGlobal(islat, xListLocal[i], yListLocal[i], Xlist[i], Ylist[i]);

}


Int_t  AliMUONSegmentationSlat::Ix()
{
// Return current pad coordinate ix during stepping
    Int_t ixl,iyl,ix,iy;
    ixl=fCurrentSlat->Ix();
    iyl=fCurrentSlat->Iy();
    
    LocalToGlobal(fSlatIndex, ixl, iyl, ix, iy);
    Int_t ixc, iyc, isc;
    Float_t xc, yc;
    GlobalToLocal(ix, iy, isc, ixc, iyc);
    Slat(isc)->GetPadC(ixc,iyc,xc,yc);
    return ix;
}


Int_t  AliMUONSegmentationSlat::Iy()
{
// Return current pad coordinate iy during stepping
    Int_t ixl,iyl,ix,iy;
    ixl=fCurrentSlat->Ix();
    iyl=fCurrentSlat->Iy();
    LocalToGlobal(fSlatIndex, ixl, iyl, ix, iy);
    return iy;
}



   // Signal Generation Condition during Stepping
Int_t AliMUONSegmentationSlat::SigGenCond(Float_t x, Float_t y, Float_t z)
{ 
//
//  True if signal generation condition fullfilled
    Float_t xlocal, ylocal;
    Int_t islat;
    GlobalToLocal(x, y, z, islat, xlocal, ylocal);
    return Slat(islat)->SigGenCond(xlocal, ylocal, z);
}

// Initialise signal generation at coord (x,y,z)
void  AliMUONSegmentationSlat::SigGenInit(Float_t x, Float_t y, Float_t z)
{
// Initialize the signal generation condition
//
    Float_t xlocal, ylocal;
    Int_t islat;

    GlobalToLocal(x, y, z, islat, xlocal, ylocal);
    Slat(islat)->SigGenInit(xlocal, ylocal, z);
}



void AliMUONSegmentationSlat::Init(Int_t chamber)
{
//    
// Initialize slat modules of quadrant +/+    
// The other three quadrants are handled through symmetry transformations
//
    printf("\n Initialise segmentation Slat \n");
//

//    Initialize Slat modules
    Int_t islat, i;
    Int_t ndiv[4];
// Pad division
    for (i=0; i<4; i++) ndiv[i]=(*fNDiv)[i];
// Half distance between slat planes
    fDz=1.76;
// Slat height    
    fSlatY=40.;
    for (i=0; i<10; i++) fSlatX[i]=0.;
    
    
// Initialize array of slats 
    fSlats  = new TObjArray(fNSlats);
// Maximum number of strips (pads) in x and y
    fNpy=0;   
    fNpx=0;
// for each slat in the quadrant (+,+)    
    for (islat=0; islat<fNSlats; islat++) {
	(*fSlats)[islat] = CreateSlatModule();

	AliMUONSegmentationSlatModule *slat =  Slat(islat);
	// Configure Slat
	slat->SetId(islat);
	
// Foward pad size
	slat->SetPadSize(fDpx, fDpy);
// Forward wire pitch
	slat->SetDAnod(fWireD);
// Foward segmentation 
	slat->SetPadDivision(ndiv);
	slat->SetPcbBoards(fPcb[islat]);
// Initialize slat module
	slat->Init(chamber);
// y-position of slat module relative to the first (closest to the beam)
	fYPosition[islat]=islat*(fSlatY-2.*fShift);
	if (TMath::Odd(islat)) fYPosition[islat] -= 2*fShift;
//
	fNpy+=slat->Npy();
	if (slat->Npx() > fNpx) fNpx=slat->Npx();
	Int_t isec;
	for (isec=0; isec< 4; isec++)
	{
	    fSlatX[islat]+=40.*fPcb[islat][isec];
	}
	
    }
// Set parent chamber number
    AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
    fChamber=&(pMUON->Chamber(chamber));
}





void AliMUONSegmentationSlat::SetNPCBperSector(Int_t *npcb)
{ 
    //  PCB distribution for station 4 (6 rows with 1+3 segmentation regions)
    for (Int_t islat=0; islat<fNSlats; islat++){ 
	fPcb[islat][0] = *(npcb + 4 * islat);
	fPcb[islat][1] = *(npcb + 4 * islat + 1);
	fPcb[islat][2] = *(npcb + 4 * islat + 2);
	fPcb[islat][3] = *(npcb + 4 * islat + 3);
    }
}


void  AliMUONSegmentationSlat::SetSlatXPositions(Float_t *xpos)
{
// Set x-positions of Slats
    for (Int_t islat=0; islat<fNSlats; islat++) fXPosition[islat]=xpos[islat];
}

AliMUONSegmentationSlatModule*  AliMUONSegmentationSlat::Slat(Int_t index) const
{ return ((AliMUONSegmentationSlatModule*) (*fSlats)[index]);}


AliMUONSegmentationSlatModule* AliMUONSegmentationSlat::
CreateSlatModule()
{
    // Factory method for slat module
    return new AliMUONSegmentationSlatModule();
}





