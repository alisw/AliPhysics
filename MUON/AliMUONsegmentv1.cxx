/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliMUONv0.h"
#include "AliMUONsegmentv1.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"

//___________________________________________
ClassImp(AliMUONsegmentationV1)

AliMUONsegmentationV1::AliMUONsegmentationV1()
  // initizalize the class with default settings
{
fNzone=1;
fDAnod=0.0;
fDpx=0.0;
fDpx=0.0; // forces crash if not initialized by user
fNZoneCut[0]=0;
fSensOffset=0;
}


void AliMUONsegmentationV1::Init(AliMUONchamber* Chamber)
{
  // valid only for T5/6
  frSensMin2 = (Chamber->frMin + 6.)*(Chamber->frMin + 6.);
  frSensMax2 = (Chamber->frMax - 6.)*(Chamber->frMax - 6.);
  fNpx=(Int_t) (Chamber->frMax/fDpx) + 1;
  fNpy=(Int_t) (Chamber->frMax/fDpy) + 1;
  //    fNwire=3;
  DefaultCut();
}

void AliMUONsegmentationV1::DefaultCut(void)
{
SetNzone(3);
AddCut(0,7*6,18*8);
AddCut(0,10*6,15*8);
AddCut(0,12*6,12*8);
AddCut(0,13*6,9*8);
AddCut(0,14*6,6*8);
AddCut(1,8*6,21*12);
AddCut(1,12*6,18*12);
AddCut(1,16*6,15*12);
AddCut(1,18*6,12*12);
AddCut(1,20*6,9*12);
AddCut(1,22*6,6*12);
SetSensOffset(3.0);
SetDAnod(0.325);
}

Int_t AliMUONsegmentationV1::GetiAnod(Float_t xhit)
{
Int_t kwire=Int_t((TMath::Abs(xhit)-fSensOffset)/fDAnod)+1;
return (xhit>0) ? kwire : -kwire ;
}

Float_t AliMUONsegmentationV1::GetAnod(Float_t xhit)
{
  Int_t kwire=Int_t((TMath::Abs(xhit)-fSensOffset)/fDAnod)+1; // to be compatible ...
    return (xhit>0) ? fDAnod*(kwire-0.5)+fSensOffset : -fDAnod*(kwire-0.5)-fSensOffset ;
}

// For chamber T5/6 p1 and p2 should be same for each zone
void AliMUONsegmentationV1::SetPADSIZ(Float_t p1, Float_t p2)
{
  fDpx=p1;
  fDpy=p2;
}

void AliMUONsegmentationV1::
    GetPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
    ix = (x>0)? Int_t((x-fSensOffset)/fDpx)+1 : Int_t((x+fSensOffset)/fDpx)-1;
    iy = (y>0)? Int_t((y-fSensOffset)/fDpy)+1 : Int_t((y+fSensOffset)/fDpy)-1;
}

void AliMUONsegmentationV1::
GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    x = (ix>0) ? (Float_t(ix)-0.5)*fDpx+fSensOffset : (Float_t(ix)+0.5)*fDpx-fSensOffset;
    y = (iy>0) ? (Float_t(iy)-0.5)*fDpy+fSensOffset : (Float_t(iy)+0.5)*fDpy-fSensOffset;
}

void AliMUONsegmentationV1::AddCut(Int_t Zone, Int_t nX, Int_t nY)
{// the pad nX,nY is last INSIDE zone Zone
if (Zone+1>=fNzone) // no cut for last Zone : it is the natural boundary of the chamber
  printf("AliMUONsegmentationV1::AddCut ==> Zone %d not allowed !\n",Zone);
fZoneX[Zone][fNZoneCut[Zone]] = nX;
fZoneY[Zone][fNZoneCut[Zone]] = nY;
fNZoneCut[Zone]++;
}

Int_t AliMUONsegmentationV1::GetZone(Float_t X, Float_t Y)
{
Int_t iX, iY;
GetPadIxy(X,Y,iX,iY);
return GetZone( iX , iY );
}

Int_t AliMUONsegmentationV1::GetZone(Int_t nX, Int_t nY)
{// Beware : first pad begins at 1 !!
Int_t aX =  TMath::Abs(nX);
Int_t aY =  TMath::Abs(nY);
Int_t zone=fNzone-1;
for (Int_t iZone=fNzone-2;iZone>=0;iZone--) 
  {
  for (Int_t iCut=0;iCut<fNZoneCut[iZone];iCut++)
    if ( aY<=fZoneY[iZone][iCut] && aX<=fZoneX[iZone][iCut] )
      {
      zone=iZone;
      break;
      } 
  }
return zone;
}


void AliMUONsegmentationV1::SetPadCoord(Int_t iX, Int_t iY)
{    
GetPadCxy(iX,iY,fx,fy);
Float_t radius2;
if ( ( (radius2=fx*fx+fy*fy) > frSensMax2 || radius2 < frSensMin2 ) 
     && MorePads() )
  NextPad();

}

void AliMUONsegmentationV1::FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{
    //
    // Find the wire position (center of charge distribution)
    Float_t x0a=GetAnod(xhit);
    //
    // and take fNsigma*sigma around this center
    Float_t x01=x0a  - dx;
    Float_t x02=x0a  + dx;
    Float_t y01=yhit - dy;
    Float_t y02=yhit + dy;

    // Do not cross over frames...
    if (x01 * x0a < 0) 
      x01 = TMath::Sign(fSensOffset, x0a);
    if (x02 * x0a < 0) 
      x02 = TMath::Sign(fSensOffset, x0a);
    if (y01 * yhit < 0) 
      y01 = TMath::Sign(fSensOffset, yhit);
    if (y02 * yhit < 0) 
      y02 = TMath::Sign(fSensOffset, yhit);
    //
    // find the pads over which the charge distributes
    GetPadIxy(x01,y01,fixmin,fiymin);
    GetPadIxy(x02,y02,fixmax,fiymax);
    
//  printf("\n FirstPad called");
//  printf("\n Hit Position %f %f",xhit,yhit);
//  printf("\n Integration limits: %i %i %i %i",fixmin,fixmax,fiymin,fiymax);
//  printf("\n Integration limits: %f %f %f %f \n",x01,x02,y01,y02);
    // 
    // Set current pad to lower left corner
    fix=fixmin;
    fiy=fiymin;
    SetPadCoord(fix,fiy);
}

void AliMUONsegmentationV1::NextPad()
{
  // 
  // Step to next pad in integration region
    if (fix != fixmax) {
	fix++;
    } else if (fiy != fiymax) {
	fix=fixmin;
	fiy++;
    } else 
	printf("\n Error: Stepping outside integration region\n ");
    SetPadCoord(fix,fiy);
}

Int_t AliMUONsegmentationV1::MorePads()
//
// Are there more pads in the integration region
{
    if (fix == fixmax && fiy == fiymax) {
	return 0;
    } else {
	return 1;	
    }
}

Int_t AliMUONsegmentationV1::Ix()
// returns the X number of pad which has to increment charge
// due to parallel read-out
{
Int_t wix = TMath::Abs(fix)-1;
Int_t wiy = TMath::Abs(fiy)-1;
Int_t zone = GetZone(fix,fiy);
switch (zone) {
 case 0: return fix;
 case 1:
   if ( wiy%3 !=1 && wix%6>2)
     return (fix>0)? fix-3 : fix+3 ;
   return fix;
 case 2:
   if ( (wiy%6 == 1 && wix%2 ==0) ||  (wiy%6 == 4 && wix%2 ==1))
     return fix;
   return fix>0? fix - ((wix%12)/4)*4 : fix + ((wix%12)/4)*4;
 default :
   printf("Couille dans AliMUONsegmentationV1::ix\n");
 }
return -1;
}

Int_t AliMUONsegmentationV1::ISector()
{
Int_t wix = TMath::Abs(fix)-1;
Int_t wiy = TMath::Abs(fiy)-1;
Int_t zone = GetZone(fix,fiy);
switch (zone) {
 case 0: return 0;
 case 1:
   if ( wiy%3 !=1 && wix%6>2)
     return 1 ;
   return 0;
 case 2:
   if ((wiy%6 == 1 && wix%2 ==0) ||  (wiy%6 == 4 && wix%2 ==1))
     return 0;
   return (wix%12)/4;
 default :
   printf("Couille dans AliMUONsegmentationV1::ISector\n");
 }
return -1;
}

void AliMUONsegmentationV1::SigGenInit(Float_t x,Float_t y,Float_t)
{
//
//  Initialises pad and wire position during stepping
    fxt =x;
    fyt =y;
    GetPadIxy(x,y,fixt,fiyt);
    fiwt= GetiAnod(x);

}

Int_t AliMUONsegmentationV1::SigGenCond(Float_t x,Float_t y,Float_t)
{
//
//  Signal will be generated if particle crosses pad boundary or
//  boundary between two wires. 
    Int_t ixt;
    Int_t iyt;
    GetPadIxy(x,y,ixt,iyt);
    Int_t iwt= GetiAnod(x);
    
    if ((ixt != fixt) || (iyt !=fiyt) || (iwt != fiwt)) {
	return 1;
    } else {
	return 0;
    }
}

void AliMUONsegmentationV1::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
    x1=fxt-fx-fDpx/2.;
    x2=x1+fDpx;
    y1=fyt-fy-fDpy/2.;
    y2=y1+fDpy;    
}

void AliMUONsegmentationV1::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])
{
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
}

void AliMUONsegmentationV1::
FitXY(AliMUONRecCluster* ,TClonesArray* )
    // Default : Centre of gravity method
{
    ;
}
