/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version LYON //
/////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

//#include "AliMUONv0.h"
#include "AliMUONchamber.h"
#include "AliMUONSegResV1.h"
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
  // beware : frMin is SENSITIVE radius by definition.
  frSensMin2 = (Chamber->RInner())*(Chamber->RInner());
  frSensMax2 = (Chamber->ROuter())*(Chamber->ROuter());
  fNpx=(Int_t) (Chamber->ROuter()/fDpx) + 1;
  fNpy=(Int_t) (Chamber->ROuter()/fDpy) + 1;
  //    fNwire=3;
  DefaultCut();
  fCorr=0;
}

void AliMUONsegmentationV1::DefaultCut(void)
{
SetNzone(3);
AddCut(0,5*6,18*8);
AddCut(0,9*6,15*8);
AddCut(0,11*6,12*8);
AddCut(0,12*6,9*8);
AddCut(0,13*6,6*8);
AddCut(1,6*6,20*12);
AddCut(1,12*6,18*12);
AddCut(1,15*6,15*12);
AddCut(1,18*6,12*12);
AddCut(1,21*6,9*12);
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
{// the pad nX,nY is last INSIDE zone Zone. First pad is labelled 1 and not 0
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

void AliMUONsegmentationV1::
SetHit(Float_t xhit, Float_t yhit)
{
    //
    // Find the wire position (center of charge distribution)
//    Float_t x0a=GetAnod(xhit);
    fxhit=xhit;
    fyhit=yhit;
}

void AliMUONsegmentationV1::
SetPad(Int_t ix, Int_t iy)
{
    GetPadCxy(ix,iy,fx,fy);
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
    fxhit=x0a;
    fyhit=yhit;
    
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

Int_t AliMUONsegmentationV1::IsParallel2(Int_t iX, Int_t iY)
// test if the pad is read in parallel for zone 2
// iX and iY are assumed to be positive and starting at 0 numbering (cF. iX)
// returns 1 or 2 if read in parallel, 
// according to the actual number in the chain, 0 else
//
// chainage is        result is 
// 1  2  3  1  2  3   1 1 1 2 2 2     y
// 7  8  9 10 11 12   0 0 0 0 0 0     ^
// 4  5  6  4  5  6   1 1 1 2 2 2     +->x
//
{
if (iY%3==1) return 0;
return (iX%6)/3+1;
}

Int_t AliMUONsegmentationV1::IsParallel3(Int_t iX, Int_t iY)
// test if the pad is read in parallel for zone 3
// iX and iY are assumed to be positive and starting at 0 numbering (cF. iX)
// returns 1,2 or 3 if read in parallel, 
// according to the actual number in the chain, 0 else
//
// chainage is                 result is
//16  2  3  1  2  3  1  2  3   0 1 1 1 2 2 2 3 3
// 7  8  9 10 11 12 13 14 15   0 0 0 0 0 0 0 0 0
// 4  5  6  4  5  6  4  5  6   1 1 1 2 2 2 3 3 3
//
{
if (iY%3==1) return 0;
return (iX%9)/3+1 - (iY%3==2 && iX%3==0);
}

Int_t AliMUONsegmentationV1::NParallel2(Int_t , Int_t iY)
// returns the number of pads connected in parallel for zone 2
// iX and iY are assumed to be positive and starting at 0 numbering (cF. iX)
//
// result is
// 2 2 2 2 2 2
// 1 1 1 1 1 1
// 2 2 2 2 2 2
//
{
if (iY%3==1) return 1;
return 2;
}

Int_t AliMUONsegmentationV1::NParallel3(Int_t iX, Int_t iY)
// test if the pad is read in parallel for zone 3
// iX and iY are assumed to be positive and starting at 0 numbering (cF. iX)
// returns 1,2 or 3 if read in parallel, 
// according to the actual number in the chain, 0 else
//
// result is 
// 1 3 3 2 3 3 2 3 3 
// 1 1 1 1 1 1 1 1 1
// 3 3 3 3 3 3 3 3 3
//
{
if (iY%3==1) return 1;
if (iY%3==2 && iX%9==0) return 1;
return 3 - (iY%3==2 && iX%3==0);
}


Int_t AliMUONsegmentationV1::Ix(Int_t trueX, Int_t trueY)
// returns the X number of pad which corresponds to the logical 
// channel, expressed in x and y.
{
Int_t wix = TMath::Abs(trueX)-1;
Int_t wiy = TMath::Abs(trueY)-1;
Int_t zone = GetZone(trueX,trueY);
Int_t par3;
switch (zone) {
 case 0: return trueX;
 case 1:
   if (IsParallel2(wix,wiy) == 2)
     return (trueX>0)? trueX-3 : trueX+3 ;
   return trueX;
 case 2:
   if ( (par3= IsParallel3(wix,wiy)) )
     return (trueX>0) ? trueX-3*(par3-1) : trueX+3*(par3-1) ;
   return trueX ;
 default :
   printf("Couille dans AliMUONsegmentationV1::ix\n");
 }
return -1;
}

Int_t AliMUONsegmentationV1::Ix()
// returns the X number of pad which has to increment charge
// due to parallel read-out
{return Ix(fix,fiy);}

Int_t AliMUONsegmentationV1::ISector()
// This function is of no use for this kind of segmentation.
{
return GetZone(fix,fiy);
}

void AliMUONsegmentationV1::SigGenInit(Float_t x,Float_t y,Float_t )
{
//
//  Initialises pad and wire position during stepping
    fxt =x;
    fyt =y;
    GetPadIxy(x,y,fixt,fiyt);
    fiwt= GetiAnod(x);

}

Int_t AliMUONsegmentationV1::SigGenCond(Float_t x,Float_t y,Float_t )
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
    x1=fxhit-fx-fDpx/2.;
    x2=x1+fDpx;
    y1=fyhit-fy-fDpy/2.;
    y2=y1+fDpy;    
}

void AliMUONsegmentationV1::GetNParallelAndOffset(Int_t iX, Int_t iY,Int_t
*Nparallel, Int_t* Offset)
{
Int_t wix = TMath::Abs(iX)-1;
Int_t wiy = TMath::Abs(iY)-1;
Int_t zone = GetZone(iX,iY);
switch (zone) {
 case 0: 
    *Nparallel=1;
    *Offset=0;
    break;
 case 1:
    *Nparallel = NParallel2(wix,wiy);
    (iX>0) ? *Offset =3 : *Offset = -3;
    if (IsParallel2(wix,wiy)>1)
      printf("GetNParallelAndOffset called for existing channel -> answer is crazy\n");
    break;
 case 2:
    *Nparallel = NParallel3(wix,wiy);
    (iX>0) ? *Offset =3 : *Offset = -3;
    if (IsParallel3(wix,wiy)>1)
      printf("GetNParallelAndOffset called for existing channel -> answer is crazy\n");
    break;
 }
}


Float_t AliMUONsegmentationV1::Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *Offset)
//
// Computes the offset for which the physical pad has the minimum distance squared
// (returned value) to the given coordinates
{
Int_t nPara,offset;
GetNParallelAndOffset(iX,iY,&nPara,&offset);
Float_t d2min=1E10;
for (Int_t i=0;i<nPara; i++)
  {
  Float_t x,y;
  GetPadCxy(iX+i*offset,iY,x,y);
  Float_t d2=(x-X)*(x-X) + (y-Y)*(y-Y);
  if ( d2min > d2)
    {
    d2min = d2;
    *Offset = i*offset;
    }
  }
return d2min; 
}

void AliMUONsegmentationV1::CleanNeighbours(Int_t* Nlist, Int_t *Xlist, 
 	Int_t *Ylist)
// In the raw neighbours list, some pads do not exist 
// and some others are read in parallel ...
// So we prune non-existing neighbours from the list (event if this should be
// at last not be a problem due to the clustering algorithm...)
{
Int_t nTot=0;
for (Int_t nList=0;nList<*Nlist;nList++)
  {
  // prune if it does not exist
  if ( Xlist[nList]==0 || Ylist[nList]==0 )
    continue;
  // compute true position
  Xlist[nTot] = Ix(Xlist[nList],Ylist[nList]) ;
  Ylist[nTot] = Ylist[nList] ;
  // and prune if it does already exist
  Int_t nTest;
  for (nTest=0;nTest<nTot; nTest++)
    {
    if ( Xlist[nTest]==Xlist[nTot] && Ylist[nTest]==Ylist[nTot])
      // we found it
      break ;
    }
  if (nTest==nTot)
    nTot++;
  }
*Nlist = nTot;
}

void AliMUONsegmentationV1::
NeighboursNonDiag(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[12], Int_t Ylist[12])
// returns the X number of pad which has to increment charge
// due to parallel read-out
{
Int_t nParallel, offset;
GetNParallelAndOffset(iX,iY,&nParallel,&offset);
//
// now fill raw list of neighbours
*Nlist=4*nParallel;
Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
if (nParallel>1) {
  Xlist[4]=Xlist[5]=iX+offset;Xlist[6]=iX+offset-1;Xlist[7]=iX+offset+1;
  Ylist[4]=iY-1;Ylist[5]=iY+1;Ylist[6]=Ylist[7]=iY;
  if (nParallel>2) {
    Xlist[8]=Xlist[9]=iX+2*offset;Xlist[10]=iX+2*offset-1;Xlist[11]=iX+2*offset+1;
    Ylist[8]=iY-1;Ylist[9]=iY+1;Ylist[10]=Ylist[11]=iY;
    }
  }
CleanNeighbours(Nlist,Xlist,Ylist);
}

void AliMUONsegmentationV1::
NeighboursDiag(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[24], Int_t Ylist[24])
// returns the X number of pad which has to increment charge
// due to parallel read-out
{
Int_t nParallel, offset;
GetNParallelAndOffset(iX,iY,&nParallel,&offset);
//
// now fill raw list of neighbours
*Nlist=0;
for (Int_t i=0;i<nParallel;i++)
  for (Int_t dx=-1;dx<2;dx++)
    for (Int_t dy=-1;dy<2;dy++)
      {
      if (dx==dy && dy==0)
        continue; 
      Xlist[*Nlist] = iX + dx + i*offset;
      Ylist[*Nlist] = iY + dy;
      (*Nlist)++;
      }
CleanNeighbours(Nlist,Xlist,Ylist);
}

void AliMUONsegmentationV1::Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, 
	Int_t Xlist[24], Int_t Ylist[24])
{NeighboursDiag(iX,iY,Nlist,Xlist,Ylist);}



void AliMUONsegmentationV1::GiveTestPoints(Int_t &n, Float_t *x, Float_t *y)
{
    n=1;
    x[0]=(TMath::Sqrt(frSensMax2)-TMath::Sqrt(frSensMin2))/2/TMath::Sqrt(2.);
    y[0]=x[0];
}

