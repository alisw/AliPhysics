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


#include "AliRICHSegResV0.h"
#include "AliRun.h"
#include "TParticle.h"
#include "TMath.h"
#include "TRandom.h"
#include "TArc.h"


ClassImp(AliRICHSegmentation)
ClassImp(AliRICHResponse)
ClassImp(AliRICHGeometry)
//___________________________________________
ClassImp(AliRICHSegmentationV0)

void AliRICHSegmentationV0::Init(AliRICHChamber* Chamber)
{
    //fNpx=(Int_t) (Chamber->ROuter()/fDpx+1);
    //fNpy=(Int_t) (Chamber->ROuter()/fDpy+1);
  fNpx=160;
  fNpy=144;
  //fNpx=80;
  //fNpy=48;
  fSector=-1;
}


Float_t AliRICHSegmentationV0::GetAnod(Float_t xhit)
{
    Float_t wire= (xhit>0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
    return fWireD*wire;
}

void AliRICHSegmentationV0::SetPadSize(Float_t p1, Float_t p2)
{
    fDpx=p1;
    fDpy=p2;
}
void AliRICHSegmentationV0::GetPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
// Please check origin of pad numbering !!!

  
  ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx);
  iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy);
  if (iy >  fNpy) iy= fNpy;
  if (iy < -fNpy) iy=-fNpy;
  if (ix >  fNpx) ix= fNpx;
  if (ix < -fNpx) ix=-fNpx;
}
void AliRICHSegmentationV0::
GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
//

  x = (ix>0) ? Float_t(ix*fDpx)-fDpx/2. : Float_t(ix*fDpx)-fDpx/2.;
  y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)-fDpy/2.;
}

void AliRICHSegmentationV0::
SetHit(Float_t xhit, Float_t yhit)
{
//
// Find the wire position (center of charge distribution)
//    Float_t x0a=GetAnod(xhit);
    fxhit=xhit;
    fyhit=yhit;
}

void AliRICHSegmentationV0::
SetPad(Int_t ix, Int_t iy)
{
    GetPadCxy(ix,iy,fx,fy);
}



void AliRICHSegmentationV0::
FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
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
    //
    // find the pads over which the charge distributes
    GetPadIxy(x01,y01,fixmin,fiymin);
    GetPadIxy(x02,y02,fixmax,fiymax);    
    // 
    // Set current pad to lower left corner
    fix=fixmin;
    fiy=fiymin;
    GetPadCxy(fix,fiy,fx,fy);
    
    //if (fSector==2)
      //printf("fix: %d, fiy: %d fx: %f, fy: %f\n",fix,fiy,fx,fy);
}

void AliRICHSegmentationV0::NextPad()
{
    //printf("\n Next Pad \n");
    
    // 
    // Step to next pad in integration region
    if (fix <= fixmax) {
//	if (fix==-1) fix++;
	fix++;
    } else if (fiy <= fiymax) {
//	if (fiy==-1) fiy++;
	fix=fixmin;
	fiy++;
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadCxy(fix,fiy,fx,fy);
}

Int_t AliRICHSegmentationV0::MorePads()

//
// Are there more pads in the integration region
{
    //printf("\n More  Pads ? \n");
  
  
  if (fix >= fixmax && fiy >= fiymax) {
    //printf("There are no more pads\n\n\n\n\n");
    return 0;
  } else {
    //printf("There are more pads\n\n");
    return 1;
  }
}

void AliRICHSegmentationV0::SigGenInit(Float_t x,Float_t y,Float_t)
{
//
//  Initialises pad and wire position during stepping
    fxt =x;
    fyt =y;
    GetPadIxy(x,y,fixt,fiyt);
    fiwt= (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}

Int_t AliRICHSegmentationV0::SigGenCond(Float_t x,Float_t y,Float_t)
{
//
//  Signal will be generated if particle crosses pad boundary or
//  boundary between two wires. 
    Int_t ixt, iyt;
    GetPadIxy(x,y,ixt,iyt);
    Int_t iwt=(x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1;
    
    if ((ixt != fixt) || (iyt !=fiyt) || (iwt != fiwt)) {
      return 1;
    } else {
      return 0;
    }
}
void AliRICHSegmentationV0::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
/*
  x1=fxt-fx-fDpx/2.;
  x2=x1+fDpx;
  y1=fyt-fy-fDpy/2.;
  y2=y1+fDpy;    
*/
  x1=fxhit-fx-fDpx/2.;
  x2=x1+fDpx;
  y1=fyhit-fy-fDpy/2.;
  y2=y1+fDpy;
}

void AliRICHSegmentationV0::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])
{
//Is used for the cluster finder, include diagonal elements
    
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
/*
    *Nlist=8;
    Xlist[0]=Xlist[1]=iX;
    Xlist[2]=iX-1;
    Xlist[3]=iX+1;
    Ylist[0]=iY-1;
    Ylist[1]=iY+1;
    Ylist[2]=Ylist[3]=iY;

   // Diagonal elements
    Xlist[4]=iX+1;
    Ylist[4]=iY+1;

    Xlist[5]=iX-1;
    Ylist[5]=iY-1;

    Xlist[6]=iX-1;
    Ylist[6]=iY+1;

    Xlist[7]=iX+1;
    Ylist[7]=iY-1;
*/
}

Float_t AliRICHSegmentationV0::Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y
, Int_t *dummy)
// Returns the square of the distance between 1 pad
// labelled by its Channel numbers and a coordinate
{
  Float_t x,y;
  GetPadCxy(iX,iY,x,y);
  return (x-X)*(x-X) + (y-Y)*(y-Y);
}


void  AliRICHSegmentationV0::GiveTestPoints(Int_t &n, Float_t *x, Float_t *y)
{
    n=1;
    x[0]=0.;
    y[0]=x[0];
}

void AliRICHSegmentationV0::Draw()
{
/*
    TArc *circle;
    Float_t scale=0.95/fRmax/2.;
    

    circle = new TArc(0.5,0.5,fRmax*scale,0.,360.);
    circle->SetFillColor(2);
    circle->Draw();

    circle = new TArc(0.5,0.5,fRmin*scale,0.,360.);
    circle->SetFillColor(1);
    circle->Draw();
*/
    ;
    
}


//___________________________________________
ClassImp(AliRICHResponseV0)

Float_t AliRICHResponseV0::IntPH(Float_t eloss)
{
    // Get number of electrons and return charge
    
    Int_t nel;
    nel= Int_t(eloss/fEIonisation);
    
    Float_t charge=0;
    if (nel == 0) nel=1;
    for (Int_t i=1;i<=nel;i++) {
	charge -= fChargeSlope*TMath::Log(gRandom->Rndm());    
    }
    return charge;
}

Float_t AliRICHResponseV0::IntPH()
{
    Float_t charge = -fChargeSlope*TMath::Log(gRandom->Rndm());
    return charge;
}



// -------------------------------------------
Float_t AliRICHResponseV0::IntXY(AliRICHSegmentation * segmentation)
{
    
    const Float_t invpitch = 1/fPitch;
    Float_t response;
//
//  Integration limits defined by segmentation model
//  
    
    Float_t xi1, xi2, yi1, yi2;
    segmentation->IntegrationLimits(xi1,xi2,yi1,yi2);

    xi1=xi1*invpitch;
    xi2=xi2*invpitch;
    yi1=yi1*invpitch;
    yi2=yi2*invpitch;

    //printf("Integration Limits: %f-%f, %f-%f\n",xi1,xi2,yi1,yi2);
    
    //printf("Invpitch:%f\n",invpitch);

    //
// The Mathieson function 
    Double_t ux1=fSqrtKx3*TMath::TanH(fKx2*xi1);
    Double_t ux2=fSqrtKx3*TMath::TanH(fKx2*xi2);
    
    Double_t uy1=fSqrtKy3*TMath::TanH(fKy2*yi1);
    Double_t uy2=fSqrtKy3*TMath::TanH(fKy2*yi2);

    //printf("Integration Data: %f-%f, %f-%f\n",ux1,ux2,uy1,uy2);
    
    //printf("%f %f %f %f\n",fSqrtKx3,fKx2,fKy4,fKx4);
    
    response=4.*fKx4*(TMath::ATan(ux2)-TMath::ATan(ux1))*fKy4*(TMath::ATan(uy2)-TMath::ATan(uy1));

    //printf("Response:%f\n",response);

    return response;       
    
}

Int_t AliRICHResponseV0::FeedBackPhotons(Float_t *source, Float_t qtot)
{
  //
  // Generate FeedBack photons
  //
  Int_t j, ipart, nt;
    
  Int_t sNfeed=0;
  
  
  // Local variables 
  Float_t cthf, ranf[2], phif, enfp = 0, sthf;
  Int_t i, ifeed;
  Float_t e1[3], e2[3], e3[3];
  Float_t vmod, uswop;
  Float_t fp, random;
  Float_t dir[3], phi;
  Int_t nfp;
  Float_t pol[3], mom[3];
  TLorentzVector position;
  //
  // Determine number of feedback photons

  //  Get weight of current particle
  TParticle *current = (TParticle*) 
    (*gAlice->Particles())[gAlice->CurrentTrack()];
    
  ifeed = Int_t(current->GetWeight()/100+0.5);
  ipart = gMC->TrackPid();
  fp = fAlphaFeedback * qtot;
  nfp = gRandom->Poisson(fp);
  
  // This call to fill the time of flight
  gMC->TrackPosition(position);
  //
  // Generate photons
  for (i = 0; i <nfp; i++) {
	
    // Direction
    gMC->Rndm(ranf, 2);
    cthf = ranf[0] * 2 - 1.;
    if (cthf < 0)  continue;
    sthf = TMath::Sqrt((1 - cthf) * (1 + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    //
    gMC->Rndm(&random, 1);
    if (random <= .57) {
      enfp = 7.5e-9;
    } else if (random <= .7) {
      enfp = 6.4e-9;
    } else {
      enfp = 7.9e-9;
    }

    dir[0] = sthf * TMath::Sin(phif);
    dir[1] = cthf;
    dir[2] = sthf * TMath::Cos(phif);
    gMC->Gdtom(dir, mom, 2);
    mom[0]*=enfp;
    mom[1]*=enfp;
    mom[2]*=enfp;
    
    // Polarisation
    e1[0] = 0;
    e1[1] = -dir[2];
    e1[2] = dir[1];
    
    e2[0] = -dir[1];
    e2[1] = dir[0];
    e2[2] = 0;
    
    e3[0] = dir[1];
    e3[1] = 0;
    e3[2] = -dir[0];
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e1[j];
      e1[j]=e3[j];
      e3[j]=uswop;
    }
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e2[j];
      e2[j]=e3[j];
      e3[j]=uswop;
    }
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e1[j]*=vmod;
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e2[j]*=vmod;
    
    gMC->Rndm(ranf, 1);
    phi = ranf[0] * 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    gMC->Gdtom(pol, pol, 2);
    
    // Put photon on the stack and label it as feedback (51, 52) 
    ++sNfeed;

    gAlice->SetTrack(Int_t(1), gAlice->CurrentTrack(), Int_t(50000051),
		     mom,source,pol,position[3],
		     "Feedback", nt, 1.);
  }
  return(sNfeed);
}

//___________________________________________
ClassImp(AliRICHGeometryV0)
