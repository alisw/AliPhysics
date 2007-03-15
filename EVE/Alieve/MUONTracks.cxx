#include "MUONTracks.h"

#include <Reve/Track.h>

#include <AliMUONTrack.h>
#include <AliMUONTrackParam.h>

#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixD.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// MUONTracks
// Produce Reve:Track from AliMUONTrack with dipole field model

ClassImp(MUONTracks)

//______________________________________________________________________
MUONTracks::MUONTracks()
{
  //
  // constructor
  //

}

//______________________________________________________________________
MUONTracks::~MUONTracks()
{
  //
  // destructor
  //

}

//______________________________________________________________________
void MUONTracks::MakeTrack(Int_t label, Reve::Track *rtrack, AliMUONTrack *mtrack)
{
  //
  // 
  //

  Double_t xv, yv;
  Float_t ax, bx, ay, by;
  Float_t xr[20], yr[20], zr[20];
  Float_t xrc[20], yrc[20], zrc[20];
  char form[1000];
    
  TMatrixD smatrix(2,2);
  TMatrixD sums(2,1);
  TMatrixD res(2,1);

  Float_t xRec, xRec0;
  Float_t yRec, yRec0;
  Float_t zRec, zRec0;
  
  Float_t zg[4] = { -1603.5, -1620.5, -1703.5, -1720.5 };

  Int_t count = 0;

  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  if (mtrack->GetMatchTrigger()) {
    //PH      track->SetName(Form("MUONTrack %2d (MT)", label));
    sprintf(form,"MUONTrack %2d (MT)", label);
    rtrack->SetName(form);
    rtrack->SetLineColor(7);
  } else {
    //PH      track->SetName(Form("MUONTrack %2d     ", label));
    sprintf(form,"MUONTrack %2d     ", label);
    rtrack->SetName(form);
    rtrack->SetLineColor(6);
  }
  
  AliMUONTrackParam *trackParam = mtrack->GetTrackParamAtVertex(); 
  xRec0  = trackParam->GetNonBendingCoor();
  yRec0  = trackParam->GetBendingCoor();
  zRec0  = trackParam->GetZ();
  
  rtrack->SetPoint(count,xRec0,yRec0,zRec0);
  count++;
  
  for (Int_t i = 0; i < 10; i++) xr[i]=yr[i]=zr[i]=0.0;
  
  Int_t nTrackHits = mtrack->GetNTrackHits();
  //printf("Nhits = %d \n",nTrackHits);
  TClonesArray* trackParamAtHit;
  for (Int_t iHit = 0; iHit < nTrackHits; iHit++){
    trackParamAtHit = mtrack->GetTrackParamAtHit();
    trackParam = (AliMUONTrackParam*) trackParamAtHit->At(iHit); 
    xRec  = trackParam->GetNonBendingCoor();
    yRec  = trackParam->GetBendingCoor();
    zRec  = trackParam->GetZ();
    
    //printf("Hit %d x %f y %f z %f \n",iHit,xRec,yRec,zRec);
    
    xr[iHit] = xRec;
    yr[iHit] = yRec;
    zr[iHit] = zRec;
    
    rtrack->SetPoint(count,xRec,yRec,zRec);
    count++;
    
  }
  
  Int_t nrc = 0;
  if (mtrack->GetMatchTrigger() && 1) {
    
    for (Int_t i = 0; i < nTrackHits; i++) {
      if (TMath::Abs(zr[i]) > 1000.0) {
	//printf("Hit %d x %f y %f z %f \n",iHit,xr[i],yr[i],zr[i]);
	xrc[nrc] = xr[i];
	yrc[nrc] = yr[i];
	zrc[nrc] = zr[i];
	nrc++;
      }
    }
    
    if (nrc < 2) return;
    
    // fit x-z
    smatrix.Zero();
    sums.Zero();
    for (Int_t i = 0; i < nrc; i++) {
      xv = (Double_t)zrc[i];
      yv = (Double_t)xrc[i];
      //printf("x-z: xv %f yv %f \n",xv,yv);
      smatrix(0,0) += 1.0;
      smatrix(1,1) += xv*xv;
      smatrix(0,1) += xv;
      smatrix(1,0) += xv;
      sums(0,0)    += yv;
      sums(1,0)    += xv*yv;
    }
    res = smatrix.Invert() * sums;
    ax = res(0,0);
    bx = res(1,0);
    
    // fit y-z
    smatrix.Zero();
    sums.Zero();
    for (Int_t i = 0; i < nrc; i++) {
      xv = (Double_t)zrc[i];
      yv = (Double_t)yrc[i];
      //printf("y-z: xv %f yv %f \n",xv,yv);
      smatrix(0,0) += 1.0;
      smatrix(1,1) += xv*xv;
      smatrix(0,1) += xv;
      smatrix(1,0) += xv;
      sums(0,0)    += yv;
      sums(1,0)    += xv*yv;
    }
    res = smatrix.Invert() * sums;
    ay = res(0,0);
    by = res(1,0);
    
    Float_t xtc, ytc, ztc;
    for (Int_t ii = 0; ii < 4; ii++) {
      
      ztc = zg[ii];
      ytc = ay+by*zg[ii];
      xtc = ax+bx*zg[ii];
      
      //printf("tc: x %f y %f z %f \n",xtc,ytc,ztc);
      
      rtrack->SetPoint(count,xtc,ytc,ztc);
      count++;
      
    }
    
  }  // end match trigger
  
}

