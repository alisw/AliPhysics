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

#include <TGraph.h>
#include <TMath.h>
#include <TF1.h>

#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"

#include "AliTRDmcmTracklet.h"

ClassImp(AliTRDmcmTracklet)

//_____________________________________________________________________________
AliTRDmcmTracklet::AliTRDmcmTracklet() 
{

  //
  // AliTRDmcmTracklet default constructor
  //

  fDetector = -1;
  fRow      = -1;

  for (Int_t time = 0; time < kNtimeBins; time++) {
    for (Int_t icl = 0; icl < kNclsPads; icl++) {
      fADC[time][icl] = 0;
    }
    for (Int_t it = 0; it < kNdict; it++) {
      fTrack[time][it] = -1;
    }
    fTime[time]   = 0;
    fCol[time]    = 0;
  }

  fNclusters  =  0;
  fN          =  0;
  fTrackLabel = -1;

  fGPos = 0;
  fGAmp = 0;

  fSlope  = 0.0;
  fOffset = 0.0;
  fTime0  = 0.0;
  fRowz   = 0.0;
  fPt     = 0.0;
  fdQdl   = 0.0;

}

//_____________________________________________________________________________
AliTRDmcmTracklet::AliTRDmcmTracklet(Int_t det, Int_t row, Int_t n) 
{

  //
  // AliTRDmcmTracklet default constructor
  //

  fDetector = det;
  fRow      = row;

  for (Int_t time = 0; time < kNtimeBins; time++) {
    for (Int_t icl = 0; icl < kNclsPads; icl++) {
      fADC[time][icl] = 0;
    }
    for (Int_t it = 0; it < kNdict; it++) {
      fTrack[time][it] = -1;
    }
    fTime[time]   = 0;
    fCol[time]    = 0;
  }

  fNclusters = 0;

  fN = n;

  fTrackLabel = -1;

  fGPos = new TGraph(0);
  fGAmp = new TGraph(0);

  fSlope  = 0.0;
  fOffset = 0.0;
  fTime0  = 0.0;
  fRowz   = 0.0;
  fPt     = 0.0;
  fdQdl   = 0.0;

}

//_____________________________________________________________________________
AliTRDmcmTracklet::~AliTRDmcmTracklet() 
{

  //
  // AliTRDmcmTracklet destructor
  //

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::MakeTrackletGraph(AliTRDgeometry *geo, Float_t field) 
{

  //
  // Tracklet graph of positions (global coordinates, rotated [cm])
  //
  
  if (!geo) {
    Error("MakeTrackletGraph","No geometry.");
    return;
  }
  
  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    Error("MakeTrackletGraph","No common params.");
    return;
  }

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    Error("MakeTrackletGraph","No instance of AliTRDcalibDB.");
    return;
  }

  Int_t iplan, icham;

  iplan = geo->GetPlane(fDetector);
  icham = geo->GetChamber(fDetector);

  AliTRDpadPlane *padPlane = commonParam->GetPadPlane(iplan,icham);

  Float_t SamplFreq = calibration->GetSamplingFrequency();

  Int_t time, col;
  Float_t amp[3];
  Float_t Xpos, Ypos, Xzero, ColSize, TimeBinSize;
  Float_t vDrift, OmegaTau, LorentzAngle, ThetaSlope;
  Float_t TiltingAngle;
  Int_t npg = 0;

  for (Int_t icl = 0; icl < fNclusters; icl++) {

    time   = GetClusterTime(icl);

    amp[0] = GetClusterADC(icl)[0];
    amp[1] = GetClusterADC(icl)[1];
    amp[2] = GetClusterADC(icl)[2];

    col    = GetClusterCol(icl);

    if (amp[0] < 0.0 || amp[1] < 0.0 || amp[2] < 0.0) continue;

    Ypos = GetClusY(amp,iplan);

    ColSize = padPlane->GetColSize(col);
    vDrift = calibration->GetVdrift(fDetector,col,fRow);
    TimeBinSize = vDrift/SamplFreq;
    
    // From v4-03-Release to HEAD28Mar06 the sign has changed from "-" to "+" 
    // due to a change in the digitizer
    OmegaTau = +TMath::Sign(1.0,(Double_t)field)*GetOmegaTau(vDrift,TMath::Abs(field));
    LorentzAngle = TMath::ATan(OmegaTau)*180.0/TMath::Pi();
    
    Xpos = (time+0.5) * TimeBinSize;
    Xpos = geo->GetTime0(iplan) - Xpos;

    Ypos = padPlane->GetColPos(col) - (Ypos + 0.5) * ColSize;

    // ExB correction
    Xzero = geo->GetTime0(iplan);
    Ypos = Ypos + (Xpos-Xzero) * OmegaTau;

    // tilted pads correction
    ThetaSlope = - padPlane->GetRowPos(fRow)/geo->GetTime0(iplan);
    TiltingAngle = padPlane->GetTiltingAngle()/180.0*TMath::Pi();
    Ypos = Ypos - (Xpos-Xzero) * ThetaSlope * TMath::Sin(TiltingAngle);

    fGPos->SetPoint(npg,(Double_t)Xpos,(Double_t)Ypos);
    npg++;

  }

  fGPos->Set(npg);
  
  fTime0 = geo->GetTime0(iplan) - AliTRDgeometry::CdrHght() - 0.5*AliTRDgeometry::CamHght();
  fRowz = 0.5*(padPlane->GetRowPos(fRow) + padPlane->GetRowPos(fRow+1));

  Double_t xMin = 0, xMax = 0, x, y;
  fGPos->GetPoint(0    ,x,y); xMax = x + 0.1;
  fGPos->GetPoint(npg-1,x,y); xMin = x - 0.1;
  
  TF1 *line = new TF1("line","[0]+x*[1]",xMin,xMax);
  fGPos->Fit(line,"WRQ0");

  fOffset = line->Eval(fTime0);
  fSlope  = TMath::ATan(line->GetParameter(1))*180.0/TMath::Pi();

  line->Delete();
  
  Float_t fX = fTime0;
  Float_t fY = fOffset;
  
  Float_t infSlope = TMath::ATan(fY/fX)/TMath::Pi()*180.0;    
  Float_t alpha = fSlope - infSlope;
  Float_t R = TMath::Sqrt(fX*fX + fY*fY)/(2.0*TMath::Sin(alpha/180.0*TMath::Pi()));

  fPt = 0.3 * field * 0.01 * R;
  
  return;

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::MakeClusAmpGraph() 
{
  //
  // Tracklet graph of cluster charges
  //

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    Error("MakeClusAmpGraph","No instance of AliTRDcalibDB.");
    return;
  }

  Int_t time;
  Float_t amp[3];
  Int_t npg = 0;
  fdQdl = 0.0;
  for (Int_t icl = 0; icl < fNclusters; icl++) {

    time   = GetClusterTime(icl);

    amp[0] = GetClusterADC(icl)[0];
    amp[1] = GetClusterADC(icl)[1];
    amp[2] = GetClusterADC(icl)[2];

    fGAmp->SetPoint(npg,(Double_t)(time+0.5),(Double_t)(amp[0]+amp[1]+amp[2]));
    npg++;

    fdQdl += amp[0]+amp[1]+amp[2];

  }

  fGAmp->Set(npg);

  fdQdl /= (Float_t)npg;

  return;

}

//_____________________________________________________________________________
Float_t AliTRDmcmTracklet::GetClusY(Float_t *adc, Int_t pla) 
{
  //
  // Cluster position in the phi direction in pad units (relative to the pad border)
  //

  Float_t Ypos = 0.0;

  Float_t A0 = adc[0];
  Float_t A1 = adc[1];
  Float_t A2 = adc[2];

  Float_t T1, T2, W1 ,W2;

  Float_t W = 1.0;  // pad units

  Float_t Sigma = 0.0;

  switch(pla) {
  case 0:
    Sigma = 0.515; break;
  case 1:
    Sigma = 0.501; break;
  case 2:
    Sigma = 0.491; break;
  case 3:
    Sigma = 0.481; break;
  case 4:
    Sigma = 0.471; break;
  case 5:
    Sigma = 0.463; break;
  default:
    Error("GetClusY","Wrong plane number.");
    return 0.0;
  }

  Sigma *= W;

  T1 = 0.0;
  W1 = 0.0;
  if( A0 > 0 ) {
    W1 = A0*A0;
    //W1 = A0;
    T1 = W1*((Sigma*Sigma)/W*TMath::Log(A1/A0)-0.5*W);
  }
  T2 = 0.0;
  W2 = 0.0;
  if( A2 > 0 ) {
    W2 = A2*A2;
    //W2 = A2;
    T2 = W2*((Sigma*Sigma)/W*TMath::Log(A2/A1)+0.5*W);
  }

  Ypos = W*(T1+T2)/(W1+W2);  // range: -0.5*W ... +0.5*W

  return Ypos;

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::CookLabel(Float_t frac) 
{
  //
  // Cook the track label from cluster labels
  //

  const Int_t kMaxTracks = 10;
  Int_t trackLabel[kMaxTracks];
  Int_t trackCount[kMaxTracks];
  for (Int_t it = 0; it < kMaxTracks; it++) {
    trackLabel[it] = -1;
    trackCount[it] = 0;
  }

  Bool_t counted;
  Int_t label, nTracks = 0;
  for (Int_t icl = 0; icl < fNclusters; icl++) {

    for (Int_t id = 0; id < kNdict; id++) {

      if (fTrack[icl][id] == -1) continue;

      label = fTrack[icl][id];

      counted = kFALSE;
      for (Int_t it = 0; it < nTracks; it++) {
	if (label == trackLabel[it]) {
	  trackCount[it]++;
	  counted = kTRUE;
	  break;
	}
      }
      if (!counted) {
	trackLabel[nTracks] = label;
	trackCount[nTracks]++;
	nTracks++;
	if (nTracks == kMaxTracks) {
	  Warning("CookLabel","Too many tracks for this tracklet.");
	  nTracks--;
	  break;
	}
      }

    }

  }

  for (Int_t it = 0; it < kMaxTracks; it++) {
    if (trackCount[it] >= (Int_t)(frac*fNclusters)) {
      fTrackLabel = trackLabel[it];
      break;
    }
  }

}

//_____________________________________________________________________________
Float_t AliTRDmcmTracklet::GetOmegaTau(Float_t vdrift, Float_t field)
{
  //
  // Returns omega*tau (tan(Lorentz-angle)) for a given drift velocity <vd> 
  // and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  //
  // Copy of the "AliTRDcalibDB" function, taking as argument the magnetic field too
  //
  
  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.004810,  0.007412,  0.010252,  0.013409,  0.016888 };
  Float_t p1[kNb] = {  0.054875,  0.081534,  0.107333,  0.131983,  0.155455 };
  Float_t p2[kNb] = { -0.008682, -0.012896, -0.016987, -0.020880, -0.024623 };
  Float_t p3[kNb] = {  0.000155,  0.000238,  0.000330,  0.000428,  0.000541 };

  Int_t ib = ((Int_t) (10 * (field - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t alphaL = p0[ib] 
      + p1[ib] * vdrift
      + p2[ib] * vdrift*vdrift
      + p3[ib] * vdrift*vdrift*vdrift;

  return TMath::Tan(alphaL);

}
