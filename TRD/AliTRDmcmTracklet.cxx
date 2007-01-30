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

///////////////////////////////////////////////////////
//                                                   //
//                                                   //
//  Tracklet object created in the local tracking    //
//                                                   //
//                                                   //
///////////////////////////////////////////////////////

#include <TGraph.h>
#include <TMath.h>
#include <TF1.h>

#include "AliLog.h"

#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDmcmTracklet.h"

ClassImp(AliTRDmcmTracklet)

//_____________________________________________________________________________
AliTRDmcmTracklet::AliTRDmcmTracklet() 
  :TObject()
  ,fDetector(-1)
  ,fRow(-1)
  ,fTrackLabel(-1)
  ,fNclusters(0)
  ,fN(0)
  ,fGPos(0)
  ,fGAmp(0)
  ,fTime0(0)
  ,fRowz(0)
  ,fSlope(0)
  ,fOffset(0)
  ,fPt(0)
  ,fdQdl(0)
{
  //
  // AliTRDmcmTracklet default constructor
  //

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

}

//_____________________________________________________________________________
AliTRDmcmTracklet::AliTRDmcmTracklet(Int_t det, Int_t row, Int_t n) 
  :TObject()
  ,fDetector(det)
  ,fRow(row)
  ,fTrackLabel(-1)
  ,fNclusters(0)
  ,fN(n)
  ,fGPos(0)
  ,fGAmp(0)
  ,fTime0(0)
  ,fRowz(0)
  ,fSlope(0)
  ,fOffset(0)
  ,fPt(0)
  ,fdQdl(0)
{
  //
  // AliTRDmcmTracklet default constructor
  //

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

  fGPos = new TGraph(0);
  fGAmp = new TGraph(0);

}

//_____________________________________________________________________________
AliTRDmcmTracklet::AliTRDmcmTracklet(const AliTRDmcmTracklet &t) 
  :TObject(t)
  ,fDetector(t.fDetector)
  ,fRow(t.fRow)
  ,fTrackLabel(t.fTrackLabel)
  ,fNclusters(t.fNclusters)
  ,fN(t.fN)
  ,fGPos(NULL)
  ,fGAmp(NULL)
  ,fTime0(t.fTime0)
  ,fRowz(t.fRowz)
  ,fSlope(t.fSlope)
  ,fOffset(t.fOffset)
  ,fPt(t.fPt)
  ,fdQdl(t.fdQdl)
{
  //
  // AliTRDmcmTracklet copy constructor
  //

  for (Int_t time = 0; time < kNtimeBins; time++) {
    for (Int_t icl = 0; icl < kNclsPads; icl++) {
      ((AliTRDmcmTracklet &) t).fADC[time][icl] = 0;
    }
    for (Int_t it = 0; it < kNdict; it++) {
      ((AliTRDmcmTracklet &) t).fTrack[time][it] = -1;
    }
    ((AliTRDmcmTracklet &) t).fTime[time]   = 0;
    ((AliTRDmcmTracklet &) t).fCol[time]    = 0;
  }

}

//_____________________________________________________________________________
AliTRDmcmTracklet::~AliTRDmcmTracklet() 
{
  //
  // AliTRDmcmTracklet destructor
  //

  if (fGPos != 0) delete fGPos;
  if (fGAmp != 0) delete fGAmp;

}

//_____________________________________________________________________________
AliTRDmcmTracklet &AliTRDmcmTracklet::operator=(const AliTRDmcmTracklet &t)
{
  //
  // Assignment operator
  //

  if (this != &t) ((AliTRDmcmTracklet &) t).Copy(*this); 
  return *this;

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::Copy(TObject &t) const
{
  //
  // Copy function
  //

  ((AliTRDmcmTracklet &) t).fDetector    = fDetector;
  ((AliTRDmcmTracklet &) t).fRow         = fRow;
  ((AliTRDmcmTracklet &) t).fTrackLabel  = fTrackLabel;
  ((AliTRDmcmTracklet &) t).fNclusters   = fNclusters;
  ((AliTRDmcmTracklet &) t).fN           = fN;
  ((AliTRDmcmTracklet &) t).fGPos        = NULL;
  ((AliTRDmcmTracklet &) t).fGAmp        = NULL;
  ((AliTRDmcmTracklet &) t).fTime0       = fTime0;
  ((AliTRDmcmTracklet &) t).fRowz        = fRowz;
  ((AliTRDmcmTracklet &) t).fSlope       = fSlope;
  ((AliTRDmcmTracklet &) t).fOffset      = fOffset;
  ((AliTRDmcmTracklet &) t).fPt          = fPt;
  ((AliTRDmcmTracklet &) t).fdQdl        = fdQdl;

  for (Int_t time = 0; time < kNtimeBins; time++) {
    for (Int_t icl = 0; icl < kNclsPads; icl++) {
      ((AliTRDmcmTracklet &) t).fADC[time][icl] = 0;
    }
    for (Int_t it = 0; it < kNdict; it++) {
      ((AliTRDmcmTracklet &) t).fTrack[time][it] = -1;
    }
    ((AliTRDmcmTracklet &) t).fTime[time]   = 0;
    ((AliTRDmcmTracklet &) t).fCol[time]    = 0;
  }

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::Reset() 
{
  //
  // Reset the tracklet information
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

  fGPos->Set(0);
  fGAmp->Set(0);

  fSlope  = 0.0;
  fOffset = 0.0;
  fTime0  = 0.0;
  fRowz   = 0.0;
  fPt     = 0.0;
  fdQdl   = 0.0;

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::AddCluster(Int_t icol, Int_t itb, Float_t *adc, Int_t *track) 
{
  //
  // Add a cluster to the tracklet
  //
 
  if (fNclusters >= kNtimeBins) return;

  for (Int_t icl = 0; icl < kNclsPads; icl++) {
    fADC[fNclusters][icl] = adc[icl]; 
  }

  fTrack[fNclusters][0] = track[0];
  fTrack[fNclusters][1] = track[1];
  fTrack[fNclusters][2] = track[2];
  fTime[fNclusters] = itb;
  fCol[fNclusters] = icol;

  fNclusters++;

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::MakeTrackletGraph(AliTRDgeometry *geo, Float_t field) 
{
  //
  // Tracklet graph of positions (global coordinates, rotated [cm])
  //
  
  if (!geo) {
    AliError("No geometry.");
    return;
  }
  
  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliError("No common parameters.");
    return;
  }

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No instance of AliTRDcalibDB.");
    return;
  }

  Int_t iplan, icham;

  iplan = geo->GetPlane(fDetector);
  icham = geo->GetChamber(fDetector);

  AliTRDpadPlane *padPlane = commonParam->GetPadPlane(iplan,icham);

  Float_t samplFreq = commonParam->GetSamplingFrequency();

  Int_t time, col;
  Float_t amp[3];
  Float_t xpos, ypos, xzero, colSize, timeBinSize;
  Float_t vDrift, omegaTau, lorentzAngle, thetaSlope;
  Float_t tiltingAngle;
  Int_t npg = 0;

  for (Int_t icl = 0; icl < fNclusters; icl++) {

    time   = GetClusterTime(icl);

    amp[0] = GetClusterADC(icl)[0];
    amp[1] = GetClusterADC(icl)[1];
    amp[2] = GetClusterADC(icl)[2];

    col    = GetClusterCol(icl);

    if (amp[0] < 0.0 || amp[1] < 0.0 || amp[2] < 0.0) continue;

    ypos = GetClusY(amp,iplan);

    colSize = padPlane->GetColSize(col);
    vDrift = calibration->GetVdrift(fDetector,col,fRow);
    timeBinSize = vDrift/samplFreq;
    
    // From v4-03-Release to HEAD28Mar06 the sign has changed from "-" to "+" 
    // due to a change in the digitizer
    omegaTau     = TMath::Sign(1.0,(Double_t)field)*GetOmegaTau(vDrift,TMath::Abs(field));
    lorentzAngle = TMath::ATan(omegaTau)*180.0/TMath::Pi();
    
    xpos = (time+0.5) * timeBinSize;
    xpos = geo->GetTime0(iplan) - xpos;

    ypos = padPlane->GetColPos(col) - (ypos + 0.5) * colSize;

    // ExB correction
    xzero = geo->GetTime0(iplan);
    ypos = ypos + (xpos-xzero) * omegaTau;

    // tilted pads correction
    thetaSlope = - padPlane->GetRowPos(fRow)/geo->GetTime0(iplan);
    tiltingAngle = padPlane->GetTiltingAngle()/180.0*TMath::Pi();
    ypos = ypos - (xpos-xzero) * thetaSlope * TMath::Sin(tiltingAngle);

    fGPos->SetPoint(npg,(Double_t)xpos,(Double_t)ypos);
    npg++;

  }

  fGPos->Set(npg);
  
  fTime0 = geo->GetTime0(iplan) - AliTRDgeometry::CdrHght() - 0.5*AliTRDgeometry::CamHght();
  fRowz = padPlane->GetRowPos(fRow) - padPlane->GetRowSize(fRow)/2.0;

  Double_t xMin = 0;
  Double_t xMax = 0;
  Double_t x, y;
  fGPos->GetPoint(0    ,x,y); 
  xMax = x + 0.1;
  fGPos->GetPoint(npg-1,x,y); 
  xMin = x - 0.1;
  
  TF1 *line = new TF1("line","[0]+x*[1]",xMin,xMax);
  fGPos->Fit(line,"WRQ0");

  fOffset = line->Eval(fTime0);
  fSlope  = TMath::ATan(line->GetParameter(1))*180.0/TMath::Pi();

  line->Delete();
  
  Float_t fx = fTime0;
  Float_t fy = fOffset;
  
  Float_t infSlope = TMath::ATan(fy/fx)/TMath::Pi()*180.0;    
  Float_t alpha    = fSlope - infSlope;
  Float_t r        = TMath::Sqrt(fx*fx + fy*fy)/(2.0*TMath::Sin(alpha/180.0*TMath::Pi()));

  fPt = 0.3 * field * 0.01 * r;
  
  return;

}

//_____________________________________________________________________________
void AliTRDmcmTracklet::MakeClusAmpGraph() 
{
  //
  // Tracklet graph of cluster charges
  //

  Int_t   time;
  Float_t amp[3];
  Int_t   npg = 0;

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
Float_t AliTRDmcmTracklet::GetClusY(Float_t *adc, Int_t pla) const
{
  //
  // Cluster position in the phi direction in pad units (relative to the pad border)
  //

  Float_t ypos = 0.0;

  Float_t a0 = adc[0];
  Float_t a1 = adc[1];
  Float_t a2 = adc[2];

  Float_t t1, t2, w1 ,w2;

  Float_t w = 1.0;  // pad units

  Float_t sigma = 0.0;

  switch(pla) {
  case 0:
    sigma = 0.515; break;
  case 1:
    sigma = 0.501; break;
  case 2:
    sigma = 0.491; break;
  case 3:
    sigma = 0.481; break;
  case 4:
    sigma = 0.471; break;
  case 5:
    sigma = 0.463; break;
  default:
    AliError("Wrong plane number.");
    return 0.0;
  }

  sigma *= w;

  t1 = 0.0;
  w1 = 0.0;
  if( a0 > 0 ) {
    w1 = a0*a0;
    //w1 = a0;
    t1 = w1*((sigma*sigma)/w*TMath::Log(a1/a0)-0.5*w);
  }
  t2 = 0.0;
  w2 = 0.0;
  if( a2 > 0 ) {
    w2 = a2*a2;
    //w2 = a2;
    t2 = w2*((sigma*sigma)/w*TMath::Log(a2/a1)+0.5*w);
  }

  ypos = w*(t1+t2)/(w1+w2);  // range: -0.5*w ... +0.5*w

  return ypos;

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
	  AliWarning("Too many tracks for this tracklet.");
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
Float_t AliTRDmcmTracklet::GetOmegaTau(Float_t vdrift, Float_t field) const
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
