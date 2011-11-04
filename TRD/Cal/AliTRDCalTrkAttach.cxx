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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Container for calibration parameters for AliTRDseedV1::AttachClusters()
// .... Longer description of content ...
//
// For calibration procedure check AliTRDtrackleOflHelper::CalibrateAttach()
// .... some reference to the calibration procedure .........
//                                                                      //
// Authors:                                                             //
//   Alex Bercuci <a.bercuci@gsi.de>                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include "TGraphErrors.h"
#include "TObjArray.h"

#include "AliLog.h"

#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"
#include "AliTRDCalTrkAttach.h"

ClassImp(AliTRDCalTrkAttach)

//______________________________________________________________________________
AliTRDCalTrkAttach::AliTRDCalTrkAttach()
  :TNamed("AliTRDCalTrkAttach", "Calibration of AliTRDseedV1::AttachClusters")
  ,fRClikeLimit(0.65)
  ,fScaleCov(2.)
  ,fLike(NULL)
{
// Default constructor

  fNsgmDy[0] = 5; fNsgmDy[1] = 7;
  fLikeMinRelDecrease[0] = .2; fLikeMinRelDecrease[1] = .3;
}

//______________________________________________________________________________
AliTRDCalTrkAttach::~AliTRDCalTrkAttach()
{
// Destructor
  if(fLike) delete fLike;
}


//______________________________________________________________________________
Double_t AliTRDCalTrkAttach::CookLikelihood(Bool_t chg, Int_t ly, Float_t pt, Float_t phiTrk, Int_t n, Double_t dyr, Double_t dphi, Double_t sr) const
{
// Calculate likelihood for a segment to belong to a tracklet
// Based on calibrated values

  if(n<4){
    AliDebug(4, Form("Failed basic cut[s] : n[%d] ...", n));
    return 0.;
  }
  //check likelihood array  
  if (!fLike || fLike->GetEntries() != AliTRDgeometry::kNlayer*kNcharge*kNcalib) {
    AliError("No usable AttachClusters calib object.");
    return 0.;
  }
  Int_t offset(ly*kNcalib*kNcharge); // offset for layer
  TGraphErrors *g(NULL);
  if(!(g = (TGraphErrors*)fLike->At(offset+Int_t(kResPos)*kNcharge+Int_t(chg)))){
    AliError("Failed retrieving AttachClusters graph.");
    return 0.;
  }
  // Interpolate p_t
  Int_t npts(g->GetN()), ip(0), jp(-1);
  Double_t x0, y0, x1, y1, dd(0.), invdx(0.), f[4]={0., 0., 0., 0.};
  for(Int_t kp(0); kp<npts; kp++){
    g->GetPoint(kp, x1, y1);
    if(x1>=pt){jp=kp; break;}
  }
  Bool_t boundary(kFALSE);
  if(jp<0){
    jp = npts-1; 
    g->GetPoint(jp, x1, y1);
    ip = npts-1;
    boundary = kTRUE;
  }else if(jp==0){ 
    ip = 0;
    boundary = kTRUE;
  }else{ 
    ip = jp-1;
    g->GetPoint(ip, x0, y0);
    invdx = 1./(x0-x1);
  }
  // process pt dependences
  // +++ process dy
  Double_t dym = boundary?y1:((pt*(y0-y1) + (x0*y1-x1*y0))*invdx),
           sym = 0.5*(g->GetErrorY(ip)+g->GetErrorY(jp));
  dd      = (dyr - dym)/sym; dd*=dd;
  f[0] = TMath::Exp(-0.5*dd);
  // +++ process dphi
  if(!(g = (TGraphErrors*)fLike->At(offset+Int_t(kResAng)*kNcharge+Int_t(chg)))){
    AliError("Failed retrieving AttachClusters graph.");
    return 0.;
  }
  g->GetPoint(ip, x0, y0);g->GetPoint(jp, x1, y1);
  Double_t dpm = boundary?y1:((pt*(y0-y1) + (x0*y1-x1*y0))*invdx),
           spm = 0.5*(g->GetErrorY(ip)+g->GetErrorY(jp));
  dd      = (dphi - dpm)/spm; dd*=dd;
  f[1] = TMath::Exp(-0.5*dd);
  // +++ process no of clusters
  if(!(g = (TGraphErrors*)fLike->At(offset+Int_t(kNclMean)*kNcharge+Int_t(chg)))){
    AliError("Failed retrieving AttachClusters graph.");
    return 0.;
  }
  g->GetPoint(ip, x0, y0);g->GetPoint(jp, x1, y1);
  Double_t nm = boundary?y1:((pt*(y0-y1) + (x0*y1-x1*y0))*invdx);
  f[2] = (nm-TMath::Abs(n-nm))/nm;
 
  // process phi dependences
  // +++ process <s>/s
  if(!(g = (TGraphErrors*)fLike->At(offset+Int_t(kSigma)*kNcharge+Int_t(chg)))){
    AliError("Failed retrieving AttachClusters graph.");
    return 0.;
  }
  // Interpolate phi [deg]
  npts=g->GetN(); jp=-1;
  for(Int_t kp(0); kp<npts; kp++){
    g->GetPoint(kp, x1, y1);
    if(x1>=phiTrk){jp=kp; break;}
  }
  if(jp<0){
    jp = npts-1; 
    g->GetPoint(jp, x1, y1);
    ip = jp;
    boundary = kTRUE;
  }else if(jp==0){ 
    ip = jp;
    boundary = kTRUE;
  }else{ 
    ip = jp-1;
    g->GetPoint(ip, x0, y0);
    invdx = 1./(x0-x1);
    boundary = kFALSE;
  }
  Double_t sm = boundary?y1:((phiTrk*(y0-y1) + (x0*y1-x1*y0))*invdx),
           ssm = 0.5*(g->GetErrorY(ip)+g->GetErrorY(jp));
  dd      = (sr - sm)/ssm; dd*=dd;
  f[3] = 1.;//TMath::Exp(-0.5*dd);

  // Calculate likelihood
  Double_t length = f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3];
  length = TMath::Sqrt(length);
  Double_t cosTht = f[0]+f[1]+f[2]+f[3];
  cosTht /= (4.*length);
  AliDebug(2, Form("Like[%5.3f] ThtLike[%6.2f](deg)\n"
    "    F_dy (%+6.2f)=%4.2f\n" 
    "    F_phi(%+6.2f)=%4.2f\n"
    "    F_ncl(%+6d)=%4.2f\n"
    "    F_<s>(%+6.2f)=%4.2f", 
    length, TMath::ACos(cosTht)*TMath::RadToDeg(),
    dyr, f[0], dphi, f[1], n, f[2], sr, f[3]));

  return length;
}

//______________________________________________________________________________
Bool_t AliTRDCalTrkAttach::LoadReferences(const Char_t *file)
{
// Load calibration data from file

  if(!file || !TFile::Open(file)){
    AliError("Parametrization file missing or unreadable.");
    return kFALSE;
  }
  TGraphErrors *g(NULL);
  Char_t co[kNcalib] = {'y', 'p', 's', 'n'},
         cs[2] = {'n', 'p'};

  if(fLike) fLike->Clear();
  else fLike = new TObjArray(AliTRDgeometry::kNlayer*kNcharge*kNcalib);
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t icalib(0); icalib<kNcalib; icalib++){
      for(Int_t isgn(0); isgn<kNcharge; isgn++){
        if(!(g = (TGraphErrors*)gFile->Get(Form("%c%c%d", co[icalib], cs[isgn], ily)))) return kFALSE;
        fLike->AddAt(g, kNcharge*(ily*kNcalib+icalib)+isgn);
      }
    }
  }
  return kTRUE;
}
