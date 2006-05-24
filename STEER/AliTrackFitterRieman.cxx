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

///////////////////////////////////////////////////////////////////////////////
//
// Class to the track points on the Riemann sphere. Inputs are
// the set of id's (volids) of the volumes in which residuals are
// calculated to construct a chi2 function to be minimized during 
// the alignment procedures. For the moment the track extrapolation is 
// taken at the space-point reference plane. The reference plane is
// found using the covariance matrix of the point
// (assuming sigma(x)=0 at the reference coordinate system.
//
// Internal usage of AliRieman class for minimization
//
//////////////////////////////////////////////////////////////////////////////

#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TArrayI.h"
#include "AliTrackFitterRieman.h"
#include "AliLog.h"
#include "TTreeStream.h"
#include "AliRieman.h"
#include "AliLog.h"

ClassImp(AliTrackFitterRieman)


AliTrackFitterRieman::AliTrackFitterRieman():
  AliTrackFitter()
{
  //
  // default constructor
  //
  fAlpha = 0.;
  fNUsed = 0;
  fConv = kFALSE;
  fMaxDelta = 3;
  fDebugStream = new TTreeSRedirector("RiemanAlignDebug.root");
  fRieman = new AliRieman(10000);  // allocate rieman
}


AliTrackFitterRieman::AliTrackFitterRieman(AliTrackPointArray *array, Bool_t owner):
  AliTrackFitter(array,owner)
{
  //
  // Constructor
  //
  fAlpha = 0.;
  fNUsed = 0;
  fConv = kFALSE;
  fMaxDelta = 3;
  if (AliLog::GetDebugLevel("","AliTrackFitterRieman")) fDebugStream = new TTreeSRedirector("RiemanAlignDebug.root");
  fRieman = new AliRieman(10000);  //allocate rieman
}

AliTrackFitterRieman::AliTrackFitterRieman(const AliTrackFitterRieman &rieman):
  AliTrackFitter(rieman)
{
  //
  // copy constructor
  //
  fAlpha = rieman.fAlpha;
  fNUsed = rieman.fNUsed;
  fConv = rieman.fConv;
  fMaxDelta = rieman.fMaxDelta;
  fRieman = new AliRieman(*(rieman.fRieman));
}

//_____________________________________________________________________________
AliTrackFitterRieman &AliTrackFitterRieman::operator =(const AliTrackFitterRieman& rieman)
{
  //
  // Assignment operator
  //
  if(this==&rieman) return *this;
  ((AliTrackFitter *)this)->operator=(rieman);

  fAlpha  = rieman.fAlpha;
  fNUsed  = rieman.fNUsed;
  fConv   = rieman.fConv;
  fMaxDelta = rieman.fMaxDelta;
  fRieman = new AliRieman(*(rieman.fRieman));
  return *this;
}


AliTrackFitterRieman::~AliTrackFitterRieman(){
  //
  //
  //
  delete fRieman;
  delete fDebugStream;
}

void AliTrackFitterRieman::Reset()
{
  // Reset the track parameters and
  // rieman sums
  //
  AliTrackFitter::Reset();
  fRieman->Reset();
  fAlpha = 0.;
  fNUsed = 0;
  fConv =kFALSE;
}

Bool_t AliTrackFitterRieman::Fit(const TArrayI *volIds,const TArrayI *volIdsFit,
				 AliAlignObj::ELayerID layerRangeMin,
			 AliAlignObj::ELayerID layerRangeMax)
{
  // Fit the track points. The method takes as an input
  // the set of id's (volids) of the volumes in which
  // one wants to calculate the residuals.
  // The following parameters are used to define the
  // range of volumes to be used in the fitting
  // As a result two AliTrackPointArray's obects are filled.
  // The first one contains the space points with
  // volume id's from volids list. The second array of points represents
  // the track extrapolations corresponding to the space points
  // in the first array. The two arrays can be used to find
  // the residuals in the volids and consequently construct a
  // chi2 function to be minimized during the alignment
  // procedures. For the moment the track extrapolation is taken
  // at the space-point reference plane. The reference plane is
  // found using the covariance matrix of the point
  // (assuming sigma(x)=0 at the reference coordinate system.

  Reset();
  Int_t npoints = fPoints->GetNPoints();
  if (fPoints && AliLog::GetDebugLevel("","AliTrackFitterRieman")>1){
    Int_t nVol    = volIds->GetSize();
    Int_t nVolFit = volIdsFit->GetSize();
    Int_t volId  = volIds->At(0);    
    (*fDebugStream)<<"PInput"<<
      "NPoints="<<npoints<<    // number of points
      "VolId="<<volId<<        // first vol ID
      "NVol="<<nVol<<          // number of volumes 
      "NvolFit="<<nVolFit<<    // number of volumes to fit
      "fPoints.="<<fPoints<<   // input points
      "\n";
  }
  if (npoints < 3) return kFALSE;

  Bool_t isAlphaCalc = kFALSE;
  AliTrackPoint p,plocal;
//   fPoints->GetPoint(p,0);
//   fAlpha = TMath::ATan2(p.GetY(),p.GetX());

  Int_t npVolId = 0;
  fNUsed = 0;
  Int_t *pindex = new Int_t[npoints];
  for (Int_t ipoint = 0; ipoint < npoints; ipoint++)
    {
      fPoints->GetPoint(p,ipoint);
      UShort_t iVolId = p.GetVolumeID();
      if (FindVolId(volIds,iVolId)) {
	pindex[npVolId] = ipoint;
	npVolId++;
      }
      if (volIdsFit != 0x0) {
	if (!FindVolId(volIdsFit,iVolId)) continue;
      }
      else {
	if (iVolId < AliAlignObj::LayerToVolUID(layerRangeMin,0) ||
	    iVolId > AliAlignObj::LayerToVolUID(layerRangeMax,
						AliAlignObj::LayerSize(layerRangeMax))) continue;
      }
      if (!isAlphaCalc) {
	fAlpha = p.GetAngle();
	isAlphaCalc = kTRUE;
      }
      plocal = p.Rotate(fAlpha);
      if (TMath::Abs(plocal.GetX())>500 || TMath::Abs(plocal.GetX())<2 || plocal.GetCov()[3]<=0 ||plocal.GetCov()[5]<=0 ){
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<</n");
	p.Dump();
	plocal.Dump();
	printf("Problematic point\n");	
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<</n");
      }else{
	AddPoint(plocal.GetX(),plocal.GetY(),plocal.GetZ(),
		 TMath::Sqrt(plocal.GetCov()[3]),TMath::Sqrt(plocal.GetCov()[5]));
      }
      // fNUsed++;  AddPoint should be responsible
    }

  if (npVolId == 0 || fNUsed < fMinNPoints) {
    delete [] pindex;
    return kFALSE;
  }

  Update();
 
  if (!fConv) {
    delete [] pindex;
    return kFALSE;
  }

  if ((fParams[0] == 0) ||
      ((-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1) <= 0)) {
    delete [] pindex;
    return kFALSE;
  }


  if (fNUsed < fMinNPoints) {
    delete [] pindex;
    return kFALSE;
  }

  fPVolId = new AliTrackPointArray(npVolId);
  fPTrack = new AliTrackPointArray(npVolId);
  AliTrackPoint p2;
  for (Int_t ipoint = 0; ipoint < npVolId; ipoint++)
    {
      Int_t index = pindex[ipoint];
      fPoints->GetPoint(p,index);
      if (GetPCA(p,p2) && (
	  TMath::Abs(p.GetX()-p2.GetX())<fMaxDelta && 
	  TMath::Abs(p.GetY()-p2.GetY())<fMaxDelta &&
	  TMath::Abs(p.GetZ()-p2.GetZ())<fMaxDelta
	  )) {
	Float_t xyz[3],xyz2[3];
	p.GetXYZ(xyz); p2.GetXYZ(xyz2);
	//	printf("residuals %f %d %d %f %f %f %f %f %f\n",fChi2,fNUsed,fConv,xyz[0],xyz[1],xyz[2],xyz2[0]-xyz[0],xyz2[1]-xyz[1],xyz2[2]-xyz[2]);	
	fPVolId->AddPoint(ipoint,&p);
	fPTrack->AddPoint(ipoint,&p2);
      }else{
	// what should be default bahavior -  
	delete [] pindex;
	delete fPVolId;
	delete fPTrack;
	fPVolId =0;
	fPTrack =0;
	return kFALSE;
      }
    }  
  
  if (AliLog::GetDebugLevel("","AliTrackFitterRieman")>0){
    //
    // Debug Info
    //
    AliTrackPointArray *lPVolId = new AliTrackPointArray(npVolId);
    AliTrackPointArray *lPTrack = new AliTrackPointArray(npVolId);
    AliRieman * residual = fRieman->MakeResiduals();
    AliTrackPoint p2local;
    for (Int_t ipoint = 0; ipoint < npVolId; ipoint++){
      Int_t index = pindex[ipoint];
      fPoints->GetPoint(p,index); 
      if (GetPCA(p,p2)) {
	plocal = p.Rotate(fAlpha);
	//	plocal.Rotate(fAlpha);
	Float_t xyz[3],xyz2[3];
	p2local = plocal;
	plocal.GetXYZ(xyz); 
	xyz2[0] = xyz[0];
	xyz2[1] = GetYat(xyz[0]); 
	xyz2[2] = GetZat(xyz[0]); 
	p2local.SetXYZ(xyz2);
	lPVolId->AddPoint(ipoint,&plocal);
	lPTrack->AddPoint(ipoint,&p2local);
      }
    }      
    //
    // debug info
    //
    Int_t nVol    = volIds->GetSize();
    Int_t nVolFit = volIdsFit->GetSize();
    Int_t volId   = volIds->At(0);   
    Int_t modId   =0;
    Int_t layer   =  AliAlignObj::VolUIDToLayer(volId,modId);
    
    (*fDebugStream)<<"Fit"<<
      "VolId="<<volId<<        // volume ID
      "Layer="<<layer<<        // layer ID
      "Module="<<modId<<       // module ID
      "NVol="<<nVol<<          // number of volumes 
      "NvolFit="<<nVolFit<<    // number of volumes to fit
      "Points0.="<<fPVolId<<    // original points
      "Points1.="<<fPTrack<<    // fitted points 
      "LPoints0.="<<lPVolId<<   // original points - local frame
      "LPoints1.="<<lPTrack<<   // fitted points   - local frame      
      "Rieman.="<<this<<        // original rieman fit
      "Res.="<<residual<<       // residuals of rieman fit  
      "\n";
    delete lPVolId;
    delete lPTrack;
    delete residual;
  }

  delete [] pindex;
  return kTRUE;
}


void AliTrackFitterRieman::AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz)
{
  //
  // add point to rieman fitter
  //
  fRieman->AddPoint(x,y,z,sy,sz);
  fNUsed = fRieman->GetN();
}



void AliTrackFitterRieman::Update(){
  //
  // 
  //
  fRieman->Update();
  fConv = kFALSE;
  if (fRieman->IsValid()){
    for (Int_t ipar=0; ipar<6; ipar++){
      fParams[ipar] = fRieman->GetParam()[ipar];
    }
    fChi2  =  fRieman->GetChi2();
    fNdf   =  fRieman->GetN()- 2;
    fNUsed =  fRieman->GetN();
    fConv = kTRUE;
  }
}

//_____________________________________________________________________________
Bool_t AliTrackFitterRieman::GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const
{
  //
  // Get the closest to a given spacepoint track trajectory point
  // Look for details in the description of the Fit() method
  //
  if (!fConv) return kFALSE;

  // First X and Y coordinates
  Double_t sin = TMath::Sin(fAlpha);
  Double_t cos = TMath::Cos(fAlpha);
  //  fParam[0]  = 1/y0
  //  fParam[1]  = -x0/y0
  //  fParam[2]  = - (R^2 - x0^2 - y0^2)/y0
  if (fParams[0] == 0) return kFALSE;
  // Track parameters in the global coordinate system
  Double_t x0 = -fParams[1]/fParams[0]*cos -         1./fParams[0]*sin;
  Double_t y0 =          1./fParams[0]*cos - fParams[1]/fParams[0]*sin;
  if ((-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1) <= 0) return kFALSE;
  Double_t r  = TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1)/
                fParams[0];

  // Define space-point refence plane
  Double_t alphap = p.GetAngle();
  Double_t sinp = TMath::Sin(alphap);
  Double_t cosp = TMath::Cos(alphap);
  Double_t x  = p.GetX()*cosp + p.GetY()*sinp;
  Double_t y  = p.GetY()*cosp - p.GetX()*sinp;
  Double_t x0p= x0*cosp + y0*sinp;
  Double_t y0p= y0*cosp - x0*sinp;
  if ((r*r - (x-x0p)*(x-x0p))<0) {
    AliWarning(Form("Track extrapolation failed ! (Track radius = %f, track circle x = %f, space-point x = %f, reference plane angle = %f\n",r,x0p,x,alphap));
    return kFALSE;
  }
  Double_t temp = TMath::Sqrt(r*r - (x-x0p)*(x-x0p));
  Double_t y1 = y0p + temp;
  Double_t y2 = y0p - temp;
  Double_t yprime = y1;
  if(TMath::Abs(y2-y) < TMath::Abs(y1-y)) yprime = y2;

  // Back to the global coordinate system
  Double_t xsecond = x*cosp - yprime*sinp;
  Double_t ysecond = yprime*cosp + x*sinp;

  // Now Z coordinate and track angles
  Double_t x2 = xsecond*cos + ysecond*sin;
  Double_t zsecond = GetZat(x2);
  Double_t dydx = GetDYat(x2);
  Double_t dzdx = GetDZat(x2);

  // Fill the cov matrix of the track extrapolation point
  Double_t cov[6] = {0,0,0,0,0,0};
  Double_t sigmax = 100*100.;
  cov[0] = sigmax;           cov[1] = sigmax*dydx;      cov[2] = sigmax*dzdx;
  cov[3] = sigmax*dydx*dydx; cov[4] = sigmax*dydx*dzdx;
  cov[5] = sigmax*dzdx*dzdx;

  Float_t  newcov[6];
  newcov[0] = cov[0]*cos*cos-
            2*cov[1]*sin*cos+
              cov[3]*sin*sin;
  newcov[1] = cov[1]*(cos*cos-sin*sin)-
             (cov[3]-cov[0])*sin*cos;
  newcov[2] = cov[2]*cos-
              cov[4]*sin;
  newcov[3] = cov[0]*sin*sin+
            2*cov[1]*sin*cos+
              cov[3]*cos*cos;
  newcov[4] = cov[4]*cos+
              cov[2]*sin;
  newcov[5] = cov[5];

  p2.SetXYZ(xsecond,ysecond,zsecond,newcov);
  if (AliLog::GetDebugLevel("","AliTrackFitterRieman")>1){
    AliTrackPoint lp0(p);
    AliTrackPoint lp2(p2);
    if (0)(*fDebugStream)<<"PCA"<<
      "P0.="<<&lp0<<
      "P2.="<<&lp2<<
      "\n";
  }
  return kTRUE;
}
