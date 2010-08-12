//========================================================================
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
//                                                                        
// Author: The ALICE Off-line Project.                                    
// Contributors are mentioned in the code where appropriate.              
//                                                                        
// Permission to use, copy, modify and distribute this software and its   
// documentation strictly for non-commercial purposes is hereby granted   
// without fee, provided that the above copyright notice appears in all   
// copies and that both the copyright notice and this permission notice   
// appear in the supporting documentation. The authors make no claims     
// about the suitability of this software for any purpose. It is          
// provided "as is" without express or implied warranty.                  
//========================================================================  
//                        
//                       Class AliEMCALTrack 
//                      ---------------------
//    A class implementing a track which is propagated to EMCAL and 
//    matches and EMCAL cluster. 
//    This track object will not update Kalman parameters, but it 
//    allows for track propagation and suitable energy loss correction,
//    even in an environment with a variable magnetic field, which is not
//    well managed in the AliExternalTrackParam class.
// ------------------------------------------------------------------------
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//=========================================================================

#include "Riostream.h"

#include "TVector3.h"

#include "AliLog.h"
#include "AliESDtrack.h" 

#include "AliEMCALTrack.h"

Bool_t AliEMCALTrack::fgUseOuterParams = kTRUE;
Bool_t AliEMCALTrack::fgCorrectForEL   = kFALSE;
Bool_t AliEMCALTrack::fgSortByPt       = kTRUE;

ClassImp(AliEMCALTrack)
//
//------------------------------------------------------------------------------
//
AliEMCALTrack::AliEMCALTrack() 
  : AliExternalTrackParam(),
    fClusterIndex(-1),
    fClusterDist(1000.0),          // default: extremely large distance
    fMass(0.13957018),             // default: pion mass
    fSeedIndex(-1),
    fSeedLabel(-1)
{
	//
	// Default constructor.
	// Sets to meaningless values the indexes corresponding to
	// ESD seed track and matched cluster.
	//

}
//
//------------------------------------------------------------------------------
//
AliEMCALTrack::AliEMCALTrack(const AliESDtrack& t) 
  : AliExternalTrackParam(),
    fClusterIndex(-1),
    fClusterDist(1000.0),
    fMass(t.GetMass()),
    fSeedIndex(-1),
    fSeedLabel(t.GetLabel())
{
	//
	// Constructor from AliESDtrack
	//

	// parameters are chosen according to static variable fUseOuterParams
	Double_t alpha=0., x=0., params[5], cov[15];
	if (fgUseOuterParams) {
	  if(t.GetOuterParam()){
	    t.GetOuterExternalParameters(alpha, x, params);
	    t.GetOuterExternalCovariance(cov);
	  }
	  else{ // no outer param available leave the default as is
	    return;
	  }
	}
	else {
		alpha = t.GetAlpha();
		t.GetExternalParameters(x, params);
		t.GetExternalCovariance(cov);
	}
	
	if (alpha < -TMath::Pi()) alpha += TMath::TwoPi();
	else if (alpha >= TMath::Pi()) alpha -= TMath::TwoPi();
	
	// set this accordingly
	Set(x, alpha, params, cov);
}
//
//------------------------------------------------------------------------------
//
AliEMCALTrack::AliEMCALTrack(const AliEMCALTrack& t) 
  : AliExternalTrackParam(t),
    fClusterIndex(t.fClusterIndex),
    fClusterDist(t.fClusterDist),
    fMass(t.fMass),
    fSeedIndex(t.fSeedIndex),
    fSeedLabel(t.fSeedLabel)
{
	//
	// Copy constructor.
	//
}
//
//------------------------------------------------------------------------------
//
AliEMCALTrack& AliEMCALTrack::operator=(const AliEMCALTrack &t)
{
	//
	// Assignment operator
	// Works like copy constructor
	//
	
	fClusterIndex = t.fClusterIndex;
	fClusterDist = t.fClusterDist;
	
	fMass = t.fMass;
	
	fSeedIndex = t.fSeedIndex;
	fSeedLabel = t.fSeedLabel;
	return *this;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTrack::Compare(const TObject *obj) const 
{
	//
	// Compare tracks.
	// How tracks are compared depends on the static flag
	// "fSortByPt" (boolean):
	// true  => tracks are compared w.r. to their transverse momentum
	// false => tracks are compared w.r. to their distance from cluster
	//
	
	AliEMCALTrack *that = (AliEMCALTrack*)obj;
	
	Double_t thisP[3], thisVal=0., thatP[3], thatVal=0.;
	
	if (fgSortByPt) {
		this->GetPxPyPz(thisP);
		that->GetPxPyPz(thatP);	
		thisVal = TMath::Sqrt(thisP[0]*thisP[0] + thisP[1]*thisP[1]);
		thatVal = TMath::Sqrt(thatP[0]*thatP[0] + thatP[1]*thatP[1]);
	}
	else {
		thisVal = this->GetClusterDist();
		thatVal = that->GetClusterDist();
	}
	
	if (thisVal > thatVal) return 1;
	else if (thisVal < thatVal) return -1;
	else return 0;
}
//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTrack::PropagateTo(Double_t xk, Double_t d, Double_t x0)
{
	//
	// Propagates a track to the plane defined by x='xk'.
	// Second parameter is the width (in units of rad. length) crossed by the track.
	// Third parameter is the reference radiation length used.
	// Track propagation includes computing energy loss (modifies curvature)
	// and multiple scattering perturbation (alters covariance matrix), if requested.
	// Method returns kFALSE when something goes wrong with computations.
	//
	// An additional operation (thanks to Yuri Belikov) is done to check
	// when the track crosses a sector boundary. If this happens, 
	// the local track reference frame is adjusted accordingly.
	//
		
	Double_t y=0.;
	Double_t field = GetBz();
	Double_t width = TMath::Pi() / 9.0; // width of TPC/TRD/EMCAL sector (= 20 deg)
	Double_t ymax  = TMath::Abs(xk * TMath::Tan(0.5 * width)); // max allowed Y in local coords at distance xk
	
	// first check: try to compute the local 'y' at the distance 'xk':
	// if this attempt fails, the propagation cannot be done
	if (!GetYAt(xk, field, y)) return kFALSE;
	
	// if is -ymax < y < ymax ==> 'direct' propagation is done;
	if (TMath::Abs(y) <= ymax) return SimplePropagation(xk, d, x0);
	
	// otherwise, try change a sector to find one where the propagation is ok
	Int_t    i=0, incr=0, istart=0, nloops=0;
	Double_t alpha = GetAlpha();
	incr = (y > ymax) ? 1 : -1;
	if (alpha < 0.0) alpha += TMath::TwoPi();
	istart = (Int_t)(alpha / width);
	for (i = istart, nloops = 0; nloops < 18; i += incr, nloops++) {
		if (i == 18) i = 0;
		if (i == -1) i = 17;
		alpha = ((Double_t)i + 0.5) * width;
		if (Rotate(alpha)) {
			if (GetYAt(xk, field, y)) {
				if (TMath::Abs(y) <= ymax) {
				  AliDebug(1,Form("Required change from sector %d to sector %d to succeed in propagation", istart, i));
					return SimplePropagation(xk, d, x0);
				}
			}
		}
	}
	
	// if the routine exits from the loop and reaches this point,
	// it means that none of the rotations succeeded
	AliWarning("Track propagation fails in every sector. Impossible to propagate.");
	return kFALSE;
}
//
//------------------------------------------------------------------------------
//
Double_t AliEMCALTrack::StraightPropagateTo(Double_t xk, Double_t &x, Double_t &y, Double_t &z)
{
	//
	// Does propagation with a straight line approximation.
	// This operation does not update track parameters, but it returns a point
	// in space, according to this propagation, which is stored in
	// the arguments #2, #3 and #4
	//
	
	Double_t oldX = GetX(), oldY = GetY(), oldZ = GetZ();
	Double_t newPos[3];
	
	newPos[0] = xk;
	newPos[1] = oldY * xk / oldX;
	newPos[2] = oldZ * xk / oldX;
	
	Double_t cs = TMath::Cos(GetAlpha()), sn = TMath::Sin(GetAlpha());
	x = newPos[0]*cs - newPos[1]*sn;
	y = newPos[0]*sn + newPos[1]*cs;
	z = newPos[2];
	
	return newPos[1];
}
//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTrack::PropagateToGlobal(Double_t x, Double_t y, Double_t z, Double_t d, Double_t x0)
{
	//
	// Propagate to a point specified with its global XYZ coordinates.
	// Here, the correct value of the 'X' parameter to be sent to "PropagateTo" is computed.
	//
	
	TVector3 vc(x, y, z);
	Double_t width = 20.0; // width of TPC/TRD/EMCAL sector
	Double_t phi   = vc.Phi() * TMath::RadToDeg();
	if (phi < 0.0) phi += 360.0;
	
	Int_t    isector = (Int_t)(phi / width);
	Double_t rotation = ((Double_t)isector + 0.5) * width;
	vc.RotateZ(-rotation * TMath::DegToRad());
	
	return PropagateTo(vc.X(), d, x0);
}
//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTrack::SimplePropagation(Double_t xk, Double_t d, Double_t x0)
{
  //
  // Recall base class method for track propagation.
  //
  
  Double_t field[3];

  GetBxByBz(field);
	
  // propagation...
  if (!AliExternalTrackParam::PropagateToBxByBz(xk, field)) return kFALSE;
  
  // EL correction is computed only if requested...
  if (!fgCorrectForEL) return kTRUE;
  return AliExternalTrackParam::CorrectForMeanMaterial(d, x0, GetMass());
}
