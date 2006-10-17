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

#include "AliLog.h"
#include "AliTracker.h"
#include "AliESDtrack.h" 

#include "AliEMCALTrack.h"

Bool_t AliEMCALTrack::fgUseOuterParams = kTRUE;

ClassImp(AliEMCALTrack)
//
//------------------------------------------------------------------------------
//
AliEMCALTrack::AliEMCALTrack() 
  : AliExternalTrackParam(),
    fClusterIndex(-1),
    fClusterDist(1000.0),
    fMass(0.0),
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
	Double_t alpha, x, params[5], cov[15];
	if (fgUseOuterParams) {
		t.GetOuterExternalParameters(alpha, x, params);
		t.GetOuterExternalCovariance(cov);
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
	// Tracks compared wrt their distance from matched point
	//
	
	AliEMCALTrack *that = (AliEMCALTrack*)obj;
	
	Double_t thisDist = GetClusterDist();
	Double_t thatDist = that->GetClusterDist();
	
	if (thisDist > thatDist) return 1;
	else if (thisDist < thatDist) return -1;
	return 0;
}
//
//------------------------------------------------------------------------------
//
Double_t AliEMCALTrack::GetBz() const 
{
	//
	// returns Bz component of the magnetic field in kG,
	// at the point where track is actually propagated
	//
	
	Double_t r[3];
	
	if (AliTracker::UniformField()) {
		return AliTracker::GetBz();
	}
	else {
		GetXYZ(r);
		return AliTracker::GetBz(r);
	}
}
//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTrack::PropagateTo(Double_t xk, Double_t eqWidth, Double_t x0)
{
	//
	// Propagates a track to the plane defined by x='xk'.
	// Second parameter represents the percent width of radiation length
	// crossed by the track.
	// Third parameter is the reference radiation length used (default: Si).
	// Track propagation includes computing energy loss (modifies curvature)
	// and multiple scattering perturbation (alters covariance matrix).
	// Method returns kFALSE when something goes wrong with computations.
	//
	
	Double_t field = GetBz();
	
	if (!AliExternalTrackParam::PropagateTo(xk, field)) return kFALSE;
	if (!AliExternalTrackParam::CorrectForMaterial(eqWidth, x0, GetMass())) return kFALSE;

	return kTRUE;            
}     
