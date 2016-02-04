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

#include "TMath.h"

#include "AliLog.h"

#include "AliMFTTrackParam.h"


/// \cond CLASSIMP
ClassImp(AliMFTTrackParam); // Class implementation in ROOT context
														/// \endcond


//=============================================================================================

AliMFTTrackParam::AliMFTTrackParam():TObject(),
fCovariances(NULL),
fPropagator(NULL),
fExtrapParameters(NULL),
fExtrapCovariances(NULL),
fSmoothParameters(NULL),
fSmoothCovariances(NULL),
fParameters(5,1),
fX(0.),
fY(0.),
fZ(0.),
fTrackChi2(0.),
fLocalChi2(0.)
{
	/// Default constructor
	fParameters.Zero();
	
}

//=============================================================================================
AliMFTTrackParam::AliMFTTrackParam(const AliMFTTrackParam& theMFTTrackParam)
: TObject(theMFTTrackParam),
fCovariances(NULL),
fPropagator(NULL),
fExtrapParameters(NULL),
fExtrapCovariances(NULL),
fSmoothParameters(NULL),
fSmoothCovariances(NULL),
fParameters(theMFTTrackParam.fParameters),
fX(theMFTTrackParam.fX),
fY(theMFTTrackParam.fY),
fZ(theMFTTrackParam.fZ),
fTrackChi2(theMFTTrackParam.fTrackChi2),
fLocalChi2(theMFTTrackParam.fLocalChi2)
{
	/// Copy constructor
	if (theMFTTrackParam.fCovariances) fCovariances = new TMatrixD(*(theMFTTrackParam.fCovariances));
	if (theMFTTrackParam.fPropagator) fPropagator = new TMatrixD(*(theMFTTrackParam.fPropagator));
	if (theMFTTrackParam.fExtrapParameters) fExtrapParameters = new TMatrixD(*(theMFTTrackParam.fExtrapParameters));
	if (theMFTTrackParam.fExtrapCovariances) fExtrapCovariances = new TMatrixD(*(theMFTTrackParam.fExtrapCovariances));
	if (theMFTTrackParam.fSmoothParameters) fSmoothParameters = new TMatrixD(*(theMFTTrackParam.fSmoothParameters));
	if (theMFTTrackParam.fSmoothCovariances) fSmoothCovariances = new TMatrixD(*(theMFTTrackParam.fSmoothCovariances));
	
	//	if(fOwnCluster) fClusterPtr = static_cast<AliMUONVCluster*>(theMFTTrackParam.fClusterPtr->Clone());
	//	else fClusterPtr = theMFTTrackParam.fClusterPtr;
}


//=============================================================================================


AliMFTTrackParam::~AliMFTTrackParam() {
	DeleteCovariances();
	
}

//__________________________________________________________________________
Double_t AliMFTTrackParam::P()  const
{
	/// return total momentum
	Double_t invPt = GetInverseTransverseMomentum();
	if(TMath::Abs(invPt)<1e-6) return 0.;
	return TMath::Sqrt(1+ 1./(GetSlopeX()*GetSlopeX() + GetSlopeY()*GetSlopeY()))/invPt;
	
}
//__________________________________________________________________________
const TMatrixD& AliMFTTrackParam::GetCovariances() const
{
	/// Return the covariance matrix (create it before if needed)
	if (!fCovariances) {
		
		fCovariances = new TMatrixD(5,5);
		fCovariances->Zero();
	}
	return *fCovariances;
}

//__________________________________________________________________________
void AliMFTTrackParam::SetCovariances(const TMatrixD& covariances)
{
	/// Set the covariance matrix
	if (fCovariances) *fCovariances = covariances;
	else fCovariances = new TMatrixD(covariances);
}

//__________________________________________________________________________
void AliMFTTrackParam::SetCovariances(const Double_t matrix[5][5])
{
	/// Set the covariance matrix
	if (fCovariances) fCovariances->SetMatrixArray(&(matrix[0][0]));
	else fCovariances = new TMatrixD(5,5,&(matrix[0][0]));
}

//__________________________________________________________________________
void AliMFTTrackParam::SetVariances(const Double_t matrix[5][5])
{
	/// Set the diagonal terms of the covariance matrix (variances)
	if (!fCovariances) fCovariances = new TMatrixD(5,5);
	fCovariances->Zero();
	for (Int_t i=0; i<5; i++) (*fCovariances)(i,i) = matrix[i][i];
}

//__________________________________________________________________________
void AliMFTTrackParam::DeleteCovariances()
{
	/// Delete the covariance matrix
	delete fCovariances;
	fCovariances = 0x0;
}
//__________________________________________________________________________
const TMatrixD& AliMFTTrackParam::GetPropagator() const
{
	/// Return the propagator (create it before if needed)
	if (!fPropagator) {
		fPropagator = new TMatrixD(5,5);
		fPropagator->UnitMatrix();
	}
	return *fPropagator;
}

//__________________________________________________________________________
void AliMFTTrackParam::ResetPropagator()
{
	/// Reset the propagator
	if (fPropagator) fPropagator->UnitMatrix();
}

//__________________________________________________________________________
void AliMFTTrackParam::UpdatePropagator(const TMatrixD& propagator)
{
	/// Update the propagator
	if (fPropagator) *fPropagator = TMatrixD(propagator,TMatrixD::kMult,*fPropagator);
	else fPropagator = new TMatrixD(propagator);
}
//__________________________________________________________________________
void AliMFTTrackParam::Print(Option_t* opt) const
{
	/// Printing TrackParam information
	/// "full" option for printing all the information about the TrackParam
	TString sopt(opt);
	sopt.ToUpper();
 
	if ( sopt.Contains("FULL") ) {
		if(TMath::Abs(fParameters(4,0))>1.e-6) {
			AliInfo(Form("\t Pt       [GeV/c] = %.2e", 1./fParameters(4,0)));
			AliInfo(Form("\t P        [GeV/c] = %.2e", 1./fParameters(4,0)*TMath::Sqrt(1.+1./(fParameters(2,0)*fParameters(2,0) + fParameters(3,0)*fParameters(3,0)))));
		}
		else AliInfo("\t Pt       [GeV/c] = Infinite");
		AliInfo(Form("\t Slope X   = %+.4f  | Theta = %f", fParameters(2,0),GetTheta()*TMath::RadToDeg()));
		AliInfo(Form("\t Slope Y   = %+.4f  | Phi   = %f", fParameters(3,0),GetPhi()*TMath::RadToDeg()));
		AliInfo(Form("\t (X,Y,Z) [cm]    = (%e, %e, %e)", fParameters(0,0), fParameters(1,0), fZ));
		AliInfo(Form("\t Chi2 Local = %e Chi2 Track : %e", fLocalChi2,fTrackChi2));
	} else {
		
	}
	
}
