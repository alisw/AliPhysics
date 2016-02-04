#ifndef AliMFTTrackParam_H
#define AliMFTTrackParam_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup MFTrec
/// \class AliMFTTrackParam
/// \brief Class holding the parameter of a MFT Standalone Track
///
///
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, IPN-Lyon
/// \date April 28th, 2015

#include "TObject.h"
#include <TMatrixD.h>
#include <TMath.h>

//=============================================================================================

class AliMFTTrackParam : public TObject {
	
public:
	
	AliMFTTrackParam();
	virtual ~AliMFTTrackParam();
	AliMFTTrackParam(const AliMFTTrackParam& theMFTTrackParam);
	
	/// set Cluster coordinates (cm)
	void     SetClusterPos(Double_t x, Double_t y, Double_t z) {fX = x;fY = y;fZ = z;}
	/// return cluster X coordinate (cm)
	Double_t GetClusterX() const {return fX;}
	/// return cluster Y coordinate (cm)
	Double_t GetClusterY() const {return fY;}
	
	/// return Z coordinate (cm)
	Double_t GetZ() const {return fZ;}
	/// set Z coordinate (cm)
	void     SetZ(Double_t z) {fZ = z;}
	/// return X coordinate (cm)
	Double_t GetX() const {return fParameters(0,0);}
	/// set X coordinate (cm)
	void     SetX(Double_t val) {fParameters(0,0) = val;}
	/// return Y coordinate (cm)
	Double_t GetY() const {return fParameters(1,0);}
	/// set Y coordinate (cm)
	void     SetY(Double_t val) {fParameters(1,0) = val;}
	/// return X slope
	Double_t GetSlopeX() const {return fParameters(2,0);}
	/// set X slope
	void     SetSlopeX(Double_t val) {fParameters(2,0) = val;}
	/// return  Y slope
	Double_t GetSlopeY() const {return fParameters(3,0);}
	/// set  Y slope
	void     SetSlopeY(Double_t val) {fParameters(3,0) = val;}
	/// return Inverse Momentum
	Double_t GetInverseTransverseMomentum() const {return fParameters(4,0);}
	/// set Inverse Momentum
	void     SetInverseTransverseMomentum(Double_t val) {fParameters(4,0) = val;}
	/// return the charge (assumed forward motion)
	Double_t GetCharge() const {return TMath::Sign(1.,fParameters(4,0));}
	/// set the charge (assumed forward motion)
	void     SetCharge(Double_t charge) {if (charge*fParameters(4,0) < 0.) fParameters(4,0) *= -1.;}
	/// return  Polar angle theta
	Double_t GetTheta() const {return TMath::ATan(TMath::Sqrt(fParameters(2,0)*fParameters(2,0) +  fParameters(3,0)*fParameters(3,0)));}
	/// return  Azimuthal angle phi
	Double_t GetPhi() const {return TMath::ATan2(fParameters(3,0),fParameters(2,0)) ;}
	
	
	/// return track parameters
	const TMatrixD& GetParameters() const {return fParameters;}
	/// set track parameters
	void            SetParameters(const TMatrixD& parameters) {fParameters = parameters;}
	/// add track parameters
	void            AddParameters(const TMatrixD& parameters) {fParameters += parameters;}
	/// return kTRUE if the covariance matrix exist, kFALSE if not
	Bool_t    CovariancesExist() const {return (fCovariances) ? kTRUE : kFALSE;}
	
	/// return the chi2 of the track when the associated cluster was attached
	Double_t GetTrackChi2() const {return fTrackChi2;}
	/// set the chi2 of the track when the associated cluster was attached
	void     SetTrackChi2(Double_t chi2) {fTrackChi2 = chi2;}
	/// return the local chi2 of the associated cluster with respect to the track
	Double_t GetLocalChi2() const {return fLocalChi2;}
	/// set the local chi2 of the associated cluster with respect to the track
	void     SetLocalChi2(Double_t chi2) {fLocalChi2 = chi2;}
	
	const TMatrixD& GetCovariances() const;
	void	SetCovariances(const TMatrixD& covariances);
	void	SetCovariances(const Double_t matrix[5][5]);
	void	SetVariances(const Double_t matrix[5][5]);
	void	DeleteCovariances();
	
	
	const TMatrixD& GetPropagator() const;
	void            ResetPropagator();
	void            UpdatePropagator(const TMatrixD& propagator);
	
	virtual void Print(Option_t* opt="") const;
	
	/// return total momentum
	Double_t P()  const;
	
	
protected:
	mutable TMatrixD *fCovariances; ///< \brief Covariance matrix of track parameters
	mutable TMatrixD *fPropagator;  ///< Jacobian used to extrapolate the track parameters and covariances to the actual z position
	mutable TMatrixD *fExtrapParameters;  //!<! Track parameters extrapolated to the actual z position (not filtered by Kalman)
	mutable TMatrixD *fExtrapCovariances; //!<! Covariance matrix extrapolated to the actual z position (not filtered by Kalman)
	
	mutable TMatrixD *fSmoothParameters;  //!<! Track parameters obtained using smoother
	mutable TMatrixD *fSmoothCovariances; //!<! Covariance matrix obtained using smoother
	
	/// Track parameters ordered as follow:      <pre>
	/// param0:   local X-coordinate of a track (cm)
	/// param1:   local Y-coordinate of a track (cm)
	/// param2:   Slope X/Z
	/// param3:   Slope Y/Z
	/// param4:   Q/p_t (1/(GeV/c)) </pre>
	TMatrixD fParameters; ///<  Track parameters
	
	Double_t fX; ///< Cluster X coordinate (cm)
	Double_t fY; ///< Cluster Y coordinate (cm)
	Double_t fZ; ///< Cluster Z coordinate (cm)
	
	Double_t fTrackChi2; ///< Chi2 of the track when the associated cluster was attached
	Double_t fLocalChi2; ///< Local chi2 of the associated cluster with respect to the track
	
	/// \cond CLASSIMP
	ClassDef(AliMFTTrackParam,1);
	/// \endcond
};

//=============================================================================================

#endif
