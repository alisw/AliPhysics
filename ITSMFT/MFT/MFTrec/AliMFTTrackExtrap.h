#ifndef AliMFTTrackExtrap_H
#define AliMFTTrackExtrap_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMFTTrackExtrap
/// \brief Track parameters in ALICE dimuon spectrometer
///
//////////////////////////////////////////////////////////////
/// Tools for track extrapolation in ALICE dimuon spectrometer
//////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMatrixD.h>

class AliMagF;
class AliMFTTrackParam;

class AliMFTTrackExtrap : public TObject 
{
 public:
	/// Constructor
  AliMFTTrackExtrap() : TObject(){};
	/// Destructor
  virtual ~AliMFTTrackExtrap(){};
  
  static void SetField();
  
  /// return kTRUE if the field is switched ON
  static Bool_t IsFieldON() {return fgFieldON;}
  
  static Double_t GetImpactParamFromBendingMomentum(Double_t bendingMomentum);
  static Double_t GetBendingMomentumFromImpactParam(Double_t impactParam);
  
  // Linearly extrapolate track parameters
  static void LinearExtrapToZ(AliMFTTrackParam* trackParam, Double_t zEnd);
  
  // Linearly extrapolate track parameters and covariances
  static void LinearExtrapToZCov(AliMFTTrackParam* trackParam, Double_t zEnd, Bool_t updatePropagator = kFALSE);
  
  // Extrapolate track parameters in magnetic field
  static Bool_t ExtrapToZ(AliMFTTrackParam *trackParam, Double_t zEnd);
  
  // Extrapolate track parameters and covariances in magnetic field
  static Bool_t ExtrapToZCov(AliMFTTrackParam* trackParam, Double_t zEnd, Bool_t updatePropagator = kFALSE);
  
  // Extrapolate track parameters to vertex, corrected for multiple scattering and energy loss effects
  // Add branson correction resolution and energy loss fluctuation to parameter covariances
  static void ExtrapToVertex(AliMFTTrackParam* trackParam,
                             Double_t xVtx, Double_t yVtx, Double_t zVtx,
                             Double_t errXVtx, Double_t errYVtx);
  
  // Extrapolate track parameters to vertex, corrected for multiple scattering effects only
  // Add branson correction resolution to parameter covariances
  static void ExtrapToVertexWithoutELoss(AliMFTTrackParam* trackParam,
					 Double_t xVtx, Double_t yVtx, Double_t zVtx,
					 Double_t errXVtx, Double_t errYVtx);
  
  // Extrapolate track parameters to vertex, corrected for energy loss effects only
  // Add dispersion due to multiple scattering and energy loss fluctuation to parameter covariances
  static void ExtrapToVertexWithoutBranson(AliMFTTrackParam* trackParam, Double_t zVtx);
  
  // Extrapolate track parameters to vertex without multiple scattering and energy loss corrections
  // Add dispersion due to multiple scattering to parameter covariances
  static void ExtrapToVertexUncorrected(AliMFTTrackParam* trackParam, Double_t zVtx);
  
  static Double_t TotalMomentumEnergyLoss(AliMFTTrackParam* trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx);
  
  static Double_t GetMCSAngle2(const AliMFTTrackParam& param, Double_t dZ, Double_t x0);
  static void     AddMCSEffect(AliMFTTrackParam *param, Double_t dZ, Double_t x0);
  
  static Bool_t ExtrapOneStepRungekutta(Double_t charge, Double_t step, const Double_t* vect, Double_t* vout);
  
	static Double_t LinearRegression(Int_t nVal, Double_t *xVal, Double_t *yVal, Double_t &p0, Double_t &p1);
	static Double_t Sagitta(Int_t nVal, Double_t *xVal, Double_t *yVal, Double_t &distL, Double_t &q2);
	static Double_t QuadraticRegression(Int_t nVal, Double_t *xVal, Double_t *yVal, Double_t &p0, Double_t &p1, Double_t &p2);
	static Double_t CircleRegression(Int_t nVal, Double_t *xVal, Double_t *yVal);

 private:
  static const Double_t fgkSimpleBPosition;     //!<! position of the dipole
  static const Double_t fgkSimpleBLength;       //!<! length of the dipole
  static       Double_t fgSimpleBValue;         //!<! magnetic field value at the centre
  static       Bool_t   fgFieldON;              //!<! kTRUE if the field is switched ON
  static const Bool_t   fgkUseHelix;		//!<! Tell whether to use Helix or not (default is Runge-Kutta)
  static const Int_t    fgkMaxStepNumber;	//!<! Maximum number of steps for track extrapolation
  static const Double_t fgkHelixStepLength;	//!<! Step lenght for track extrapolation (used in Helix)
  static const Double_t fgkRungeKuttaMaxResidue;//!<! Maximal distance (in Z) to destination to stop the track extrapolation (used in Runge-Kutta)
  
  // Functions

  /// Not implemented
  AliMFTTrackExtrap(const AliMFTTrackExtrap& trackExtrap);
  /// Not implemented
  AliMFTTrackExtrap& operator=(const AliMFTTrackExtrap& trackExtrap);

  static Bool_t ExtrapToZHelix(AliMFTTrackParam *trackParam, Double_t Z);
  static Bool_t ExtrapToZRungekutta(AliMFTTrackParam *trackParam, Double_t Z);
  
  static void ConvertTrackParamForExtrap(AliMFTTrackParam* trackParam, Double_t forwardBackward, Double_t *v3);
  static void RecoverTrackParam(Double_t *v3, Double_t Charge, AliMFTTrackParam* trackParam);
  
  static void ExtrapToVertex(AliMFTTrackParam* trackParam,
                             Double_t xVtx, Double_t yVtx, Double_t zVtx,
                             Double_t errXVtx, Double_t errYVtx,
                             Bool_t correctForMCS, Bool_t correctForEnergyLoss);
  
  static void AddMCSEffectInAbsorber(AliMFTTrackParam* trackParam, Double_t signedPathLength, Double_t f0, Double_t f1, Double_t f2);
  static void CorrectMCSEffectInAbsorber(AliMFTTrackParam* param,
                                         Double_t xVtx, Double_t yVtx, Double_t zVtx,
                                         Double_t errXVtx, Double_t errYVtx,
                                         Double_t absZBeg, Double_t pathLength, Double_t f0, Double_t f1, Double_t f2);
  static void CorrectELossEffectInAbsorber(AliMFTTrackParam* param, Double_t eLoss, Double_t sigmaELoss2);
  static Bool_t GetAbsorberCorrectionParam(Double_t trackXYZIn[3], Double_t trackXYZOut[3], Double_t pTotal,
                                           Double_t &pathLength, Double_t &f0, Double_t &f1, Double_t &f2,
                                           Double_t &meanRho, Double_t &totalELoss, Double_t &sigmaELoss2);
  
  static Double_t BetheBloch(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicZ, Double_t atomicZoverA);
  static Double_t EnergyLossFluctuation(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicZoverA);
  
  static void Cov2CovP(const TMatrixD &param, TMatrixD &cov);
  static void CovP2Cov(const TMatrixD &param, TMatrixD &cov);
  
  static void ExtrapOneStepHelix(Double_t charge, Double_t step, const Double_t *vect, Double_t *vout);
  static void ExtrapOneStepHelix3(Double_t field, Double_t step, const Double_t *vect, Double_t *vout);
  
  ClassDef(AliMFTTrackExtrap, 0) // Tools for track extrapolation in ALICE dimuon spectrometer
};
	
#endif
