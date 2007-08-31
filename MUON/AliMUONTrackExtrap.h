#ifndef ALIMUONTRACKEXTRAP_H
#define ALIMUONTRACKEXTRAP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrackExtrap
/// \brief Track parameters in ALICE dimuon spectrometer
///
//////////////////////////////////////////////////////////////
/// Tools for track extrapolation in ALICE dimuon spectrometer
//////////////////////////////////////////////////////////////

#include <TObject.h>

class AliMagF;
class AliMUONTrackParam;

class AliMUONTrackExtrap : public TObject 
{
 public:
	/// Constructor
  AliMUONTrackExtrap() : TObject(){};
	/// Destructor
  virtual ~AliMUONTrackExtrap(){};
  
	/// set field map
  static void SetField(const AliMagF* magField) {fgkField = magField;}
  
  static Double_t GetImpactParamFromBendingMomentum(Double_t bendingMomentum);
  static Double_t GetBendingMomentumFromImpactParam(Double_t impactParam);
  
  static void LinearExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd);
  static void ExtrapToZ(AliMUONTrackParam *trackParam, Double_t zEnd);
  static void ExtrapToZCov(AliMUONTrackParam* trackParam, Double_t zEnd, Bool_t updatePropagator = kFALSE);
  static void ExtrapToStation(AliMUONTrackParam *trackParamIn, Int_t station, AliMUONTrackParam *trackParamOut);
  static void ExtrapToVertexUncorrected(AliMUONTrackParam* trackParam, Double_t zVtx);
  static void ExtrapToVertex(AliMUONTrackParam *trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx,
  			     Bool_t CorrectForMCS = kTRUE, Bool_t CorrectForEnergyLoss = kTRUE);
  
  static Double_t TotalMomentumEnergyLoss(AliMUONTrackParam* trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx);
  
  static Double_t GetMCSAngle2(const AliMUONTrackParam& param, Double_t dZ, Double_t x0);
  static void     AddMCSEffect(AliMUONTrackParam *param, Double_t dZ, Double_t x0);
  
  static void ExtrapOneStepRungekutta(Double_t charge, Double_t step, Double_t* vect, Double_t* vout);
  
  
 private:
  static const AliMagF* fgkField;		//!< field map
  static const Bool_t   fgkUseHelix;		//!< Tell whether to use Helix or not (default is Runge-Kutta)
  static const Int_t    fgkMaxStepNumber;	//!< Maximum number of steps for track extrapolation
  static const Double_t fgkHelixStepLength;	//!< Step lenght for track extrapolation (used in Helix)
  static const Double_t fgkRungeKuttaMaxResidue;//!< Maximal distance (in Z) to destination to stop the track extrapolation (used in Runge-Kutta)
  
  // Functions

  /// Not implemented
  AliMUONTrackExtrap(const AliMUONTrackExtrap& trackExtrap);
  /// Not implemented
  AliMUONTrackExtrap& operator=(const AliMUONTrackExtrap& trackExtrap);

  static void ExtrapToZHelix(AliMUONTrackParam *trackParam, Double_t Z);
  static void ExtrapToZRungekutta(AliMUONTrackParam *trackParam, Double_t Z);
  
  static void ConvertTrackParamForExtrap(AliMUONTrackParam* trackParam, Double_t forwardBackward, Double_t *v3);
  static void RecoverTrackParam(Double_t *v3, Double_t Charge, AliMUONTrackParam* trackParam);
  
  static void AddMCSEffectInAbsorber(AliMUONTrackParam* trackParam, Double_t pathLength, Double_t f0, Double_t f1, Double_t f2);
  static void GetAbsorberCorrectionParam(Double_t trackXYZIn[3], Double_t trackXYZOut[3], Double_t pTotal, Double_t &pathLength,
					 Double_t &f0, Double_t &f1, Double_t &f2, Double_t &meanRho, Double_t &totalELoss);
  
  static Double_t BetheBloch(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ);
  
  static void ExtrapOneStepHelix(Double_t charge, Double_t step, Double_t *vect, Double_t *vout);
  static void ExtrapOneStepHelix3(Double_t field, Double_t step, Double_t *vect, Double_t *vout);

  static void GetField(Double_t *Position, Double_t *Field);
  
  ClassDef(AliMUONTrackExtrap, 0) // Tools for track extrapolation in ALICE dimuon spectrometer
};
	
#endif
