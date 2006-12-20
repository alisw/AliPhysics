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
  
  static void ExtrapToZ(AliMUONTrackParam *trackParam, Double_t Z);
  static void ExtrapToZCov(AliMUONTrackParam* trackParam, Double_t zEnd);
  static void ExtrapToStation(AliMUONTrackParam *trackParamIn, Int_t station, AliMUONTrackParam *trackParamOut);
  static void ExtrapToVertexUncorrected(AliMUONTrackParam* trackParam, Double_t zVtx);
  static void ExtrapToVertex(AliMUONTrackParam *trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx);
  
  static void AddMCSEffectInTrackParamCov(AliMUONTrackParam *param, Double_t dZ, Double_t x0);
  
  static void ExtrapOneStepRungekutta(Double_t charge, Double_t step, Double_t* vect, Double_t* vout);
  
  
 private:
  static const AliMagF* fgkField;	  //!< field map
  static const Int_t    fgkMaxStepNumber; //!< Maximum number of steps for track extrapolation
  static const Double_t fgkStepLength;	  //!< Step lenght for track extrapolation

  // Functions
  AliMUONTrackExtrap(const AliMUONTrackExtrap& trackExtrap);
  AliMUONTrackExtrap& operator=(const AliMUONTrackExtrap& trackExtrap);

  static void ConvertTrackParamForExtrap(AliMUONTrackParam* trackParam, Double_t *v3, Double_t forwardBackward);
  static void RecoverTrackParam(Double_t *v3, Double_t Charge, AliMUONTrackParam* trackParam);
  
  static void BransonCorrection(AliMUONTrackParam *trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx);
  static Double_t TotalMomentumEnergyLoss(Double_t thetaLimit, Double_t pTotal, Double_t theta);
  static void FieldCorrection(AliMUONTrackParam *trackParam, Double_t Z);
  
  static void ExtrapOneStepHelix(Double_t charge, Double_t step, Double_t *vect, Double_t *vout);
  static void ExtrapOneStepHelix3(Double_t field, Double_t step, Double_t *vect, Double_t *vout);

  static void GetField(Double_t *Position, Double_t *Field);
  
  ClassDef(AliMUONTrackExtrap, 0) // Tools for track extrapolation in ALICE dimuon spectrometer
};
	
#endif
