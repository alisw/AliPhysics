#ifndef ALIMUONTRACKPARAM_H
#define ALIMUONTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

///////////////////////////////////////////////////
// Track parameters in ALICE dimuon spectrometer
///////////////////////////////////////////////////

#include <TObject.h>

class AliMUONTrackParam : public TObject {
 public:
  AliMUONTrackParam(){
    fInverseBendingMomentum = 0;
    fBendingSlope = 0;
    fNonBendingSlope = 0;
    fZ = 0;
    fBendingCoor = 0;
    fNonBendingCoor = 0;
    // Constructor
  } // Constructor
  virtual ~AliMUONTrackParam(){} // Destructor
  
  AliMUONTrackParam(const AliMUONTrackParam& rhs);// copy constructor (should be added per default !)
  AliMUONTrackParam& operator=(const  AliMUONTrackParam& rhs);// (should be added per default !)
  // Get and Set methods for data
  Double_t GetInverseBendingMomentum(void) const {return fInverseBendingMomentum;}
  void SetInverseBendingMomentum(Double_t InverseBendingMomentum) {fInverseBendingMomentum = InverseBendingMomentum;}
  Double_t GetBendingSlope(void) const {return fBendingSlope;}
  void SetBendingSlope(Double_t BendingSlope) {fBendingSlope = BendingSlope;}
  Double_t GetNonBendingSlope(void) const {return fNonBendingSlope;}
  void SetNonBendingSlope(Double_t NonBendingSlope) {fNonBendingSlope = NonBendingSlope;}
  Double_t GetZ(void) const {return fZ;}
  void SetZ(Double_t Z) {fZ = Z;}
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
  void SetBendingCoor(Double_t BendingCoor) {fBendingCoor = BendingCoor;}
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
  void SetNonBendingCoor(Double_t NonBendingCoor) {fNonBendingCoor = NonBendingCoor;}

  void ExtrapToZ(Double_t Z);
  void ExtrapToStation(Int_t Station, AliMUONTrackParam *TrackParam);
  void ExtrapToVertex();  // extrapolation to vertex through the absorber
  void BransonCorrection(); // makes Branson correction
  // returns total momentum after energy loss correction in the absorber
  Double_t TotalMomentumEnergyLoss(Double_t thetaLimit, Double_t pTotal, Double_t theta);
  void FieldCorrection(Double_t Z); // makes simple magnetic field correction through the absorber 

 protected:
 private:
  Double_t fInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  Double_t fBendingSlope; // Bending slope (cm ** -1)
  Double_t fNonBendingSlope; // Non bending slope (cm ** -1)
  Double_t fZ; // Z coordinate (cm)
  Double_t fBendingCoor; // bending coordinate (cm)
  Double_t fNonBendingCoor; // non bending coordinate (cm)

  void SetGeant3Parameters(Double_t *VGeant3, Double_t ForwardBackward);
  void GetFromGeant3Parameters(Double_t *VGeant3, Double_t Charge);

  ClassDef(AliMUONTrackParam, 1) // Track parameters in ALICE dimuon spectrometer
    };
	
#endif
