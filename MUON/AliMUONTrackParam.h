#ifndef ALIMUONTRACKPARAM_H
#define ALIMUONTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

#include <TROOT.h>

class AliMUONHitForRec;
class AliMUONSegment;

class AliMUONTrackParam : public TObject {
 public:
  AliMUONTrackParam(){
    // Constructor
    ;} // Constructor
  virtual ~AliMUONTrackParam(){
    // Destructor
    ;} // Destructor

  // Get and Set methods for data
  Double_t GetInverseBendingMomentum(void);
  void SetInverseBendingMomentum(Double_t InverseBendingMomentum);
  Double_t GetBendingSlope(void);
  void SetBendingSlope(Double_t BendingSlope);
  Double_t GetNonBendingSlope(void);
  void SetNonBendingSlope(Double_t NonBendingSlope);
  Double_t GetZ(void);
  void SetZ(Double_t Z);
  Double_t GetBendingCoor(void);
  void SetBendingCoor(Double_t BendingCoor);
  Double_t GetNonBendingCoor(void);
  void SetNonBendingCoor(Double_t NonBendingCoor);

  void ExtrapToZ(Double_t Z);
  void ExtrapToStation(Int_t Station, AliMUONTrackParam *TrackParam);
  void ExtrapToVertex();  // extrapolation to vertex through the absorber
  void BransonCorrection(); // makes Branson correction
  Double_t TotalMomentumEnergyLoss(Double_t rLimit, Double_t pTotal, Double_t theta, Double_t xEndAbsorber, Double_t yEndAbsorber); // returns total momentum after energy loss correction in the absorber

 protected:
 private:
  Double_t fInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1)
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
