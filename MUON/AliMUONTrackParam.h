#ifndef ALIMUONTRACKPARAM_H
#define ALIMUONTRACKPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrackParam
/// \brief Track parameters in ALICE dimuon spectrometer
///
////////////////////////////////////////////////////
/// Track parameters in ALICE dimuon spectrometer
////////////////////////////////////////////////////

#include <TObject.h>
#include <TMatrixDfwd.h>
#include "AliMUONHitForRec.h"

class AliESDMuonTrack;

class AliMUONTrackParam : public TObject 
{
 public:
  AliMUONTrackParam(); // Constructor
  virtual ~AliMUONTrackParam(); // Destructor
  
  AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam);
  AliMUONTrackParam& operator=(const  AliMUONTrackParam& theMUONTrackParam);

  void GetParamFrom(const AliESDMuonTrack& esdMuonTrack);
  void SetParamFor(AliESDMuonTrack& esdMuonTrack);

  // Get and Set methods for data
	/// return inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  Double_t GetInverseBendingMomentum(void) const {return fInverseBendingMomentum;}
	/// set inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  void     SetInverseBendingMomentum(Double_t inverseBendingMomentum) {fInverseBendingMomentum = inverseBendingMomentum;}
	/// return bending slope (cm ** -1)
  Double_t GetBendingSlope(void) const {return fBendingSlope;}
	/// set bending slope (cm ** -1)
  void     SetBendingSlope(Double_t bendingSlope) {fBendingSlope = bendingSlope;}
	/// return non bending slope (cm ** -1)
  Double_t GetNonBendingSlope(void) const {return fNonBendingSlope;}
	/// set non bending slope (cm ** -1)
  void     SetNonBendingSlope(Double_t nonBendingSlope) {fNonBendingSlope = nonBendingSlope;}
	/// return Z coordinate (cm)
  Double_t GetZ(void) const {return fZ;}
	/// set Z coordinate (cm)
  void     SetZ(Double_t z) {fZ = z;}
	/// return bending coordinate (cm)
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
	/// set bending coordinate (cm)
  void     SetBendingCoor(Double_t bendingCoor) {fBendingCoor = bendingCoor;}
	/// return non bending coordinate (cm)
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
	/// set non bending coordinate (cm)
  void     SetNonBendingCoor(Double_t nonBendingCoor) {fNonBendingCoor = nonBendingCoor;}
  
  void     SetTrackParam(AliMUONTrackParam& theMUONTrackParam);
  
  AliMUONHitForRec* GetHitForRecPtr(void) const;
	/// set pointeur to associated HitForRec
  void     SetHitForRecPtr(AliMUONHitForRec* hitForRec) {fHitForRecPtr = hitForRec;}
  
  Double_t Px() const;  // return px
  Double_t Py() const;  // return py
  Double_t Pz() const;  // return pz
  Double_t P()  const;  // return total momentum

	/// return kTRUE if the covariance matrix exist, kFALSE if not
  Bool_t    CovariancesExist(void) {return (fCovariances) ? kTRUE : kFALSE;}
  TMatrixD* GetCovariances(void);
  void      SetCovariances(TMatrixD* covariances);
  void      SetCovariances(Double_t matrix[5][5]);
  void      SetVariances(Double_t matrix[5][5]);
  void      DeleteCovariances(void);
  
  void EvalCovariances(AliMUONHitForRec* hit2);

	/// necessary for sorting TClonesArray of TrackHit's
  Bool_t IsSortable () const {return kTRUE;}
  Int_t Compare(const TObject* trackParam) const;

  virtual void Print(Option_t* opt="") const;
 

 private:
  // Parameters
  Double_t fNonBendingCoor; ///< Non bending coordinate (cm)
  Double_t fNonBendingSlope; ///< Non bending slope (cm ** -1)
  Double_t fBendingCoor; ///< Bending coordinate (cm)
  Double_t fBendingSlope; ///< Bending slope (cm ** -1)
  Double_t fInverseBendingMomentum; ///< Inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  Double_t fZ; ///< Z coordinate (cm)
  
  /// Covariance matrix of track parameters, ordered as follow:
  ///    <X,X>      <X,SlopeX>        <X,Y>      <X,SlopeY>       <X,InvP_yz>
  /// <X,SlopeX>  <SlopeX,SlopeX>  <Y,SlopeX>  <SlopeX,SlopeY>  <SlopeX,InvP_yz>
  ///    <X,Y>      <Y,SlopeX>        <Y,Y>      <Y,SlopeY>       <Y,InvP_yz>
  /// <X,SlopeY>  <SlopeX,SlopeY>  <Y,SlopeY>  <SlopeY,SlopeY>  <SlopeY,InvP_yz>
  /// <X,InvP_yz> <SlopeX,InvP_yz> <Y,InvP_yz> <SlopeY,InvP_yz> <InvP_yz,InvP_yz>
  TMatrixD *fCovariances; //!
  
  AliMUONHitForRec *fHitForRecPtr; //!< Pointer to associated HitForRec if any
  
  ClassDef(AliMUONTrackParam, 2) // Track parameters in ALICE dimuon spectrometer
};
	
#endif
