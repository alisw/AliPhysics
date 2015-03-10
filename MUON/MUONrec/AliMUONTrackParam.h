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
#include <TMatrixD.h>

class AliMUONVCluster;

class AliMUONTrackParam : public TObject 
{
 public:
  AliMUONTrackParam(); // Constructor
  virtual ~AliMUONTrackParam(); // Destructor
  
  AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam);
  AliMUONTrackParam& operator=(const  AliMUONTrackParam& theMUONTrackParam);

  // Get and Set methods for data
	/// return Z coordinate (cm)
  Double_t GetZ() const {return fZ;}
	/// set Z coordinate (cm)
  void     SetZ(Double_t z) {fZ = z;}
	/// return non bending coordinate (cm)
  Double_t GetNonBendingCoor() const {return fParameters(0,0);}
	/// set non bending coordinate (cm)
  void     SetNonBendingCoor(Double_t nonBendingCoor) {fParameters(0,0) = nonBendingCoor;}
	/// return non bending slope (cm ** -1)
  Double_t GetNonBendingSlope() const {return fParameters(1,0);}
	/// set non bending slope (cm ** -1)
  void     SetNonBendingSlope(Double_t nonBendingSlope) {fParameters(1,0) = nonBendingSlope;}
	/// return bending coordinate (cm)
  Double_t GetBendingCoor() const {return fParameters(2,0);}
	/// set bending coordinate (cm)
  void     SetBendingCoor(Double_t bendingCoor) {fParameters(2,0) = bendingCoor;}
	/// return bending slope (cm ** -1)
  Double_t GetBendingSlope() const {return fParameters(3,0);}
	/// set bending slope (cm ** -1)
  void     SetBendingSlope(Double_t bendingSlope) {fParameters(3,0) = bendingSlope;}
	/// return inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  Double_t GetInverseBendingMomentum() const {return fParameters(4,0);}
	/// set inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)
  void     SetInverseBendingMomentum(Double_t inverseBendingMomentum) {fParameters(4,0) = inverseBendingMomentum;}
	/// return the charge (assumed forward motion)
  Double_t GetCharge() const {return TMath::Sign(1.,fParameters(4,0));}
	/// set the charge (assumed forward motion)
  void     SetCharge(Double_t charge) {if (charge*fParameters(4,0) < 0.) fParameters(4,0) *= -1.;}
  
  	/// return track parameters
  const TMatrixD& GetParameters() const {return fParameters;}
  	/// set track parameters
  void            SetParameters(const TMatrixD& parameters) {fParameters = parameters;}
  	/// add track parameters
  void            AddParameters(const TMatrixD& parameters) {fParameters += parameters;}
  
  Double_t Px() const;  // return px
  Double_t Py() const;  // return py
  Double_t Pz() const;  // return pz
  Double_t P()  const;  // return total momentum

	/// return kTRUE if the covariance matrix exist, kFALSE if not
  Bool_t    CovariancesExist() const {return (fCovariances) ? kTRUE : kFALSE;}
  
  const TMatrixD& GetCovariances() const;
  void            SetCovariances(const TMatrixD& covariances);
  void            SetCovariances(const Double_t matrix[5][5]);
  void            SetVariances(const Double_t matrix[5][5]);
  void            DeleteCovariances();
  
  const TMatrixD& GetPropagator() const;
  void            ResetPropagator();
  void            UpdatePropagator(const TMatrixD& propagator);
  
  const TMatrixD& GetExtrapParameters() const;
  void            SetExtrapParameters(const TMatrixD& parameters);
  
  const TMatrixD& GetExtrapCovariances() const;
  void            SetExtrapCovariances(const TMatrixD& covariances);
  
  const TMatrixD& GetSmoothParameters() const;
  void            SetSmoothParameters(const TMatrixD& parameters);
  
  const TMatrixD& GetSmoothCovariances() const;
  void            SetSmoothCovariances(const TMatrixD& covariances);
  
        /// get pointeur to associated cluster
  AliMUONVCluster* GetClusterPtr() const {return fClusterPtr;}
  void             SetClusterPtr(AliMUONVCluster* cluster, Bool_t owner = kFALSE);
  
  	/// return kTRUE if the associated cluster can be removed from the track it belongs to
  Bool_t IsRemovable() const {return fRemovable;}
  	/// set the flag telling whether the associated cluster can be removed from the track it belongs to or not
  void   SetRemovable(Bool_t removable) {fRemovable = removable;}
  
  /// return the chi2 of the track when the associated cluster was attached
  Double_t GetTrackChi2() const {return fTrackChi2;}
  	/// set the chi2 of the track when the associated cluster was attached
  void     SetTrackChi2(Double_t chi2) {fTrackChi2 = chi2;}
  	/// return the local chi2 of the associated cluster with respect to the track
  Double_t GetLocalChi2() const {return fLocalChi2;}
  	/// set the local chi2 of the associated cluster with respect to the track
  void     SetLocalChi2(Double_t chi2) {fLocalChi2 = chi2;}
  
	/// necessary for sorting TClonesArray of AliMUONTrackParam
  Bool_t IsSortable () const {return kTRUE;}
  Int_t Compare(const TObject* trackParam) const;

  Bool_t CompatibleTrackParam(const AliMUONTrackParam &trackParam, Double_t sigma2Cut, Double_t &normChi2) const;

  virtual void Print(Option_t* opt="") const;
 
  virtual void Clear(Option_t* opt="");

 private:
  
  Double_t fZ; ///< Z coordinate (cm)
  
  /// Track parameters ordered as follow:      <pre>
  /// X       = Non bending coordinate   (cm)
  /// SlopeX  = Non bending slope        (cm ** -1)
  /// Y       = Bending coordinate       (cm)
  /// SlopeY  = Bending slope            (cm ** -1)
  /// InvP_yz = Inverse bending momentum (GeV/c ** -1) times the charge (assumed forward motion)  </pre>
  TMatrixD fParameters; ///< \brief Track parameters 
  
  /// Covariance matrix of track parameters, ordered as follow:      <pre>
  ///    <X,X>      <X,SlopeX>        <X,Y>      <X,SlopeY>       <X,InvP_yz>
  /// <X,SlopeX>  <SlopeX,SlopeX>  <Y,SlopeX>  <SlopeX,SlopeY>  <SlopeX,InvP_yz>
  ///    <X,Y>      <Y,SlopeX>        <Y,Y>      <Y,SlopeY>       <Y,InvP_yz>
  /// <X,SlopeY>  <SlopeX,SlopeY>  <Y,SlopeY>  <SlopeY,SlopeY>  <SlopeY,InvP_yz>
  /// <X,InvP_yz> <SlopeX,InvP_yz> <Y,InvP_yz> <SlopeY,InvP_yz> <InvP_yz,InvP_yz>  </pre>
  mutable TMatrixD *fCovariances; ///< \brief Covariance matrix of track parameters 
  
  mutable TMatrixD *fPropagator;        //!<! Jacobian used to extrapolate the track parameters and covariances to the actual z position
  mutable TMatrixD *fExtrapParameters;  //!<! Track parameters extrapolated to the actual z position (not filtered by Kalman)
  mutable TMatrixD *fExtrapCovariances; //!<! Covariance matrix extrapolated to the actual z position (not filtered by Kalman)
  
  mutable TMatrixD *fSmoothParameters;  //!<! Track parameters obtained using smoother
  mutable TMatrixD *fSmoothCovariances; //!<! Covariance matrix obtained using smoother
  
  AliMUONVCluster *fClusterPtr; //!<! Pointer to associated cluster if any
  Bool_t           fOwnCluster; //!<! Ownership of the associated cluster
  
  Bool_t fRemovable; //!<! kTRUE if the associated cluster can be removed from the track it belongs to
  
  Double_t fTrackChi2; //!<! Chi2 of the track when the associated cluster was attached
  Double_t fLocalChi2; //!<! Local chi2 of the associated cluster with respect to the track
  
  ClassDef(AliMUONTrackParam, 4) // Track parameters in ALICE dimuon spectrometer
};
	
#endif
