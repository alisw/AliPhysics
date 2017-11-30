#ifndef ALIMUONTRACKSMEARING_H
#define ALIMUONTRACKSMEARING_H

/// \class AliMuonTrackSmearing
/// \brief Parameterized muon resolution
///
/// The class allow propagate the measured cluster resolution
/// to the track momentum
///
/// \author Philippe Pillot <pillot@subatech.in2p3.fr>, Subatech
/// \author Diego Stocco <stocco@subatech.in2p3.fr>, Subatech
/// \date Jan 26, 2017

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class TF1;
//class AliVParticle;
//class AliMCEvent;
class TRootIOCtor;

class AliMuonTrackSmearing : public TObject {

public:
  AliMuonTrackSmearing ( TRootIOCtor* ioCtor );
  AliMuonTrackSmearing ( Int_t chosenFunc );
  virtual ~AliMuonTrackSmearing ();

  enum {
    kCrystalBall, // Crystal Ball function
    kBreitWigner, // Breight Wigner function
    kGaus         // Gaussian function
  };

//  void ClearRecoTrackList ();

//  AliVParticle* GetRecoTrack ( const AliVParticle* genParticle );
//  AliVParticle* GetRecoTrack ( const AliVParticle* recoTrack, const AliMCEvent* mcEvent );
  TLorentzVector GetRecoTrack ( Double_t pGen, Double_t etaGen, Double_t phiGen, Double_t chargeGen, Double_t &recoCharge, Double_t &rAbs );

  virtual void Print ( Option_t* option = "" ) const;

  /// Get the chosen function for smearing
  Int_t GetChosenFunc () { return fChosenFunc; }

  /// Get chamber resolution used during tracking
  Double_t GetSigmaTrk () { return fSigmaTrk; }
  
  /// Get sigma cut used during tracking
  Double_t GetSigmaTrkCut () { return fSigmaTrkCut; }
  
  /// Set chamber resolution along x on station 1
  void SetSigmaxChSt1 ( Double_t sigmaxChSt1 ) { fSigmaxChSt1 = sigmaxChSt1; }
  
  /// Set chamber resolution along y on station 1
  void SetSigmayChSt1 ( Double_t sigmayChSt1 ) { fSigmayChSt1 = sigmayChSt1; }
  
  /// Set overall chamber resolution along y
  void SetSigmayCh ( Double_t sigmayCh ) { fSigmayCh = sigmayCh; }
  
  // Set the tail parameters of the Crystal Ball functions
  void SetCrystalBallParams ( Double_t xChSt1Par1, Double_t xChSt1Par2, Double_t yChSt1Par1, Double_t yChSt1Par2, Double_t yChPar1, Double_t yChPar2 );

  /// Set systematic shift in alignment.
  /// The shift applies to mu plus. A shift with the opposite sign will be applied to mu minus
  void SetNSigmaShift ( Double_t nSigmaShift ) { fNSigmaShift = nSigmaShift; }

  /// Set chamber resolution used during tracking
  void SetSigmaTrk ( Double_t sigmaTrk ) { fSigmaTrk = sigmaTrk; }

  /// Set sigma cut used during tracking
  void SetSigmaTrkCut ( Double_t sigmaTrkCut ) { fSigmaTrkCut = sigmaTrkCut; }

  /// Set flag to tune kalman filter
  void SetTuneKalman ( Bool_t tuneKalman ) { fTuneKalman = tuneKalman; }
  
  // function computing the expected momentum resolution versus p given the smearing settings
  Double_t PResVsP( const Double_t *x, const Double_t *par );

  // function computing the expected slope resolution versus p given the smearing settings
  Double_t SlopeResVsP( const Double_t *x, const Double_t *par );
  
private:
  AliMuonTrackSmearing(const AliMuonTrackSmearing&);
  AliMuonTrackSmearing& operator=(const AliMuonTrackSmearing&);

  void SetupDefaultValues ();

  void ComputeRecoTrack ( Double_t pGen, Double_t etaGen, Double_t phiGen, Double_t chargeGen );
//  AliVParticle* GetRecoTrack ( const AliVParticle* genParticle, const AliVParticle* recoTrack );

  Double_t CrystalBallSymmetric ( Double_t *xx,Double_t *par );

  Double_t ELoss ( Double_t momentum, Double_t theta ) const;
  Double_t ELossFluctuation2 ( Double_t momentum, Double_t rhoZoverA ) const;
  Double_t FWHMELoss2 ( Double_t momentum, Double_t theta ) const;
  Double_t GenRndBreitWigner ( Double_t mean, Double_t sigma, Double_t max ) const;
  Double_t GenRndCrystalBall ( Double_t mean, Double_t sigma, Double_t tail1, Double_t tail2, Double_t max );
  Double_t GenRndGaus ( Double_t mean, Double_t sigma ) const;

  Double_t MCS2 ( Double_t momentum, Double_t dZ, Double_t x0 ) const;
  Double_t PToThetaDev ( Double_t momentum ) const;
  Double_t SigmaSlopeFromMCSInAbs2 ( Double_t momentum, Double_t theta ) const;
  Double_t SigmaSlopeFromMCSInCh2 ( Double_t momentum, Bool_t at1stCl, Double_t zB ) const;
  Double_t SigmaSlopeFromRes2 ( Bool_t bendingDir, Bool_t at1stCl, Double_t zB ) const;
  Double_t SigmaThetaDevFromMCS2 ( Double_t momentum ) const;
  Double_t SigmaThetaDevFromRes2() const;
  Double_t ThetaDevToP ( Double_t thetaDev ) const;

  Int_t fChosenFunc; ///< chosen function
  Double_t fSigmaTrk; ///< chamber resolution used during tracking
  Double_t fSigmaTrkCut; ///< sigma cut used during tracking
  Double_t fSigmaxChSt1; ///< x cluster sigma in station 1
  Double_t fSigmayChSt1; ///< y cluster sigma in station 1
  Double_t fSigmayCh; ///< y cluster sigma in average
  Double_t fNSigmaShift; ///< global shift in alignment
  Double_t fZB02; ///< branson plan position angle 0-2
  Double_t fZB23; ///< branson plan position angle 2-3
  Double_t fZB310; ///< branson plan position angle 3-10
  Double_t fMuMass; ///< muon mass
  Bool_t fTuneKalman; ///< tune the parameterization of MCS and energy loss to fit the momentum and angular resolution given by: kTRUE:  the Kalman filter ; kFALSE: the performance task
  std::vector<Double_t> fCrystalBallTails; ///< Tail parameters of Crystal ball functions
  TF1* fCrystalBall; //!<! Crystal Ball function
  Double_t fRecoCharge; //!<! Reconstructed track charge
  Double_t fRAbs; //!<! Reconstructed track transverse position at the end of the absorber (cm)
  TLorentzVector fRecoTrack; //!<! Reconstructed track from cluster resolution
//  std::vector<AliVParticle*> fRecoTrackList; //!<! Bookkeeping of produced tracks


  ClassDef(AliMuonTrackSmearing, 1); // Trigger chamber efficiencies
  /// \endcond
};

#endif
