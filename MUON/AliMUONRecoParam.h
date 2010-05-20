#ifndef ALIMUONRECOPARAM_H
#define ALIMUONRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONRecoParam
/// \brief Class with MUON reconstruction parameters
///
//  Author: Philippe Pillot

#include "AliDetectorRecoParam.h"
#include "TString.h"
#include <TMath.h>

class AliMUONRecoParam : public AliDetectorRecoParam
{
 public: 
  AliMUONRecoParam();
  virtual ~AliMUONRecoParam();
  
  static AliMUONRecoParam *GetLowFluxParam();
  static AliMUONRecoParam *GetHighFluxParam();
  static AliMUONRecoParam *GetCosmicParam();
  static AliMUONRecoParam *GetCalibrationParam();
  
  /// set the calibration mode (see GetCalibrationMode() for possible modes)
  void SetCalibrationMode(Option_t* mode) { fCalibrationMode = mode; fCalibrationMode.ToUpper();}
  
  Option_t* GetCalibrationMode() const;
  
  /// set the clustering (pre-clustering) mode
  void      SetClusteringMode(Option_t* mode) {fClusteringMode = mode; fClusteringMode.ToUpper();}
  /// get the clustering (pre-clustering) mode
  Option_t* GetClusteringMode() const {return fClusteringMode.Data();}
  
  /// Get the (truncated) average of sigmas of pedestal measurements, i.e. noise, of pads
  Double_t AverageNoisePadCharge() const { return fAverageNoisePadCharge; }
  /// Set the average of sigmas of pedestal measurements, i.e. noise, of pads
  void AverageNoisePadCharge(Double_t noise) { fAverageNoisePadCharge = noise; }
  
  /// Get the lowest charge we allow for pads
  Double_t LowestPadCharge() const { return fChargeSigmaCut*fAverageNoisePadCharge; }

  /// Get the cut applied to cut on cluster charge (the charge is cut if below fClusterChargeCut*LowestPadCharge())
  Double_t ClusterChargeCut() const { return fClusterChargeCut; }
  /// Set the cut applied to cut on cluster charge (the charge is cut if below fClusterChargeCut*LowestPadCharge())
  void ClusterChargeCut(Double_t n) { fClusterChargeCut=n; }
  
  /// Get the lowest possible cluster charge
  Double_t LowestClusterCharge() const { return ClusterChargeCut()*LowestPadCharge(); }
     
  /// set the tracking mode
  void      SetTrackingMode(Option_t* mode) {fTrackingMode = mode; fTrackingMode.ToUpper();}
  /// get the tracking mode
  Option_t* GetTrackingMode() const {return fTrackingMode.Data();}
  
  /// switch on/off the combined cluster/track reconstruction
  void      CombineClusterTrackReco(Bool_t flag) {fCombinedClusterTrackReco = flag;}
  /// return kTRUE/kFALSE if the combined cluster/track reconstruction is on/off
  Bool_t    CombineClusterTrackReco() const {return fCombinedClusterTrackReco;}
  
  /// save all cluster info (including pads) in ESD, for the given percentage of events
  void      SaveFullClusterInESD(Bool_t flag, Double_t percentOfEvent = 100.) {fSaveFullClusterInESD = flag;
                                 fPercentOfFullClusterInESD = (fSaveFullClusterInESD) ? percentOfEvent : 0.;}
  /// return kTRUE/kFALSE depending on whether we save all cluster info in ESD or not
  Bool_t    SaveFullClusterInESD() const {return fSaveFullClusterInESD;}
  /// return the percentage of events for which all cluster info are stored in ESD
  Double_t  GetPercentOfFullClusterInESD() const {return fPercentOfFullClusterInESD;}
  
  /// set the minimum value (GeV/c) of momentum in bending plane
  void     SetMinBendingMomentum(Double_t val) {fMinBendingMomentum = val;}
  /// return the minimum value (GeV/c) of momentum in bending plane
  Double_t GetMinBendingMomentum() const {return fMinBendingMomentum;}
  /// set the maximum value (GeV/c) of momentum in bending plane
  void     SetMaxBendingMomentum(Double_t val) {fMaxBendingMomentum = val;}
  /// return the maximum value (GeV/c) of momentum in bending plane
  Double_t GetMaxBendingMomentum() const {return fMaxBendingMomentum;}
  
  /// set the maximum value of the non bending slope
  void     SetMaxNonBendingSlope(Double_t val) {fMaxNonBendingSlope = val;}
  /// return the maximum value of the non bending slope
  Double_t GetMaxNonBendingSlope() const {return fMaxNonBendingSlope;}
  /// set the maximum value of the bending slope
  void     SetMaxBendingSlope(Double_t val) {fMaxBendingSlope = val;}
  /// return the maximum value of the bending slope
  Double_t GetMaxBendingSlope() const {return fMaxBendingSlope;}
  
  /// switch on/off the track selection according to their slope (instead of their impact parameter)
  void     SelectOnTrackSlope(Bool_t flag) {fSelectTrackOnSlope = flag;}
  /// return kTRUE/kFALSE if tracks are selected according to their slope/impact parameter
  Bool_t   SelectOnTrackSlope() const {return fSelectTrackOnSlope;}
  
  /// set the vertex dispersion (cm) in non bending plane
  void     SetNonBendingVertexDispersion(Double_t val) {fNonBendingVertexDispersion = val;} 
  /// return the vertex dispersion (cm) in non bending plane
  Double_t GetNonBendingVertexDispersion() const {return fNonBendingVertexDispersion;}
  /// set the vertex dispersion (cm) in bending plane
  void     SetBendingVertexDispersion(Double_t val) {fBendingVertexDispersion = val;} 
  /// return the vertex dispersion (cm) in bending plane
  Double_t GetBendingVertexDispersion() const {return fBendingVertexDispersion;}
  
  /// set the maximum distance to the track to search for compatible cluster(s) in non bending direction
  void     SetMaxNonBendingDistanceToTrack(Double_t val) {fMaxNonBendingDistanceToTrack = val;} 
  /// return the maximum distance to the track to search for compatible cluster(s) in non bending direction
  Double_t GetMaxNonBendingDistanceToTrack() const {return fMaxNonBendingDistanceToTrack;}
  /// set the maximum distance to the track to search for compatible cluster(s) in bending direction
  void     SetMaxBendingDistanceToTrack(Double_t val) {fMaxBendingDistanceToTrack = val;} 
  /// return the maximum distance to the track to search for compatible cluster(s) in bending direction
  Double_t GetMaxBendingDistanceToTrack() const {return fMaxBendingDistanceToTrack;}
  
  /// set the cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking
  void     SetSigmaCutForTracking(Double_t val) {fSigmaCutForTracking = val;} 
  /// return the cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking
  Double_t GetSigmaCutForTracking() const {return fSigmaCutForTracking;}
  
  /// switch on/off the track improvement and keep the default cut in sigma to apply on cluster (local chi2)
  void     ImproveTracks(Bool_t flag) {fImproveTracks = flag;} 
  /// switch on/off the track improvement and set the cut in sigma to apply on cluster (local chi2)
  void     ImproveTracks(Bool_t flag, Double_t sigmaCut) {fImproveTracks = flag; fSigmaCutForImprovement = sigmaCut;} 
  /// return kTRUE/kFALSE if the track improvement is switch on/off
  Bool_t   ImproveTracks() const {return fImproveTracks;}
  /// return the cut in sigma to apply on cluster (local chi2) during track improvement
  Double_t GetSigmaCutForImprovement() const {return fSigmaCutForImprovement;}
  
  /// set the cut in sigma to apply on track during trigger hit pattern search
  void     SetSigmaCutForTrigger(Double_t val) {fSigmaCutForTrigger = val; fMaxNormChi2MatchTrigger = val*val;} 
  /// return the cut in sigma to apply on track during trigger hit pattern search
  Double_t GetSigmaCutForTrigger() const {return fSigmaCutForTrigger;}
  /// set the cut in strips to apply on trigger track during trigger chamber efficiency
  void     SetStripCutForTrigger(Double_t val) {fStripCutForTrigger = val;} 
  /// return the cut in strips to apply on trigger track during trigger chamber efficiency
  Double_t GetStripCutForTrigger() const {return fStripCutForTrigger;}
  /// set the maximum search area in strips to apply on trigger track during trigger chamber efficiency
  void     SetMaxStripAreaForTrigger(Double_t val) {fMaxStripAreaForTrigger = val;} 
  /// return the maximum search area in strips to apply on trigger track during trigger chamber efficiency
  Double_t GetMaxStripAreaForTrigger() const {return fMaxStripAreaForTrigger;}
  
  /// return the maximum normalized chi2 of tracking/trigger track matching
  Double_t GetMaxNormChi2MatchTrigger() const {return fMaxNormChi2MatchTrigger;}
  
  /// switch on/off the tracking of all the possible candidates (track only the best one if switched off)
  void     TrackAllTracks(Bool_t flag) {fTrackAllTracks = flag;} 
  /// return kTRUE/kFALSE if the tracking of all the possible candidates is switched on/off
  Bool_t   TrackAllTracks() const {return fTrackAllTracks;}
  
  /// switch on/off the recovering of tracks being lost during reconstruction
  void     RecoverTracks(Bool_t flag) {fRecoverTracks = flag;} 
  /// return kTRUE/kFALSE if the recovering of tracks being lost during reconstruction is switched on/off
  Bool_t   RecoverTracks() const {return fRecoverTracks;}
  
  /// switch on/off the fast building of track candidates (assuming linear propagation between stations 4 and 5)
  void     MakeTrackCandidatesFast(Bool_t flag) {fMakeTrackCandidatesFast = flag;} 
  /// return kTRUE/kFALSE if the fast building of track candidates is switched on/off
  Bool_t   MakeTrackCandidatesFast() const {return fMakeTrackCandidatesFast;}
  
  /// switch on/off the building of track candidates starting from 1 cluster in each of the stations 4 and 5
  void     MakeMoreTrackCandidates(Bool_t flag) {fMakeMoreTrackCandidates = flag;} 
  /// return kTRUE/kFALSE if the building of extra track candidates is switched on/off
  Bool_t   MakeMoreTrackCandidates() const {return fMakeMoreTrackCandidates;}
  
  /// switch on/off the completion of reconstructed track
  void     ComplementTracks(Bool_t flag) {fComplementTracks = flag;} 
  /// return kTRUE/kFALSE if completion of the reconstructed track is switched on/off
  Bool_t   ComplementTracks() const {return fComplementTracks;}
  
  /// remove tracks sharing cluster in stations 1 or 2
  void     RemoveConnectedTracksInSt12(Bool_t flag) {fRemoveConnectedTracksInSt12 = flag;} 
  /// return kTRUE/kFALSE whether tracks sharing cluster in station 1 and 2 must be removed or not
  Bool_t   RemoveConnectedTracksInSt12() const {return fRemoveConnectedTracksInSt12;}
  
  /// switch on/off the use of the smoother
  void     UseSmoother(Bool_t flag) {fUseSmoother = flag;} 
  /// return kTRUE/kFALSE if the use of the smoother is switched on/off
  Bool_t   UseSmoother() const {return fUseSmoother;}
  
  /// switch on/off a chamber in the reconstruction
  void     UseChamber(Int_t iCh, Bool_t flag) {if (iCh >= 0 && iCh < 10) fUseChamber[iCh] = flag;}
  /// return kTRUE/kFALSE whether the chamber must be used or not
  Bool_t   UseChamber(Int_t iCh) const {return (iCh >= 0 && iCh < 10) ? fUseChamber[iCh] : kFALSE;}
  
  /// request or not at least one cluster in the station to validate the track
  void     RequestStation(Int_t iSt, Bool_t flag) {if (iSt >= 0 && iSt < 5) fRequestStation[iSt] = flag;}
  /// return kTRUE/kFALSE whether at least one cluster is requested in the station to validate the track
  Bool_t   RequestStation(Int_t iSt) const {return (iSt >= 0 && iSt < 5) ? fRequestStation[iSt] : kFALSE;}
  /// return an integer where first 5 bits are set to 1 if the corresponding station is requested
  UInt_t   RequestedStationMask() const;
  
  /// set the bypassSt45 value
  void   BypassSt45(Bool_t st4, Bool_t st5);
  
  /// return kTRUE if we should replace clusters in St 4 and 5 by generated clusters from trigger tracks
  Bool_t BypassSt45() const { return fBypassSt45==45; }
  
  /// return kTRUE if we should replace clusters in St 4 by generated clusters from trigger tracks
  Bool_t BypassSt4() const { return BypassSt45() || fBypassSt45==4 ; }
  
  /// return kTRUE if we should replace clusters in St 5 by generated clusters from trigger tracks
  Bool_t BypassSt5() const { return BypassSt45() || fBypassSt45==5 ; }
  
  /// Set Low and High threshold for St12 HV
  void    SetHVSt12Limits(float low, float high) { fHVSt12Limits[0]=low; fHVSt12Limits[1]=high; }
  /// Retrieve low limit for St12's HV
  Float_t HVSt12LowLimit() const { return fHVSt12Limits[0]; }
  /// Retrieve high limit for St12's HV
  Float_t HVSt12HighLimit() const { return fHVSt12Limits[1]; }
  
  /// Set Low and High threshold for St345 HV
  void    SetHVSt345Limits(float low, float high) { fHVSt345Limits[0]=low; fHVSt345Limits[1]=high; } 
  /// Retrieve low limit for St345's HV
  Float_t HVSt345LowLimit() const { return fHVSt345Limits[0]; }
  /// Retrieve high limit for St345's HV
  Float_t HVSt345HighLimit() const { return fHVSt345Limits[1]; }
  
  /// Set Low and High threshold for pedestal mean
  void    SetPedMeanLimits(float low, float high) { fPedMeanLimits[0]=low; fPedMeanLimits[1]=high; }
  /// Retrieve low limit of ped mean
  Float_t PedMeanLowLimit() const { return fPedMeanLimits[0]; }
  /// Retrieve high limit of ped mean
  Float_t PedMeanHighLimit() const { return fPedMeanLimits[1]; }
  
  /// Set Low and High threshold for pedestal sigma 
  void    SetPedSigmaLimits(float low, float high) { fPedSigmaLimits[0]=low; fPedSigmaLimits[1]=high; }
  /// Retrieve low limit of ped sigma
  Float_t PedSigmaLowLimit() const { return fPedSigmaLimits[0]; }
  /// Retrieve high limit of ped sigma
  Float_t PedSigmaHighLimit() const { return fPedSigmaLimits[1]; }
  
  /// Set Low and High threshold for gain a0 term
  void    SetGainA1Limits(float low, float high) { fGainA1Limits[0]=low; fGainA1Limits[1]=high; }
  /// Retrieve low limit of a1 (linear term) gain parameter
  Float_t GainA1LowLimit() const { return fGainA1Limits[0]; }
  /// Retrieve high limit of a1 (linear term) gain parameter
  Float_t GainA1HighLimit() const { return fGainA1Limits[1]; }
  
  /// Set Low and High threshold for gain a1 term
  void    SetGainA2Limits(float low, float high) { fGainA2Limits[0]=low; fGainA2Limits[1]=high; }
  /// Retrieve low limit of a2 (quadratic term) gain parameter
  Float_t GainA2LowLimit() const { return fGainA2Limits[0]; }
  /// Retrieve high limit of a2 (quadratic term) gain parameter
  Float_t GainA2HighLimit() const { return fGainA2Limits[1]; }
  
  /// Set Low and High threshold for gain threshold term
  void    SetGainThresLimits(float low, float high) { fGainThresLimits[0]=low; fGainThresLimits[1]=high; }
  /// Retrieve low limit on threshold gain parameter
  Float_t GainThresLowLimit() const { return fGainThresLimits[0]; }
  /// Retrieve high limit on threshold gain parameter
  Float_t GainThresHighLimit() const { return fGainThresLimits[1]; }
  
  /// Set the goodness mask (see AliMUONPadStatusMapMaker)
  void   SetPadGoodnessMask(UInt_t mask) { fPadGoodnessMask=mask; }
  /// Get the goodness mask
  UInt_t PadGoodnessMask() const { return fPadGoodnessMask; }
  
  /// Number of sigma cut we must apply when cutting on adc-ped
  Double_t ChargeSigmaCut() const { return fChargeSigmaCut; }
  
  /// Number of sigma cut we must apply when cutting on adc-ped
  void ChargeSigmaCut(Double_t value) { fChargeSigmaCut=value; }
  
  /// Set the default non bending resolution of chamber iCh
  void     SetDefaultNonBendingReso(Int_t iCh, Double_t val) {if (iCh >= 0 && iCh < 10) fDefaultNonBendingReso[iCh] = val;}
  /// Get the default non bending resolution of chamber iCh
  Double_t GetDefaultNonBendingReso(Int_t iCh) const {return (iCh >= 0 && iCh < 10) ? fDefaultNonBendingReso[iCh] : FLT_MAX;}
  /// Set the default bending resolution of chamber iCh
  void     SetDefaultBendingReso(Int_t iCh, Double_t val) {if (iCh >= 0 && iCh < 10) fDefaultBendingReso[iCh] = val;}
  /// Get the default bending resolution of chamber iCh
  Double_t GetDefaultBendingReso(Int_t iCh) const {return (iCh >= 0 && iCh < 10) ? fDefaultBendingReso[iCh] : FLT_MAX;}
  
  /// Set the maximum number of trigger tracks above which the tracking is cancelled
  void SetMaxTriggerTracks(Int_t maxTriggerTracks) {fMaxTriggerTracks = maxTriggerTracks;}
  /// Get the maximum number of trigger tracks above which the tracking is cancelled
  Int_t GetMaxTriggerTracks() const {return fMaxTriggerTracks;}
  
  /// Set the maximum number of track candidates above which the tracking abort
  void SetMaxTrackCandidates(Int_t maxTrackCandidates) {fMaxTrackCandidates = maxTrackCandidates;}
  /// Get the maximum number of track candidates above which the tracking abort
  Int_t GetMaxTrackCandidates() const {return fMaxTrackCandidates;}
  
  /// Set the limits for the acceptable manu occupancy
  void SetManuOccupancyLimits(float low, float high) { fManuOccupancyLimits[0]=low; fManuOccupancyLimits[1]=high; }
  /// Retrieve low value of manu occupancy limit
  Float_t ManuOccupancyLowLimit() const { return fManuOccupancyLimits[0]; }
  /// Retrieve high value of manu occupancy limit
  Float_t ManuOccupancyHighLimit() const { return fManuOccupancyLimits[1]; }

  /// Set the limits for the acceptable bp occupancy
  void SetBuspatchOccupancyLimits(float low, float high) { fBuspatchOccupancyLimits[0]=low; fBuspatchOccupancyLimits[1]=high; }
  /// Retrieve low value of bp occupancy limit
  Float_t BuspatchOccupancyLowLimit() const { return fBuspatchOccupancyLimits[0]; }
  /// Retrieve high value of bp occupancy limit
  Float_t BuspatchOccupancyHighLimit() const { return fBuspatchOccupancyLimits[1]; }

  /// Set the limits for the acceptable DE occupancy
  void SetDEOccupancyLimits(float low, float high) { fDEOccupancyLimits[0]=low; fDEOccupancyLimits[1]=high; }
  /// Retrieve low value of DE occupancy limit
  Float_t DEOccupancyLowLimit() const { return fDEOccupancyLimits[0]; }
  /// Retrieve high value of DE occupancy limit
  Float_t DEOccupancyHighLimit() const { return fDEOccupancyLimits[1]; }
  
  /// Set the fraction of buspatches outside the occupancy limits
  void SetFractionOfBuspatchOutsideOccupancyLimit(float v) { fFractionOfBuspatchOutsideOccupancyLimit = v; }
  /// Get the fraction of buspatches outside the occupancy limits
  Float_t FractionOfBuspatchOutsideOccupancyLimit() const { return fFractionOfBuspatchOutsideOccupancyLimit; }

  virtual void Print(Option_t *option = "") const;
  
  /// Get the max event size (soft limit)
  virtual Double_t EventSizeSoftLimit() const { return fEventSizeSoftLimit; }
  
  /// Get the max event size (hard limit)
  virtual Double_t EventSizeHardLimit() const { return fEventSizeHardLimit; }

  /// Set the max event size limits
  virtual void SetEventSizeLimits(Double_t soft, Double_t hard) { fEventSizeSoftLimit=soft; fEventSizeHardLimit=hard; }
  
  /// Get the percentage of token lost error we allow
  virtual Double_t TokenLostLimit() const { return fTokenLostLimit; }

  /// Set the percentage of token lost error we allow
  virtual void SetTokenLostLimit(Double_t limit) { fTokenLostLimit = limit; }

  /// Whether or not we try to recover corrupted raw data
  virtual Bool_t TryRecover() const { return fTryRecover; }

  /// Set the try recover corrupted raw data (use kTRUE only if you know what you are doing. Should be left to kFALSE by default)
  virtual void TryRecover(Bool_t flag) { fTryRecover = flag; }

  /// Create object ready to be put in OCDB
  static TObjArray* Create(const char* settings);
  
  /// Show what is the OCDB for that run
  static void Show(Int_t runNumber, const char* ocdbPath="raw://");
  
private:
  
  void SetDefaultLimits();
  
 private:
  
  /// clustering mode:  NOCLUSTERING, PRECLUSTER, PRECLUSTERV2, PRECLUSTERV3, COG, <pre>
  ///                   SIMPLEFIT, SIMPLEFITV3, MLEM:DRAW, MLEM, MLEMV2, MLEMV3   </pre>
  TString fClusteringMode; ///< \brief name of the clustering (+ pre-clustering) mode
  
  /// tracking mode: ORIGINAL, KALMAN
  TString fTrackingMode; ///< \brief name of the tracking mode
  
  Double32_t fMinBendingMomentum; ///< minimum value (GeV/c) of momentum in bending plane
  Double32_t fMaxBendingMomentum; ///< maximum value (GeV/c) of momentum in bending plane
  Double32_t fMaxNonBendingSlope; ///< maximum value of the non bending slope
  Double32_t fMaxBendingSlope;    ///< maximum value of the bending slope (used only if B = 0)
  
  Double32_t fNonBendingVertexDispersion; ///< vertex dispersion (cm) in non bending plane (used for original tracking only)
  Double32_t fBendingVertexDispersion;    ///< vertex dispersion (cm) in bending plane (used for original tracking only)
  
  Double32_t fMaxNonBendingDistanceToTrack; ///< maximum distance to the track to search for compatible cluster(s) in non bending direction
  Double32_t fMaxBendingDistanceToTrack;    ///< maximum distance to the track to search for compatible cluster(s) in bending direction
  
  Double32_t fSigmaCutForTracking; ///< cut in sigma to apply on cluster (local chi2) and track (global chi2) during tracking

  Double32_t fSigmaCutForImprovement; ///< cut in sigma to apply on cluster (local chi2) during track improvement
  
  Double32_t fSigmaCutForTrigger; ///< cut in sigma to apply on track during trigger hit pattern search

  Double32_t fStripCutForTrigger; ///< cut in strips to apply on trigger track during trigger chamber efficiency

  Double32_t fMaxStripAreaForTrigger; ///< max. search area in strips to apply on trigger track during trigger chamber efficiency
  
  Double32_t fMaxNormChi2MatchTrigger; ///< maximum normalized chi2 of tracking/trigger track matching
  
  Double32_t fPercentOfFullClusterInESD; ///< percentage of events for which all cluster info are stored in ESD
  
  Bool_t     fCombinedClusterTrackReco; ///< switch on/off the combined cluster/track reconstruction
  
  Bool_t     fTrackAllTracks; ///< kTRUE to track all the possible candidates; kFALSE to track only the best ones
  
  Bool_t     fRecoverTracks; ///< kTRUE to try to recover the tracks getting lost during reconstruction
  
  Bool_t     fMakeTrackCandidatesFast; ///< kTRUE to make candidate tracks assuming linear propagation between stations 4 and 5
  
  Bool_t     fMakeMoreTrackCandidates; ///< kTRUE to make candidate tracks starting from 1 cluster in each of the stations 4 and 5
  
  Bool_t     fComplementTracks; ///< kTRUE to try to complete the reconstructed tracks by adding missing clusters
  
  Bool_t     fImproveTracks; ///< kTRUE to try to improve the reconstructed tracks by removing bad clusters
  
  Bool_t     fUseSmoother; ///< kTRUE to use the smoother to compute track parameters/covariances and local chi2 at each cluster (used for Kalman tracking only)
  
  Bool_t     fSaveFullClusterInESD; ///< kTRUE to save all cluster info (including pads) in ESD
  
  /// calibration mode:  GAIN, NOGAIN, GAINCONSTANTCAPA, INJECTIONGAIN
  TString    fCalibrationMode; ///<\brief calibration mode
  
  Int_t      fBypassSt45; ///< non-zero to use trigger tracks to generate "fake" clusters in St 4 and 5. Can be 0, 4, 5 or 45 only
  
  Bool_t     fUseChamber[10]; ///< kTRUE to use the chamber i in the tracking algorithm
  
  Bool_t     fRequestStation[5]; ///< kTRUE to request at least one cluster in station i to validate the track
  
  Double32_t fGainA1Limits[2]; ///< Low and High threshold for gain a0 parameter
  Double32_t fGainA2Limits[2]; ///< Low and High threshold for gain a1 parameter
  Double32_t fGainThresLimits[2]; ///< Low and High threshold for gain threshold parameter
  Double32_t fHVSt12Limits[2]; ///< Low and High threshold for St12 HV
  Double32_t fHVSt345Limits[2]; ///< Low and High threshold for St345 HV
  Double32_t fPedMeanLimits[2]; ///< Low and High threshold for pedestal mean
  Double32_t fPedSigmaLimits[2]; ///< Low and High threshold for pedestal sigma
  
  UInt_t     fPadGoodnessMask; ///< goodness mask (see AliMUONPadStatusMaker)
  
  Double32_t fChargeSigmaCut; ///< number of sigma to cut on adc-ped 
  
  Double32_t fDefaultNonBendingReso[10]; ///< default chamber resolution in the non-bending direction
  Double32_t fDefaultBendingReso[10]; ///< default chamber resolution in the bending direction
  
  Bool_t     fRemoveConnectedTracksInSt12; ///< kTRUE to remove tracks sharing cluster in station 1 and 2
  
  Int_t      fMaxTriggerTracks; ///< maximum number of trigger tracks above which the tracking is cancelled
  Int_t      fMaxTrackCandidates; ///< maximum number of track candidates above which the tracking abort
  
  Bool_t     fSelectTrackOnSlope; ///< select track candidates according to their slope (instead of their impact parameter)
  
  Double32_t fManuOccupancyLimits[2]; ///< low and high thresholds for manu occupancy cut
  Double32_t fBuspatchOccupancyLimits[2]; ///< low and high thresholds for bus patch occupancy cut
  Double32_t fDEOccupancyLimits[2]; ///< low and high thresholds for DE occupancy cut

  Double32_t fMissingPadFractionLimit; ///< DEPRECATED
  Double32_t fFractionOfBuspatchOutsideOccupancyLimit; ///< above this limit, we consider we have too many buspatches out of the allowed occupancy range

  Double32_t fAverageNoisePadCharge; ///< the (truncated, typically at 10%) mean of the sigma of the pedestals, in femto-coulomb
  Double32_t fClusterChargeCut; ///< the cluster is cut if its charge is below fClusterChargeCut*LowestPadCharge()
  
  Double32_t fEventSizeSoftLimit; ///< (soft) limit on mean event size per event (KB)
  Double32_t fEventSizeHardLimit; ///< (hard) limit on mean event size per event (KB)
  
  Double32_t fTokenLostLimit; ///< limit on the fraction of token lost error per event we allow
  
  Bool_t     fTryRecover; ///< try to recover corrupted raw data
  
  // functions
  void SetLowFluxParam();
  void SetHighFluxParam();
  void SetCosmicParam();
  void SetCalibrationParam();
  
  ClassDef(AliMUONRecoParam,168) // MUON reco parameters
  // we're at 167 not because we had that many versions, but because at some point (version 15->16)
  // 166 was committed by error, and we did not to go reverse afterwards...
};

#endif

