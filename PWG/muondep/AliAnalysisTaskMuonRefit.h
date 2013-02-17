#ifndef ALIANALYSISTASKMUONREFIT_H
#define ALIANALYSISTASKMUONREFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup muondep
/// \class AliAnalysisTaskMuonRefit
/// \brief intermediate task for refitting
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>
#include "AliAnalysisTaskSE.h"

class AliMUONVCluster;
class AliMUONGeometryTransformer;
class AliMUONESDInterface;
class AliMUONRefitter;
class AliMUONTriggerCircuit;
class AliMUONTrackHitPattern;

class AliAnalysisTaskMuonRefit : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonRefit();
  AliAnalysisTaskMuonRefit(const char *name);
  virtual ~AliAnalysisTaskMuonRefit();
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  // reset the cluster resolution (by default use the one set in the recoParam)
  void ResetClusterResolution(Int_t chId, Double_t valNB, Double_t valB);
  void ResetClusterResolution(Double_t valNB[10], Double_t valB[10]);
  
  /// Enable track improvement (clusters/tracks that do not pass the sigma cut will be removed)
  /// If sigmaCut < 0: use the one set in the recoParam
  /// By default the "improveTrack" flag and the sigma cut are taken from the recoParam
  void ImproveTracks(Bool_t flag = kTRUE, Double_t sigmaCut = -1.) {fImproveTracks = flag; fSigmaCut = sigmaCut;} 
  
  /// Set the sigma cut for tracker/trigger track matching (by default use the one set in the recoParam)
  void SetSigmaCutForTrigger(Double_t val) {fSigmaCutForTrigger = val;}
  
  /// Set the OCDB path to the alignment file used in the reco (if not set use default storage)
  void SetAlignStorage(const char* ocdbPath) { fOldAlignStorage = ocdbPath; }
  
  // set the OCDB path to the alignment files to be used to realign the clusters before refitting
  void ReAlign(const char* oldAlignStorage = 0x0, const char* newAlignStorage = "");
  
  // Set the magnetic field map to be used during refitting
  void SetFieldPath(const char* field) { fField = field; }
  
  /// set the flag to remove mono-cathod clusters (either considering all the pads or only the ones directly below)
  void RemoveMonoCathodClusters(Bool_t flag = kTRUE, Bool_t checkAllPads = kTRUE) {fRemoveMonoCathCl = flag; fCheckAllPads = checkAllPads;}
  
  /// Tag the refitted tracks rejected by improvement (or no longer matched) instead of removing them (or the trigger part)
  /// Eventually keep the old track parameters if the refitted track no longer match the trigger
  void TagBadTracks(Bool_t tag = kTRUE, Bool_t keep = kTRUE) {fTagBadTracks = tag; fKeepOldParam = keep;}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonRefit(const AliAnalysisTaskMuonRefit& rhs);
  /// Not implemented
  AliAnalysisTaskMuonRefit& operator = (const AliAnalysisTaskMuonRefit& rhs);
  
  void ModifyCluster(AliMUONVCluster& cl);
  void CheckPads(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const;
  void CheckPadsBelow(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const;
  Bool_t SetMagField() const;
  
private:
  
  Double_t fClusterResNB[10]; ///< cluster resolution in non-bending direction
  Double_t fClusterResB[10];  ///< cluster resolution in bending direction
  
  TString  fDefaultStorage;     ///< location of the default OCDB storage
  Bool_t   fImproveTracks;      ///< enable track improvement
  Double_t fSigmaCut;           ///< sigma cut for track improvement
  Double_t fSigmaCutForTrigger; ///< sigma cut for tracker/trigger track matching
  Bool_t   fReAlign;            ///< flag telling wether to re-align the spectrometer or not before refitting
  TString  fOldAlignStorage;    ///< location of the OCDB storage where to find old MUON/Align/Data (use the default one if empty)
  TString  fNewAlignStorage;    ///< location of the OCDB storage where to find new MUON/Align/Data (use the default one if empty)
  AliMUONGeometryTransformer *fOldGeoTransformer; //!< geometry transformer used to recontruct the present data
  AliMUONGeometryTransformer *fNewGeoTransformer; //!< new geometry transformer containing the new alignment to be applied
  TString  fField;              ///< magnetic field map to be used during refitting
  Bool_t   fRemoveMonoCathCl;   ///< remove or not the mono-cathod clusters
  Bool_t   fCheckAllPads;       ///< use all pads or only the ones directly below the cluster to look for mono-cathods
  Bool_t   fTagBadTracks;       ///< tag the bad tracks instead of removing them
  Bool_t   fKeepOldParam;       ///< keep old track parameters if the refitted track no longer match the trigger

  AliMUONTriggerCircuit  *fTriggerCircuit;  //!< trigger circuit
  AliMUONESDInterface    *fESDInterface;    //!< esd interface to recover muon objects
  AliMUONRefitter        *fRefitter;        //!< refitter object
  AliMUONTrackHitPattern *fTrackHitPattern; //!< object to perform the tracker/trigger track matching
  
  ClassDef(AliAnalysisTaskMuonRefit, 2); // track refitter
};

//________________________________________________________________________
inline void AliAnalysisTaskMuonRefit::ResetClusterResolution(Int_t chId, Double_t valNB, Double_t valB)
{
  /// set chamber non-bending and bending resolutions
  if (chId < 0 || chId >= 10) return;
  fClusterResNB[chId] = valNB;
  fClusterResB[chId] = valB;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonRefit::ResetClusterResolution(Double_t valNB[10], Double_t valB[10])
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 10; i++) {
    fClusterResNB[i] = valNB[i];
    fClusterResB[i] = valB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonRefit::ReAlign(const char* oldAlignStorage, const char* newAlignStorage)
{
  /// Set the flag to activate the re-alignment and the specific storage where to find the old/new alignment data.
  /// If oldAlignStorage = 0x0 we assume the spectrometer was not aligned before (default geometry)
  /// If old(new)AlignStorage = "" we assume the old(new) alignment data are in the default storage
  if (oldAlignStorage) fOldAlignStorage = oldAlignStorage;
  else fOldAlignStorage = "none";
  fNewAlignStorage = newAlignStorage;
  fReAlign = kTRUE;
}

#endif

