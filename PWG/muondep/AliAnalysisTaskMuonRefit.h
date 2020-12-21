#ifndef ALIANALYSISTASKMUONREFIT_H
#define ALIANALYSISTASKMUONREFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup pwg_muondep_misc
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
  
  // Set OCDB path + version/subversion to find the alignment file used in the reco (if not set use default storage)
  void SetAlignStorage(const char* ocdbPath, Int_t version = -1, Int_t subVersion = -1);

  // re-align using default storages
  void ReAlignFromDefaultStorage();

  // Re-align clusters before refitting and set OCDB paths + versions/subversions to find the old/new alignment files
  void ReAlign(const char* oldAlignStorage = 0x0, Int_t oldVersion = -1, Int_t oldSubVersion = -1,
               const char* newAlignStorage = "", Int_t newVersion = -1, Int_t newSubVersion = -1);
  
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

  /// load Align storage version from ESD
  Bool_t GetAlignStorageFromESD();

private:
  
  Double_t fClusterResNB[10]; ///< cluster resolution in non-bending direction
  Double_t fClusterResB[10];  ///< cluster resolution in bending direction
  
  TString  fDefaultStorage;     ///< location of the default OCDB storage
  Bool_t   fImproveTracks;      ///< enable track improvement
  Double_t fSigmaCut;           ///< sigma cut for track improvement
  Double_t fSigmaCutForTrigger; ///< sigma cut for tracker/trigger track matching
  Bool_t   fReAlign;            ///< flag telling wether to re-align the spectrometer or not before refitting
  Bool_t   fUseDefaultAlignStorage; ///< flag telling wether to use default alignment storage for refitting. Old is set to what was used for the ESD. New is latest in the default storage
  TString  fOldAlignStorage;    ///< location of the OCDB storage where to find old MUON/Align/Data (use the default one if empty)
  Int_t    fOldAlignVersion;    ///< specific version of the old MUON/Align/Data/object to load
  Int_t    fOldAlignSubVersion; ///< specific subversion of the old MUON/Align/Data/object to load
  TString  fNewAlignStorage;    ///< location of the OCDB storage where to find new MUON/Align/Data (use the default one if empty)
  Int_t    fNewAlignVersion;    ///< specific version of the new MUON/Align/Data/object to load
  Int_t    fNewAlignSubVersion; ///< specific subversion of the new MUON/Align/Data/object to load
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
  
  ClassDef(AliAnalysisTaskMuonRefit, 3); // track refitter
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
inline void AliAnalysisTaskMuonRefit::SetAlignStorage(const char* ocdbPath, Int_t version, Int_t subVersion)
{
  /// Set the OCDB path + version/subversion to find the alignment file used in the reco.
  /// If ocdbPath = 0x0: do not apply any alignment (default geometry)
  /// If ocdbPath = "" : assume the alignment data are in the default storage
  /// If version = subversion = -1 the lastest object is loaded
  if (ocdbPath) {
    fNewAlignStorage = ocdbPath;
    fNewAlignVersion = version;
    fNewAlignSubVersion = subVersion;
  } else {
    fNewAlignStorage = "none";
    fNewAlignVersion = -1;
    fNewAlignSubVersion = -1;
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonRefit::ReAlignFromDefaultStorage()
{
  fUseDefaultAlignStorage = kTRUE;
  fReAlign = kTRUE;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonRefit::ReAlign(const char* oldAlignStorage, Int_t oldVersion, Int_t oldSubVersion,
                                              const char* newAlignStorage, Int_t newVersion, Int_t newSubVersion)
{
  /// Set the flag to activate the re-alignment and set the specific storages where to find
  /// the old/new alignment files of specified version/subversion.
  /// If old(new)AlignStorage = 0x0: do not apply any alignment (default geometry)
  /// If old(new)AlignStorage = "" : assume the old(new) alignment data are in the default storage
  /// If version = subversion = -1 the lastest object is loaded
  if (oldAlignStorage) {
    fOldAlignStorage = oldAlignStorage;
    fOldAlignVersion = oldVersion;
    fOldAlignSubVersion = oldSubVersion;
  } else {
    fOldAlignStorage = "none";
    fOldAlignVersion = -1;
    fOldAlignSubVersion = -1;
  }
  if (newAlignStorage) {
    fNewAlignStorage = newAlignStorage;
    fNewAlignVersion = newVersion;
    fNewAlignSubVersion = newSubVersion;
  } else {
    fNewAlignStorage = "none";
    fNewAlignVersion = -1;
    fNewAlignSubVersion = -1;
  }
  fReAlign = kTRUE;
}

#endif

