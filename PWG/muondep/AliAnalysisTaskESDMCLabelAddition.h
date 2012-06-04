#ifndef ALIANALYSISTASKESDMCLABELADDITION_H
#define ALIANALYSISTASKESDMCLABELADDITION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

#include <TString.h>
#include "AliAnalysisTaskSE.h"

class AliMUONTrack;
class AliMUONVTrackStore;

class AliAnalysisTaskESDMCLabelAddition : public AliAnalysisTaskSE
{
  
public:
  AliAnalysisTaskESDMCLabelAddition();
  AliAnalysisTaskESDMCLabelAddition(const char* name);
  virtual ~AliAnalysisTaskESDMCLabelAddition() {;}
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  /// Set the sigma cut to associate clusters with TrackRefs by position (instead of using recoParam)
  void SetExternalTrkSigmaCut(Double_t cut) { fExternalTrkSigmaCut = cut; }
  
  /// Set the sigma cut to associate trigger to triggerable track by position (instead of using recoParam)
  void SetExternalTrgSigmaCut(Double_t cut) { fExternalTrgSigmaCut = cut; }
  
  
private:
  
  AliAnalysisTaskESDMCLabelAddition(const AliAnalysisTaskESDMCLabelAddition&);
  AliAnalysisTaskESDMCLabelAddition& operator=(const AliAnalysisTaskESDMCLabelAddition&);
  
  // Check whether this combination of clusters correspond to a decaying particle or not
  Int_t IsDecay(Int_t nClusters, Int_t *chId, Int_t *labels, Bool_t &isReconstructible, Int_t &lastCh) const;
  
  // Try to match clusters between track and trackRef and add the corresponding MC labels to the arrays
  void AddCompatibleClusters(const AliMUONTrack &track, const AliMUONTrack &trackRef,
			     TArrayI *labels, Int_t *nLabels) const;
  
  // Check whether this track correspond to a decaying particle by using cluster MC labels
  Int_t IsDecayByLabel(const AliMUONTrack &track, Bool_t &isReconstructible, Int_t &lastCh) const;
  
  // Check whether this track correspond to a decaying particle by comparing clusters position
  Int_t IsDecayByPosition(const AliMUONTrack &track, const AliMUONVTrackStore &trackRefStore,
			  Bool_t &isReconstructible, Int_t &lastCh) const;
  
  TString  fDefaultStorage;       ///< location of the default OCDB storage
  UInt_t   fRequestedStationMask; //!< mask of requested stations
  Bool_t   fRequest2ChInSameSt45; //!< 2 fired chambers requested in the same station (4 or 5) or not
  Double_t fExternalTrkSigmaCut;  ///< sigma cut to associate clusters with TrackRefs (instead of using recoParam)
  Double_t fSigmaCut;             //!< sigma cut to associate clusters with TrackRefs
  Double_t fExternalTrgSigmaCut;  ///< sigma cut to associate trigger to triggerable track (instead of using recoParam)
  Double_t fSigmaCutTrig;         //!< sigma cut to associate trigger to triggerable track
      
  ClassDef(AliAnalysisTaskESDMCLabelAddition, 3); // Analysis task for standard ESD filtering
  
};

#endif

