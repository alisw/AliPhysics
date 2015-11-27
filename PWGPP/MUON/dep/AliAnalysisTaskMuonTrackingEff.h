#ifndef ALIANALYSISTASKMUONTRACKINGEFF_H
#define ALIANALYSISTASKMUONTRACKINGEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup base
/// \class AliAnalysisTaskMuonTrackingEff
/// \brief tracking chamber efficiency from ESD data
//Author: Nicolas LE BRIS - SUBATECH Nantes

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"
#include "TString.h"

class TList;
class TObjArray;
class AliCounterCollection;
class AliMUONGeometryTransformer;
class AliMUONTrack;
class AliMUONTrackParam;
class AliMpArea;
class AliMpPad;
class AliMUON2DMap;


class AliAnalysisTaskMuonTrackingEff : public AliAnalysisTaskSE
{
  
 public:
  
  AliAnalysisTaskMuonTrackingEff();
  AliAnalysisTaskMuonTrackingEff(TString name);
  virtual ~AliAnalysisTaskMuonTrackingEff();

  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fOCDBpath = ocdbPath; }
  
  /// Set the OCDB path to the alignment file used in the reco (if not set use default storage)
  void SetAlignStorage(const char* ocdbPath) { fAlignOCDBpath = ocdbPath; }
  
  /// Set the OCDB path to the recoParam file used in the reco (if not set use default storage)
  void SetRecoParamStorage(const char* ocdbPath) { fRecoParamOCDBpath = ocdbPath; }
  
  /// Select tracks in the given centrality range
  void SelectCentrality(Double_t min, Double_t max) {fCentMin = min; fCentMax = max;}
  
  // set cuts to select tracks to be considered
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  
  // set default cuts to select tracks to be considered
  void SetDefaultMuonTrackCuts(Bool_t isMC);
  
  // get cuts to select tracks to be considered
  AliMuonTrackCuts* MuonTrackCuts() {return fMuonTrackCuts;}
  
  /// set the muon low pT cut
  void SetMuonPtCut(Double_t cut) {fPtCut = cut;}
  
  /// set the flag to select tracks using MC label
  void UseMCLabel(Bool_t flag = kTRUE) { fUseMCLabel = flag; }
  
  /// enable the display in the terminate
  void EnableDisplay(Bool_t flag = kTRUE) { fEnableDisplay = flag; }
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void NotifyRun();
  virtual void Terminate(Option_t *);
  
  
 private:
  
  /// Not implemented
  AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& rhs);
  /// Not implemented
  AliAnalysisTaskMuonTrackingEff& operator = (const AliAnalysisTaskMuonTrackingEff& rhs);
  
  // Identify clusters/chambers that can be removed from the track
  Bool_t TagRemovableClusters(AliMUONTrack &track, Bool_t removableChambers[10]);
  
  // Find which detection elements should have been hit and record the missing clusters
  void FindAndRecordMissingClusters(AliMUONTrackParam &param, Int_t chamber, Double_t trackInfo[6]);
  
  // Find the intersection point between the track (assuming straight line) and the DE in the global frame
  void Intersect(AliMUONTrackParam &param, Int_t deId, Double_t p[3]);
  
  // Check whether (global) area overlaps with the given DE
  Bool_t OverlapDE(AliMpArea &area, Int_t deId);
  
  // Register the cluster in the given stores
  void RecordCluster(Int_t chamber, Int_t deId, AliMpPad pad[2], Double_t trackInfo[6],
		     TString clusterKey, TList *chamberHistList, Bool_t recordChamber);
  
  /// Look for pads at the cluster's location
  Bool_t FindPads(Int_t deId, Double_t pos[3], AliMpPad pad[2]);
  
  
private:
  
  static const Int_t fgkNofDE[11];  ///< Total number of detection elements in each chamber
  static const Int_t fgkNofBusPath; ///< Total number of bus patches
  static const Int_t fgkNofManu;    ///< Total number of manus
  
  Bool_t   fOCDBLoaded;             //!< Determine if the OCDB and =geometry have been loaded
  TString  fOCDBpath;               ///< OCDB path
  TString  fAlignOCDBpath;          ///< OCDB path to the alignment file
  TString  fRecoParamOCDBpath;      ///< OCDB path to the recoParam file
  Double_t fCentMin;                ///< select centrality > fCentMin
  Double_t fCentMax;                ///< select centrality <= fCentMax
  AliMuonTrackCuts* fMuonTrackCuts; ///< cuts to select tracks to be considered
  Double_t fPtCut;                  ///< cut on minimum pt
  Bool_t   fUseMCLabel;             ///< select tracks using MC label
  Bool_t   fEnableDisplay;          ///< enable the display in the terminate

  AliMUONGeometryTransformer *fTransformer; //!< Transformer object
  
  TObjArray *fDEPlanes; //!< vectors (x0, y0, z0, a, b, c) defining the plane of each DE in the global frame
  
  AliCounterCollection* fClusters; //!< detected (all), accepted (for efficiency calculation) and expected clusters
  AliCounterCollection* fEvents; //!< number of analyzed events
  TList* fChamberTDHistList; //!< List of histograms of the tracks detected in the chambers.
  TList* fChamberTTHistList; //!< List of histograms of the tracks which have passed through the chambers.
  TList* fChamberSDHistList; //!< List of histograms of the tracks only detected by one chamber of the station.
  TList* fExtraHistList;     //!< List of extra histograms.


  ClassDef(AliAnalysisTaskMuonTrackingEff, 5)
  
};


//________________________________________________________________________
inline void AliAnalysisTaskMuonTrackingEff::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonTrackingEff::SetDefaultMuonTrackCuts(Bool_t isMC)
{
  /// set default cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts("stdCuts", "stdCuts");
  fMuonTrackCuts->SetAllowDefaultParams();
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
                                AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetIsMC(isMC);
}

#endif

