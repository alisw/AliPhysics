#ifndef ALIANALYSISTASKMUONRESOLUTION_H
#define ALIANALYSISTASKMUONRESOLUTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup muondep
/// \class AliAnalysisTaskMuonResolution
/// \brief Muon spectrometer resolution
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>
#include <TMatrixD.h>
#include <TF1.h>
#include "AliMUONConstants.h"
#include "AliAnalysisTaskSE.h"

class TH1;
class TH2;
class TGraphErrors;
class TObjArray;
class TList;
class AliMUONTrack;
class AliMUONTrackParam;
class AliMUONGeometryTransformer;

class AliAnalysisTaskMuonResolution : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonResolution();
  AliAnalysisTaskMuonResolution(const char *name);
  virtual ~AliAnalysisTaskMuonResolution();
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  void SetStartingResolution(Int_t chId, Double_t valNB, Double_t valB);
  void SetStartingResolution(Double_t valNB[10], Double_t valB[10]);
  void GetStartingResolution(Double_t valNB[10], Double_t valB[10]) const;
  
  /// set the minimum momentum value of the tracks used to compute the resolution
  void SetMinMomentum(Double_t val) { fMinMomentum = val; }
  
  /// set the flag to use only tracks passing the physics selection
  void SelectPhysics(Bool_t flag = kTRUE) {fSelectPhysics = flag;}
  
  /// set the flag to use only tracks matched with trigger or not
  void MatchTrigger(Bool_t flag = kTRUE) { fMatchTrig = flag; }
  
  /// set the flag to use only tracks passing the acceptance cuts (Rabs, eta)
  void ApplyAccCut(Bool_t flag = kTRUE) { fApplyAccCut = flag; }
  
  /// Select events belonging to at least one of the trigger classes selected by the mask to fill histograms:
  /// - if the physics selection is used, apply the mask to the trigger word returned by the physics selection
  /// - if not, apply the mask to the trigger word built by looking for triggers listed in "fSelectTriggerClass"
  void SelectTrigger(Bool_t flag = kTRUE, UInt_t mask = AliVEvent::kMUON) {fSelectTrigger = flag; fTriggerMask = mask;}
  
  /// set the extrapolation mode to get the track parameters and covariances at a given cluster:
  /// 0 = extrapolate from the closest cluster; 1 = extrapolate from the previous cluster except between stations 2-3-4
  void SetExtrapMode(Int_t val) { fExtrapMode = val; }
  
  /// set the flag to add or not the systematic shifts of the residuals to the resolution
  void CorrectForSystematics(Bool_t flag = kTRUE) { fCorrectForSystematics = flag; }
  
  void ReAlign(const char* oldAlignStorage = 0x0, const char* newAlignStorage = "");
  
  /// return the list of summary canvases
  TObjArray* GetCanvases() {return fCanvases;}
  
  /// set the flag to show the progression bar
  void ShowProgressBar(Bool_t flag = kTRUE) {fShowProgressBar = flag;}
  
  /// set the flag to print the cluster resolution per chamber/DE
  void PrintClusterRes(Bool_t perCh = kTRUE, Bool_t perDE = kFALSE) {fPrintClResPerCh = perCh; fPrintClResPerDE = perDE;}
  
  void FitResiduals(Bool_t flag = kTRUE);
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *);
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonResolution(const AliAnalysisTaskMuonResolution& rhs);
  /// Not implemented
  AliAnalysisTaskMuonResolution& operator = (const AliAnalysisTaskMuonResolution& rhs);
  
  void ModifyClusters(AliMUONTrack& track);
  void Zoom(TH1* h, Double_t fractionCut = 0.01);
  void ZoomLeft(TH1* h, Double_t fractionCut = 0.02);
  void ZoomRight(TH1* h, Double_t fractionCut = 0.02);
  void GetMean(TH1* h, Double_t& mean, Double_t& meanErr, TGraphErrors* g = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE, Bool_t enableFit = kTRUE);
  void GetRMS(TH1* h, Double_t& rms, Double_t& rmsErr, TGraphErrors* g = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE);
  void FillSigmaClusterVsP(const TH2* hIn, const TH2* hOut, TGraphErrors* g, Bool_t zoom = kTRUE);
  void Cov2CovP(const AliMUONTrackParam &param, TMatrixD &covP);
  UInt_t BuildTriggerWord(const TString& FiredTriggerClasses);
  
private:
  
  enum eResiduals {
    kResidualPerChClusterIn             = 0,  ///< cluster-track residual-X/Y distribution per chamber (cluster attached to the track)
    kResidualPerChClusterOut            = 2,  ///< cluster-track residual-X/Y distribution per chamber (cluster not attached to the track)
    kTrackResPerCh                      = 4,  ///< track resolution-X/Y per chamber
    kMCSPerCh                           = 6,  ///< MCS X/Y-dispersion of extrapolated track per chamber
    kClusterRes2PerCh                   = 8,  ///< cluster X/Y-resolution per chamber
    kResidualPerDEClusterIn             = 10, ///< cluster-track residual-X/Y distribution per DE (cluster attached to the track)
    kResidualPerDEClusterOut            = 12, ///< cluster-track residual-X/Y distribution per DE (cluster not attached to the track)
    kTrackResPerDE                      = 14, ///< track resolution-X/Y per DE
    kMCSPerDE                           = 16, ///< MCS X/Y-dispersion of extrapolated track per DE
    kResidualPerHalfChClusterIn         = 18, ///< cluster-track residual-X/Y distribution per half chamber (cluster attached to the track)
    kResidualPerHalfChClusterOut        = 20, ///< cluster-track residual-X/Y distribution per half chamber (cluster not attached to the track)
    kTrackResPerHalfCh                  = 22, ///< track resolution-X/Y per half chamber
    kMCSPerHalfCh                       = 24, ///< MCS X/Y-dispersion of extrapolated track per half chamber
    kLocalChi2PerCh                     = 26, ///< local chi2-X/Y/total distribution per chamber
    kLocalChi2PerDE                     = 29  ///< local chi2-X/Y/total distribution per DE
  };
  
  enum eResidualsVsP {
    kResidualInChVsPClusterIn           = 0,  ///< cluster-track residual-X/Y distribution in chamber i versus momentum (cluster attached to the track)
    kResidualInChVsPClusterOut          = 20  ///< cluster-track residual-X/Y distribution in chamber i versus momentum (cluster not attached to the track)
  };
  
  enum eLocalChi2 {
    kLocalChi2PerChMean                 = 0,  ///< local chi2-X/Y/total per chamber: mean
    kLocalChi2PerDEMean                 = 3   ///< local chi2-X/Y/total per DE: mean
  };
  
  enum eChamberRes {
    kResidualPerChMeanClusterIn         = 0,  ///< cluster-track residual-X/Y per chamber: mean (cluster in)
    kResidualPerChMeanClusterOut        = 2,  ///< cluster-track residual-X/Y per chamber: mean (cluster out)
    kResidualPerChSigmaClusterIn        = 4,  ///< cluster-track residual-X/Y per chamber: sigma (cluster in)
    kResidualPerChSigmaClusterOut       = 6,  ///< cluster-track residual-X/Y per chamber: sigma (cluster out)
    kResidualPerChDispersionClusterOut  = 8,  ///< cluster-track residual-X/Y per chamber: dispersion (cluster out)
    kCombinedResidualPerChSigma         = 10, ///< combined cluster-track residual-X/Y per chamber
    kCombinedResidualSigmaVsP           = 12, ///< cluster X/Y-resolution per chamber versus momentum
    kTrackResPerChMean                  = 14, ///< track X/Y-resolution per chamber
    kMCSPerChMean                       = 16, ///< MCS X/Y-dispersion of extrapolated track per chamber
    kClusterResPerCh                    = 18, ///< cluster X/Y-resolution per chamber
    kCalcClusterResPerCh                = 20, ///< calculated cluster X/Y-resolution per chamber
    kResidualPerDEMeanClusterIn         = 22, ///< cluster-track residual-X/Y per DE: mean (cluster in)
    kResidualPerDEMeanClusterOut        = 24, ///< cluster-track residual-X/Y per DE: mean (cluster out)
    kCombinedResidualPerDESigma         = 26, ///< combined cluster-track residual-X/Y per DE
    kClusterResPerDE                    = 28, ///< cluster X/Y-resolution per DE
    kResidualPerHalfChMeanClusterIn     = 30, ///< cluster-track residual-X/Y per half chamber: mean (cluster in)
    kResidualPerHalfChMeanClusterOut    = 32, ///< cluster-track residual-X/Y per half chamber: mean (cluster out)
    kCombinedResidualPerHalfChSigma     = 34, ///< combined cluster-track residual-X/Y per half chamber
    kClusterResPerHalfCh                = 36  ///< cluster X/Y-resolution per half chamber
  };
  
  enum eTrackRes {
    kUncorrPRes                         = 0,  ///< muon momentum reconstructed resolution at first cluster vs p
    kPRes                               = 1,  ///< muon momentum reconstructed resolution at vertex vs p
    kUncorrPtRes                        = 2,  ///< muon transverse momentum reconstructed resolution at first cluster vs p
    kPtRes                              = 3,  ///< muon transverse momentum reconstructed resolution at vertex vs p
    kUncorrSlopeRes                     = 4,  ///< muon slope-X/Y reconstructed resolution at first cluster vs p
    kSlopeRes                           = 6   ///< muon slope-X/Y reconstructed resolution at vertex vs p
  };
  
  enum eCanvases {
    kResPerCh                           = 0,  ///< summary canvas
    kResPerChVsP                        = 1,  ///< summary canvas
    kResPerDE                           = 2,  ///< summary canvas
    kResPerHalfCh                       = 3   ///< summary canvas
  };
  
  static const Int_t fgkMinEntries; ///< minimum number of entries needed to compute resolution
  
  TObjArray*  fResiduals;    //!< List of residual histos
  TObjArray*  fResidualsVsP; //!< List of residual vs. p histos
  TObjArray*  fLocalChi2;    //!< List of plots related to local chi2 per chamber/DE
  TObjArray*  fChamberRes;   //!< List of plots related to chamber/DE resolution
  TObjArray*  fTrackRes;     //!< List of plots related to track resolution (p, pT, ...)
  TObjArray*  fCanvases;     //!< List of canvases summarizing the results
  
  Double_t fClusterResNB[10]; ///< cluster resolution in non-bending direction
  Double_t fClusterResB[10];  ///< cluster resolution in bending direction
  
  TString  fDefaultStorage;        ///< location of the default OCDB storage
  Int_t    fNEvents;               //!< number of processed events
  Bool_t   fShowProgressBar;       ///< show the progression bar
  Bool_t   fPrintClResPerCh;       ///< print the cluster resolution per chamber
  Bool_t   fPrintClResPerDE;       ///< print the cluster resolution per DE
  TF1*     fGaus;                  ///< gaussian function to fit the residuals
  Double_t fMinMomentum;           ///< use only tracks with momentum higher than this value
  Bool_t   fSelectPhysics;         ///< use only tracks passing the physics selection
  Bool_t   fMatchTrig;             ///< use only tracks matched with trigger
  Bool_t   fApplyAccCut;           ///< use only tracks passing the acceptance cuts (Rabs, eta)
  Bool_t   fSelectTrigger;         ///< use only tracks passing the trigger selection
  UInt_t   fTriggerMask;           ///< trigger mask to be used when selecting tracks
  Int_t    fExtrapMode;            ///< extrapolation mode to get the track parameters and covariances at a given cluster
  Bool_t   fCorrectForSystematics; ///< add or not the systematic shifts of the residuals to the resolution
  Bool_t   fOCDBLoaded;            //!< flag telling if the OCDB has been properly loaded or not
  Int_t    fNDE;                   //!< total number of DE
  Int_t    fDEIndices[1100];       //!< index of DE in histograms refered by ID
  Int_t    fDEIds[200];            //!< ID of DE refered by index in histograms
  Bool_t   fReAlign;               ///< flag telling wether to re-align the spectrometer or not before computing resolution
  TString  fOldAlignStorage;       ///< location of the OCDB storage where to find old MUON/Align/Data (use the default one if empty)
  TString  fNewAlignStorage;       ///< location of the OCDB storage where to find new MUON/Align/Data (use the default one if empty)
  AliMUONGeometryTransformer* fOldGeoTransformer; //!< geometry transformer used to recontruct the present data
  AliMUONGeometryTransformer* fNewGeoTransformer; //!< new geometry transformer containing the new alignment to be applied
  
  TList* fSelectTriggerClass; //!< list of trigger class that can be selected to fill histograms
  
  ClassDef(AliAnalysisTaskMuonResolution, 2); // chamber resolution analysis
};

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetStartingResolution(Int_t chId, Double_t valNB, Double_t valB)
{
  /// set chamber non-bending and bending resolutions
  if (chId < 0 || chId >= AliMUONConstants::NTrackingCh()) return;
  fClusterResNB[chId] = valNB;
  fClusterResB[chId] = valB;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetStartingResolution(Double_t valNB[10], Double_t valB[10])
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    fClusterResNB[i] = valNB[i];
    fClusterResB[i] = valB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::GetStartingResolution(Double_t valNB[10], Double_t valB[10]) const
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    valNB[i] = fClusterResNB[i];
    valB[i] = fClusterResB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::ReAlign(const char* oldAlignStorage, const char* newAlignStorage)
{
  /// Set the flag to activate the re-alignment and the specific storage where to find the old/new alignment data.
  /// If oldAlignStorage = 0x0 we assume the spectrometer was not aligned before (default geometry)
  /// If old(new)AlignStorage = "" we assume the old(new) alignment data are in the default storage
  if (oldAlignStorage) fOldAlignStorage = oldAlignStorage;
  else fOldAlignStorage = "none";
  fNewAlignStorage = newAlignStorage;
  fReAlign = kTRUE;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::FitResiduals(Bool_t flag)
{
  /// set gaussian function to fit the residual distribution to extract the mean and the dispersion.
  /// if not set: take the mean and the RMS of the distribution
  if (fGaus) delete fGaus;
  if (flag) fGaus = new TF1("fGaus","gaus");
  else fGaus = NULL;
}

#endif

