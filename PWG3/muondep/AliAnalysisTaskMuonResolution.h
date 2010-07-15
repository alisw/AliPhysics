#ifndef ALIANALYSISTASKMUONRESOLUTION_H
#define ALIANALYSISTASKMUONRESOLUTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskMuonResolution
/// \brief Muon spectrometer resolution
//Author: Philippe Pillot - SUBATECH Nantes

#include <TString.h>
#include <TMatrixD.h>
#include "AliMUONConstants.h"

class TH1;
class TH2;
class TGraphErrors;
class TObjArray;
class AliMUONTrack;
class AliMUONTrackParam;
class AliMUONGeometryTransformer;

class AliAnalysisTaskMuonResolution : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonResolution();
  AliAnalysisTaskMuonResolution(const char *name);
  virtual ~AliAnalysisTaskMuonResolution();
  
  void SetStartingResolution(Int_t chId, Double_t valNB, Double_t valB);
  void SetStartingResolution(Double_t valNB[10], Double_t valB[10]);
  void GetStartingResolution(Double_t valNB[10], Double_t valB[10]);
  
  /// set the minimum momentum value of the tracks used to compute the resolution
  void SetMinMomentum(Double_t val) { fMinMomentum = val; }
  
  /// set the flag to use only tracks passing the physics selection
  void SelectPhysics(Bool_t flag = kTRUE) {fSelectPhysics = flag;}
  
  /// set the flag to use only tracks matched with trigger or not
  void MatchTrigger(Bool_t flag = kTRUE) { fMatchTrig = flag; }
  
  /// set the extrapolation mode to get the track parameters and covariances at a given cluster:
  /// 0 = extrapolate from the closest cluster; 1 = extrapolate from the previous cluster except between stations 2-3-4
  void SetExtrapMode(Int_t val) { fExtrapMode = val; }
  
  /// set the flag to add or not the systematic shifts of the residuals to the resolution
  void CorrectForSystematics(Bool_t flag = kTRUE) { fCorrectForSystematics = flag; }
  
  void ReAlign(const char* oldAlignStorage = 0x0, const char* newAlignStorage = "");
  
  /// return the list of summary canvases
  TObjArray* GetCanvases() {return fCanvases;}
  
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
  void GetMean(TH1* h, Double_t& mean, Double_t& meanErr, TGraphErrors* g = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE);
  void GetRMS(TH1* h, Double_t& rms, Double_t& rmsErr, TGraphErrors* g = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE);
  void FillSigmaClusterVsP(TH2* hIn, TH2* hOut, TGraphErrors* g, Bool_t zoom = kTRUE);
  void Cov2CovP(const AliMUONTrackParam &param, TMatrixD &covP);
  
private:
  
  enum outputIndices {
    kResidualPerCh_ClusterIn            = 0,  ///< cluster-track residual-X/Y distribution per chamber (cluster attached to the track)
    kResidualPerCh_ClusterOut           = 2,  ///< cluster-track residual-X/Y distribution per chamber (cluster not attached to the track)
    kTrackResPerCh                      = 4,  ///< track resolution-X/Y per chamber
    kMCSPerCh                           = 6,  ///< MCS X/Y-dispersion of extrapolated track per chamber
    kClusterRes2PerCh                   = 8,  ///< cluster X/Y-resolution per chamber
    kResidualInChVsP_ClusterIn          = 10, ///< cluster-track residual-X/Y distribution in chamber i versus momentum (cluster attached to the track)
    kResidualInChVsP_ClusterOut         = 30, ///< cluster-track residual-X/Y distribution in chamber i versus momentum (cluster not attached to the track)
    kResidualPerDE_ClusterIn            = 50, ///< cluster-track residual-X/Y distribution per DE (cluster attached to the track)
    kResidualPerDE_ClusterOut           = 52, ///< cluster-track residual-X/Y distribution per DE (cluster not attached to the track)
    kTrackResPerDE                      = 54, ///< track resolution-X/Y per DE
    kMCSPerDE                           = 56, ///< MCS X/Y-dispersion of extrapolated track per DE
    kResidualPerHalfCh_ClusterIn        = 58, ///< cluster-track residual-X/Y distribution per half chamber (cluster attached to the track)
    kResidualPerHalfCh_ClusterOut       = 60, ///< cluster-track residual-X/Y distribution per half chamber (cluster not attached to the track)
    kTrackResPerHalfCh                  = 62, ///< track resolution-X/Y per half chamber
    kMCSPerHalfCh                       = 64, ///< MCS X/Y-dispersion of extrapolated track per half chamber
    
    kLocalChi2PerCh                     = 100, ///< local chi2-X/Y/total distribution per chamber
    kLocalChi2PerDE                     = 103, ///< local chi2-X/Y/total distribution per DE
    kLocalChi2PerChMean                 = 106, ///< local chi2-X/Y/total per chamber: mean
    kLocalChi2PerDEMean                 = 109, ///< local chi2-X/Y/total per DE: mean
    
    kResidualPerChMean_ClusterIn        = 150, ///< cluster-track residual-X/Y per chamber: mean (cluster in)
    kResidualPerChMean_ClusterOut       = 152, ///< cluster-track residual-X/Y per chamber: mean (cluster out)
    kResidualPerChSigma_ClusterIn       = 154, ///< cluster-track residual-X/Y per chamber: sigma (cluster in)
    kResidualPerChSigma_ClusterOut      = 156, ///< cluster-track residual-X/Y per chamber: sigma (cluster out)
    kResidualPerChDispersion_ClusterOut = 158, ///< cluster-track residual-X/Y per chamber: dispersion (cluster out)
    kCombinedResidualPerChSigma         = 160, ///< combined cluster-track residual-X/Y per chamber
    kCombinedResidualSigmaVsP           = 162, ///< cluster X/Y-resolution per chamber versus momentum
    kTrackResPerChMean                  = 164, ///< track X/Y-resolution per chamber
    kMCSPerChMean                       = 166, ///< MCS X/Y-dispersion of extrapolated track per chamber
    kClusterResPerCh                    = 168, ///< cluster X/Y-resolution per chamber
    kCalcClusterResPerCh                = 170, ///< calculated cluster X/Y-resolution per chamber
    kResidualPerDEMean_ClusterIn        = 172, ///< cluster-track residual-X/Y per DE: mean (cluster in)
    kResidualPerDEMean_ClusterOut       = 174, ///< cluster-track residual-X/Y per DE: mean (cluster out)
    kCombinedResidualPerDESigma         = 176, ///< combined cluster-track residual-X/Y per DE
    kClusterResPerDE                    = 178, ///< cluster X/Y-resolution per DE
    kResidualPerHalfChMean_ClusterIn    = 180, ///< cluster-track residual-X/Y per half chamber: mean (cluster in)
    kResidualPerHalfChMean_ClusterOut   = 182, ///< cluster-track residual-X/Y per half chamber: mean (cluster out)
    kCombinedResidualPerHalfChSigma     = 184, ///< combined cluster-track residual-X/Y per half chamber
    kClusterResPerHalfCh                = 186, ///< cluster X/Y-resolution per half chamber
    
    kUncorrPRes                         = 250, ///< muon momentum reconstructed resolution at first cluster vs p
    kPRes                               = 251, ///< muon momentum reconstructed resolution at vertex vs p
    kUncorrPtRes                        = 252, ///< muon transverse momentum reconstructed resolution at first cluster vs p
    kPtRes                              = 253, ///< muon transverse momentum reconstructed resolution at vertex vs p
    kUncorrSlopeRes                     = 254, ///< muon slope-X/Y reconstructed resolution at first cluster vs p
    kSlopeRes                           = 256, ///< muon slope-X/Y reconstructed resolution at vertex vs p
    
    kResPerCh                           = 300, ///< summary canvas
    kResPerChVsP                        = 301, ///< summary canvas
    kResPerDE                           = 302, ///< summary canvas
    kResPerHalfCh                       = 303  ///< summary canvas
  };
  
  static const Int_t fgkMinEntries; ///< minimum number of entries needed to compute resolution
  
  TObjArray*  fResiduals;    //!< List of residual histos
  TObjArray*  fResidualsVsP; //!< List of residual vs. p histos
  TObjArray*  fLocalChi2;    //!< List of plots related to local chi2 per chamber/DE
  TObjArray*  fChamberRes;   //!< List of plots related to chamber/DE resolution
  TObjArray*  fTrackRes;     //!< List of plots related to track resolution (p, pT, ...)
  TObjArray*  fCanvases;     //!< List of canvases summarizing te results
  
  Double_t fClusterResNB[10]; ///< cluster resolution in non-bending direction
  Double_t fClusterResB[10];  ///< cluster resolution in bending direction
  
  Int_t    fNEvents;               //!< number of processed events
  Double_t fMinMomentum;           ///< use only tracks with momentum higher than this value
  Bool_t   fSelectPhysics;         ///< use only tracks passing the physics selection
  Bool_t   fMatchTrig;             ///< use only tracks matched with trigger
  /// extrapolation mode to get the track parameters and covariances at a given cluster:
  /// 0 = extrapolate from the closest cluster; 1 = extrapolate from the previous cluster except between stations 2-3-4
  Int_t    fExtrapMode;
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
  
  ClassDef(AliAnalysisTaskMuonResolution, 1); // chamber resolution analysis
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
inline void AliAnalysisTaskMuonResolution::GetStartingResolution(Double_t valNB[10], Double_t valB[10])
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

#endif

