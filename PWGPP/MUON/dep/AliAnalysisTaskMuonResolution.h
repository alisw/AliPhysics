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
#include "AliAnalysisTaskSE.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"

class TH1;
class TH2;
class TGraphErrors;
class TObjArray;
class TList;
class AliMUONTrack;
class AliMUONTrackParam;
class AliMUONGeometryTransformer;
class AliMUONVCluster;

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
  
  void SetHalfChShift(Int_t hchId, Double_t valNB, Double_t valB);
  void SetHalfChShift(Double_t valNB[20], Double_t valB[20]);
  void GetHalfChShift(Double_t valNB[20], Double_t valB[20]) const;
  void ShiftHalfCh(Bool_t flag = kTRUE) { fShiftHalfCh = flag; }
  void PrintHalfChShift(Bool_t flag = kTRUE) { fPrintHalfChShift = flag; }
  
  void SetDEShift(Int_t iDE, Double_t valNB, Double_t valB);
  void SetDEShift(Double_t valNB[200], Double_t valB[200]);
  void GetDEShift(Double_t valNB[200], Double_t valB[200]) const;
  void ShiftDE(Bool_t flag = kTRUE) { fShiftDE = flag; }
  void PrintDEShift(Bool_t flag = kTRUE) { fPrintDEShift = flag; }
  
  /// set the minimum momentum value of the tracks used to compute the resolution
  void SetMinMomentum(Double_t val) { fMinMomentum = val; }
  
  /// set the minimum pT value of the tracks used to compute the resolution
  void SetMinPt(Double_t val) { fMinPt = val; }
  
  /// set the sign of the tracks used to compute the resolution
  void SetMuonSign(Short_t sign) { fSign = sign; }
  
  // set standard cuts to select events to be considered
  void SetMuonEventCuts(AliMuonEventCuts &eventCuts);
  
  // set standard cuts to select tracks to be considered
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  
  /// select only tracks with MC label or not
  void UseMCLabel(Bool_t flag = kTRUE) { fUseMCLabel = flag; }

  /// set the extrapolation mode to get the track parameters and covariances at a given cluster:
  /// 0 = extrapolate from the closest cluster; 1 = extrapolate from the previous cluster except between stations 2-3-4
  void SetExtrapMode(Int_t val) { fExtrapMode = val; }
  
  /// set the flag to add or not the systematic shifts of the residuals to the resolution
  void CorrectForSystematics(Bool_t flag = kTRUE) { fCorrectForSystematics = flag; }
  
  void SetAlignStorage(const char* ocdbPath, Int_t version = -1, Int_t subVersion = -1);
  
  void ReAlign(const char* oldAlignStorage = 0x0, Int_t oldVersion = -1, Int_t oldSubVersion = -1,
               const char* newAlignStorage = "", Int_t newVersion = -1, Int_t newSubVersion = -1);
  
  /// return the list of summary canvases
  TObjArray* GetCanvases() {return fCanvases;}
  
  /// set the flag to show the progression bar
  void ShowProgressBar(Bool_t flag = kTRUE) {fShowProgressBar = flag;}
  
  /// set the flag to print the cluster resolution per chamber/DE
  void PrintClusterRes(Bool_t perCh = kTRUE, Bool_t perDE = kFALSE) {fPrintClResPerCh = perCh; fPrintClResPerDE = perDE;}
  
  void FitResiduals(Bool_t flag = kTRUE);
  
  /// set the flag to remove mono-cathod clusters (either considering all the pads or only the ones directly below)
  void RemoveMonoCathodClusters(Bool_t flag = kTRUE, Bool_t checkAllPads = kTRUE) {fRemoveMonoCathCl = flag; fCheckAllPads = checkAllPads;}
  
  /// set the flag to improve the track before measuring the resolution
  void ImproveTracks(Bool_t flag = kTRUE) {fImproveTracks = flag;}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual Bool_t UserNotify();
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
  void GetMeanRMS(TH1* h, Double_t& mean, Double_t& meanErr,Double_t& rms, Double_t& rmsErr, TGraphErrors* gMean = 0x0,
		  TGraphErrors* gRMS = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE, Bool_t enableFit = kTRUE);
  void FillMeanSigmaClusterVsX(const TH2* hIn, const TH2* hOut, TGraphErrors* gMean, TGraphErrors* gSigma);
  void Cov2CovP(const AliMUONTrackParam &param, TMatrixD &covP);
  void CheckPads(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const;
  void CheckPadsBelow(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const;
  
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
    kResidualInChVsPClusterOut          = 20, ///< cluster-track residual-X/Y distribution in chamber i versus momentum (cluster not attached to the track)
    kResidualVsPClusterIn               = 40, ///< cluster-track residual-X/Y distribution integrated over chambers versus momentum (cluster attached to the track)
    kResidualVsPClusterOut              = 42  ///< cluster-track residual-X/Y distribution integrated over chambers versus momentum (cluster not attached to the track)
  };
  
  enum eResidualsVsCent {
    kResidualInChVsCentClusterIn        = 0,  ///< cluster-track residual-X/Y distribution in chamber i versus centrality (cluster attached to the track)
    kResidualInChVsCentClusterOut       = 20, ///< cluster-track residual-X/Y distribution in chamber i versus centrality (cluster not attached to the track)
    kResidualVsCentClusterIn            = 40, ///< cluster-track residual-X/Y distribution integrated over chambers versus centrality (cluster attached to the track)
    kResidualVsCentClusterOut           = 42  ///< cluster-track residual-X/Y distribution integrated over chambers versus centrality (cluster not attached to the track)
  };
  
  enum eResidualsVsAngle {
    kResidualInChVsAngleClusterIn       = 0,  ///< cluster-track residual-X/Y distribution in chamber i versus track angle in X/Y direction (cluster attached to the track)
    kResidualInChVsAngleClusterOut      = 20, ///< cluster-track residual-X/Y distribution in chamber i versus track angle in X/Y direction (cluster not attached to the track)
    kResidualVsAngleClusterIn           = 40, ///< cluster-track residual-X/Y distribution integrated over chambers versus track angle in X/Y direction (cluster attached to the track)
    kResidualVsAngleClusterOut          = 42, ///< cluster-track residual-X/Y distribution integrated over chambers versus track angle in X/Y direction (cluster not attached to the track)
    kResidualInHalfChVsAngleClusterIn   = 44  ///< cluster-track residual-X/Y distribution in half-chamber i versus track angle in X/Y direction (cluster attached to the track)
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
    kTrackResPerChMean                  = 12, ///< track X/Y-resolution per chamber
    kMCSPerChMean                       = 14, ///< MCS X/Y-dispersion of extrapolated track per chamber
    kClusterResPerCh                    = 16, ///< cluster X/Y-resolution per chamber
    kCalcClusterResPerCh                = 18, ///< calculated cluster X/Y-resolution per chamber
    kResidualPerDEMeanClusterIn         = 20, ///< cluster-track residual-X/Y per DE: mean (cluster in)
    kResidualPerDEMeanClusterOut        = 22, ///< cluster-track residual-X/Y per DE: mean (cluster out)
    kCombinedResidualPerDESigma         = 24, ///< combined cluster-track residual-X/Y per DE
    kClusterResPerDE                    = 26, ///< cluster X/Y-resolution per DE
    kResidualPerHalfChMeanClusterIn     = 28, ///< cluster-track residual-X/Y per half chamber: mean (cluster in)
    kResidualPerHalfChMeanClusterOut    = 30, ///< cluster-track residual-X/Y per half chamber: mean (cluster out)
    kCombinedResidualPerHalfChSigma     = 32, ///< combined cluster-track residual-X/Y per half chamber
    kClusterResPerHalfCh                = 34, ///< cluster X/Y-resolution per half chamber
    kResidualMeanClusterInVsP           = 36, ///< cluster-track residual-X/Y per chamber versus momentum: mean (cluster in)
    kCombinedResidualSigmaVsP           = 38, ///< cluster X/Y-resolution per chamber versus momentum
    kCombinedResidualAllChSigmaVsP      = 40, ///< cluster X/Y-resolution integrated over chambers versus momentum
    kResidualMeanClusterInVsCent        = 42, ///< cluster-track residual-X/Y per chamber versus centrality: mean (cluster in)
    kCombinedResidualSigmaVsCent        = 44, ///< cluster X/Y-resolution per chamber versus centrality
    kCombinedResidualAllChSigmaVsCent   = 46, ///< cluster X/Y-resolution integrated over chambers versus centrality
    kResidualMeanClusterInVsAngle       = 48, ///< cluster-track residual-X/Y per chamber versus track angle in X/Y direction: mean (cluster in)
    kCombinedResidualSigmaVsAngle       = 50, ///< cluster X/Y-resolution per chamber versus track angle in X/Y direction
    kCombinedResidualAllChSigmaVsAngle  = 52, ///< cluster X/Y-resolution integrated over chambers versus track angle in X/Y direction
    kHChResidualMeanClusterInVsAngle    = 54  ///< cluster-track residual-X/Y per half-chamber versus track angle in X/Y direction: mean (cluster in)
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
    kResPerDE                           = 1,  ///< summary canvas
    kResPerHalfCh                       = 2,  ///< summary canvas
    kResPerChVsP                        = 3,  ///< summary canvas
    kResPerChVsCent                     = 4,  ///< summary canvas
    kResPerChVsAngle                    = 5,  ///< summary canvas
    kShiftPerChVsP                      = 6,  ///< summary canvas
    kShiftPerChVsCent                   = 7,  ///< summary canvas
    kShiftPerChVsAngle                  = 8,  ///< summary canvas
    kDetailResPerCh                     = 9,  ///< summary canvas
    kDetailResPerHalfCh                 = 13, ///< summary canvas
    kDetailResPerDE                     = 17  ///< summary canvas
  };
  
  static const Int_t fgkMinEntries; ///< minimum number of entries needed to compute resolution
  
  TObjArray*  fResiduals;       //!< List of residual histos
  TObjArray*  fResidualsVsP;    //!< List of residual vs. p histos
  TObjArray*  fResidualsVsCent; //!< List of residual vs. centrality histos
  TObjArray*  fResidualsVsAngle;//!< List of residual vs. track angle histos
  TObjArray*  fLocalChi2;       //!< List of plots related to local chi2 per chamber/DE
  TObjArray*  fChamberRes;      //!< List of plots related to chamber/DE resolution
  TObjArray*  fTrackRes;        //!< List of plots related to track resolution (p, pT, ...)
  TObjArray*  fCanvases;        //!< List of canvases summarizing the results
  TObjArray*  fTmpHists;        //!< List of temporary histograms
  
  Double_t fClusterResNB[10];   ///< cluster resolution in non-bending direction
  Double_t fClusterResB[10];    ///< cluster resolution in bending direction
  
  Double_t fHalfChShiftNB[20];  ///< half-chamber deplacements in non-bending direction
  Double_t fHalfChShiftB[20];   ///< half-chamber deplacements in bending direction
  
  Double_t fDEShiftNB[200];     ///< DE deplacements in non-bending direction
  Double_t fDEShiftB[200];      ///< DE deplacements in bending direction
  
  TString  fDefaultStorage;        ///< location of the default OCDB storage
  Int_t    fNEvents;               //!< number of processed events
  Bool_t   fShowProgressBar;       ///< show the progression bar
  Bool_t   fPrintClResPerCh;       ///< print the cluster resolution per chamber
  Bool_t   fPrintClResPerDE;       ///< print the cluster resolution per DE
  TF1*     fGaus;                  ///< gaussian function to fit the residuals
  Double_t fMinMomentum;           ///< use only tracks with momentum higher than this value
  Double_t fMinPt;                 ///< use only tracks with pT higher than this value
  Short_t  fSign;                  ///< use only tracks of this sign
  Bool_t   fUseMCLabel;            ///< use only tracks with MC label or not
  Int_t    fExtrapMode;            ///< extrapolation mode to get the track parameters and covariances at a given cluster
  Bool_t   fCorrectForSystematics; ///< add or not the systematic shifts of the residuals to the resolution
  Bool_t   fRemoveMonoCathCl;      ///< remove or not the mono-cathod clusters
  Bool_t   fCheckAllPads;          ///< use all pads or only the ones directly below the cluster to look for mono-cathods
  Bool_t   fImproveTracks;         ///< flag telling whether to improve or not the track before measuring the resolution
  Bool_t   fShiftHalfCh;           ///< flag telling wether to displace half-chambers by fHalfChShift(N)B[i] or not
  Bool_t   fPrintHalfChShift;      ///< print the half-chamber displacements
  Bool_t   fShiftDE;               ///< flag telling wether to displace DEs by fDEShift(N)B[i] or not
  Bool_t   fPrintDEShift;          ///< print the DE displacements
  Bool_t   fOCDBLoaded;            //!< flag telling if the OCDB has been properly loaded or not
  Int_t    fNDE;                   //!< total number of DE
  Int_t    fDEIndices[1100];       //!< index of DE in histograms refered by ID
  Int_t    fDEIds[200];            //!< ID of DE refered by index in histograms
  Bool_t   fReAlign;               ///< flag telling whether to re-align the spectrometer or not before computing resolution
  TString  fOldAlignStorage;       ///< location of the OCDB storage where to find old MUON/Align/Data (use the default one if empty)
  Int_t    fOldAlignVersion;       ///< specific version of the old MUON/Align/Data/object to load
  Int_t    fOldAlignSubVersion;    ///< specific subversion of the old MUON/Align/Data/object to load
  TString  fNewAlignStorage;       ///< location of the OCDB storage where to find new MUON/Align/Data (use the default one if empty)
  Int_t    fNewAlignVersion;       ///< specific version of the new MUON/Align/Data/object to load
  Int_t    fNewAlignSubVersion;    ///< specific subversion of the new MUON/Align/Data/object to load
  AliMUONGeometryTransformer* fOldGeoTransformer; //!< geometry transformer used to recontruct the present data
  AliMUONGeometryTransformer* fNewGeoTransformer; //!< new geometry transformer containing the new alignment to be applied
  
  AliMuonEventCuts* fMuonEventCuts; ///< cuts to select events to be considered
  AliMuonTrackCuts* fMuonTrackCuts; ///< cuts to select tracks to be considered
  
  ClassDef(AliAnalysisTaskMuonResolution, 7); // chamber resolution analysis
};

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetStartingResolution(Int_t chId, Double_t valNB, Double_t valB)
{
  /// set chamber non-bending and bending resolutions
  if (chId < 0 || chId >= 10) return;
  fClusterResNB[chId] = valNB;
  fClusterResB[chId] = valB;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetStartingResolution(Double_t valNB[10], Double_t valB[10])
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 10; i++) {
    fClusterResNB[i] = valNB[i];
    fClusterResB[i] = valB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::GetStartingResolution(Double_t valNB[10], Double_t valB[10]) const
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 10; i++) {
    valNB[i] = fClusterResNB[i];
    valB[i] = fClusterResB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetHalfChShift(Int_t hchId, Double_t valNB, Double_t valB)
{
  /// set chamber non-bending and bending resolutions
  if (hchId < 0 || hchId >= 20) return;
  fHalfChShiftNB[hchId] = valNB;
  fHalfChShiftB[hchId] = valB;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetHalfChShift(Double_t valNB[20], Double_t valB[20])
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 20; i++) {
    fHalfChShiftNB[i] = valNB[i];
    fHalfChShiftB[i] = valB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::GetHalfChShift(Double_t valNB[20], Double_t valB[20]) const
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 20; i++) {
    valNB[i] = fHalfChShiftNB[i];
    valB[i] = fHalfChShiftB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetDEShift(Int_t iDE, Double_t valNB, Double_t valB)
{
  /// set chamber non-bending and bending resolutions
  if (iDE < 0 || iDE >= 200) return;
  fDEShiftNB[iDE] = valNB;
  fDEShiftB[iDE] = valB;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetDEShift(Double_t valNB[200], Double_t valB[200])
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 200; i++) {
    fDEShiftNB[i] = valNB[i];
    fDEShiftB[i] = valB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::GetDEShift(Double_t valNB[200], Double_t valB[200]) const
{
  /// set chambers non-bending and bending resolutions
  for (Int_t i = 0; i < 200; i++) {
    valNB[i] = fDEShiftNB[i];
    valB[i] = fDEShiftB[i];
  }
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetAlignStorage(const char* ocdbPath, Int_t version, Int_t subVersion)
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
inline void AliAnalysisTaskMuonResolution::ReAlign(const char* oldAlignStorage, Int_t oldVersion, Int_t oldSubVersion,
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

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::FitResiduals(Bool_t flag)
{
  /// set gaussian function to fit the residual distribution to extract the mean and the dispersion.
  /// if not set: take the mean and the RMS of the distribution
  delete fGaus;
  if (flag) fGaus = new TF1("fGaus","gaus");
  else fGaus = NULL;
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetMuonEventCuts(AliMuonEventCuts &eventCuts)
{
  /// set standard cuts to select events to be considered
  delete fMuonEventCuts;
  fMuonEventCuts = new AliMuonEventCuts(eventCuts);
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonResolution::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set standard cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}

#endif

