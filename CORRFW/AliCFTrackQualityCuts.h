/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// The class AliCFTrackQualityCuts is designed to select reconstructed tracks
// of high quality and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, reconstructed tracks, it is applicable on
// ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - number of clusters in the TPC
// - number of clusters in the ITS
// - chi2 / cluster in the TPC
// - chi2 / cluster in the ITS
// - covariance matrix diagonal elements
// - track status (cf AliESDtrack.h)
//
// The cut values for these cuts are set with the corresponding set functions.
// All cut classes provided by the correction framework are supposed to be
// added in the Analysis Framwork's class AliAnalysisFilter and applied by
// the filter via a loop.
//
// author: I. Kraus (Ingrid.Kraus@cern.ch)
// idea taken form
// AliESDtrackCuts writte by Jan Fiete Grosse-Oetringhaus and
// AliRsnDaughterCut class written by A. Pulvirenti.

#ifndef ALICFTRACKQUALITYCUTS_H
#define ALICFTRACKQUALITYCUTS_H

#include "AliCFCutBase.h"

class TH2F ;
class TH1F ;
class TBits;
class AliESDtrack ;

class AliCFTrackQualityCuts : public AliCFCutBase
{
 public :
  AliCFTrackQualityCuts() ;
  AliCFTrackQualityCuts(Char_t* name, Char_t* title) ;
  AliCFTrackQualityCuts(const AliCFTrackQualityCuts& c) ;
  AliCFTrackQualityCuts& operator=(const AliCFTrackQualityCuts& c) ;
  ~AliCFTrackQualityCuts();
  void Copy(TObject &c) const;

  Bool_t IsSelected(TObject* obj);
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  // cut value setter
  void SetMinNClusterTPC(Int_t cluster=-1)		{fMinNClusterTPC = cluster;}
  void SetMinNClusterITS(Int_t cluster=-1)		{fMinNClusterITS = cluster;}
  void SetMaxChi2PerClusterTPC(Double_t chi=1e99)	{fMaxChi2PerClusterTPC = chi;}
  void SetMaxChi2PerClusterITS(Double_t chi=1e99)	{fMaxChi2PerClusterITS = chi;}
  void SetMaxCovDiagonalElements(Float_t c1=1e99, Float_t c2=1e99, Float_t c3=1e99, Float_t c4=1e99, Float_t c5=1e99) 
{fCovariance11Max=c1;fCovariance22Max=c2;fCovariance33Max=c3;fCovariance44Max=c4;fCovariance55Max=c5;}
  void SetStatus(ULong_t status=0) {fStatus = status ;}

  // QA histograms
  void DrawHistograms(Bool_t drawLogScale=kTRUE);
  void SaveHistograms(const Char_t* dir = 0);
  void AddQAHistograms(TList *qaList);
  // QA histogram setter
  // please use indices from the enumeration below
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins);
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax);

  // indeces/counters for single selections
  enum { 
    kCutClusterTPC=0,	// number of clusters in TPC
    kCutClusterITS,	// number of clusters in ITS
    kCutChi2TPC,	// chi2 per cluster in TPC
    kCutChi2ITS,	// chi2 per cluster in ITS
    kCutCovElement11,	// diagonal element 11 of covariance matrix
    kCutCovElement22,	// diagonal element 22 of covariance matrix
    kCutCovElement33,	// diagonal element 33 of covariance matrix
    kCutCovElement44,	// diagonal element 44 of covariance matrix
    kCutCovElement55,	// diagonal element 55 of covariance matrix
    kCutStatus,         // track status
    kNCuts=10,		// number of single selections
    kNStepQA=2,		// number of QA steps (before/after the cuts)
    kNHist=9		// number of QA histograms
  };

 private:
  void SelectionBitMap(TObject* obj);
  void DefineHistograms(); 		// books histograms and TList
  void Initialise();			// sets everything to 0
  void FillHistograms(TObject* obj, Bool_t b);
					// Fills histograms before and after cuts
  Double_t fMinNClusterTPC;		// min number of clusters in TPC
  Double_t fMinNClusterITS;		// min number of clusters in ITS
  Double_t fMaxChi2PerClusterTPC;	// max chi2 per clusters in TPC
  Double_t fMaxChi2PerClusterITS;	// max chi2 per clusters in ITS

  Double_t fCovariance11Max ;		// max covariance matrix element 11
  Double_t fCovariance22Max ;		// max covariance matrix element 22
  Double_t fCovariance33Max ;		// max covariance matrix element 33
  Double_t fCovariance44Max ;		// max covariance matrix element 44
  Double_t fCovariance55Max ;		// max covariance matrix element 55

  ULong_t fStatus;    // track status

  TH1F* fhCutStatistics;		// Histogram: statistics of what cuts the tracks did not survive
  TH2F* fhCutCorrelation;		// Histogram: 2d statistics plot

  TH1F* fhQA[kNHist][kNStepQA];		// QA Histograms
  TBits *fBitmap ; 			// stores single selection decisions

  // QA histogram setters
  Int_t fhNBinsClusterTPC;		// number of bins+1: cluster TPC
  Int_t fhNBinsClusterITS;		// number of bins+1: cluster ITS
  Int_t fhNBinsChi2TPC;			// number of bins+1: chi2 per cluster TPC
  Int_t fhNBinsChi2ITS;			// number of bins+1: chi2 per cluster ITS
  Int_t fhNBinsCovariance11;		// number of bins+1: covariance matrix element 11
  Int_t fhNBinsCovariance22;		// number of bins+1: covariance matrix element 22
  Int_t fhNBinsCovariance33;		// number of bins+1: covariance matrix element 33
  Int_t fhNBinsCovariance44;		// number of bins+1: covariance matrix element 44
  Int_t fhNBinsCovariance55;		// number of bins+1: covariance matrix element 55

  Double_t *fhBinLimClusterTPC;	//[fhNBinsClusterTPC] bin limits: cluster TPC
  Double_t *fhBinLimClusterITS;	//[fhNBinsClusterITS] bin limits: cluster ITS
  Double_t *fhBinLimChi2TPC;	//[fhNBinsChi2TPC] bin limits: chi2 per cluster TPC
  Double_t *fhBinLimChi2ITS;	//[fhNBinsChi2ITS] bin limits: chi2 per cluster ITS
  Double_t *fhBinLimCovariance11;//[fhNBinsCovariance11] bin limits: covariance matrix element 11
  Double_t *fhBinLimCovariance22;//[fhNBinsCovariance22] bin limits: covariance matrix element 22
  Double_t *fhBinLimCovariance33;//[fhNBinsCovariance33] bin limits: covariance matrix element 33
  Double_t *fhBinLimCovariance44;//[fhNBinsCovariance44] bin limits: covariance matrix element 44
  Double_t *fhBinLimCovariance55;//[fhNBinsCovariance55] bin limits: covariance matrix element 55

  ClassDef(AliCFTrackQualityCuts,3);
};

#endif
