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

// The class AliCFTrackKineCuts is designed to select both generated 
// and reconstructed tracks of a given range in momentum space,
// electric charge and azimuthal emission angle phi 
// and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, generated and reconstructed tracks, it is 
// applicable on ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - total momentum
// - pt
// - px
// - py
// - pz
// - eta
// - rapidity
// - phi
// - charge
// - is charged
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

#ifndef ALICFTRACKKINECUTS_H
#define ALICFTRACKKINECUTS_H

#include "AliCFCutBase.h"

class TH2;
class TBits;
class AliVParticle;

class AliCFTrackKineCuts : public AliCFCutBase
{
 public :
  AliCFTrackKineCuts() ;
  AliCFTrackKineCuts(Char_t* name, Char_t* title) ;
  AliCFTrackKineCuts(const AliCFTrackKineCuts& c) ;
  AliCFTrackKineCuts& operator=(const AliCFTrackKineCuts& c) ;
  ~AliCFTrackKineCuts();
  void Copy(TObject &c) const;

  void GetBitMap(TObject* obj, TBits *bitmap);
  Bool_t IsSelected(TObject* obj);
  void Init();

  // cut value setter
  void SetMomentumRange(Double_t momentumMin=0., Double_t momentumMax=1e99) {fMomentumMin=momentumMin; fMomentumMax=momentumMax;}
  void SetPtRange(Double_t ptMin=0., Double_t ptMax=1e99) {fPtMin=ptMin; fPtMax=ptMax;}
  void SetPxRange(Double_t pxMin=-1e99, Double_t pxMax=1e99) {fPxMin=pxMin; fPxMax=pxMax;}
  void SetPyRange(Double_t pyMin=-1e99, Double_t pyMax=1e99) {fPyMin=pyMin; fPyMax=pyMax;}
  void SetPzRange(Double_t pzMin=-1e99, Double_t pzMax=1e99) {fPzMin=pzMin; fPzMax=pzMax;}
  void SetEtaRange(Double_t etaMin=-1e99, Double_t etaMax=1e99) {fEtaMin=etaMin; fEtaMax=etaMax;}
  void SetRapidityRange(Double_t rapMin=-1e99, Double_t rapMax=1e99) {fRapidityMin=rapMin; fRapidityMax=rapMax;} 
  void SetPhiRange(Double_t phiMin=-10., Double_t phiMax=10.) {fPhiMin=phiMin; fPhiMax=phiMax;}
  void SetChargeRec(Double_t charge=10.) {fCharge=charge;}
  void SetChargeMC(Double_t charge=10.) {fCharge=charge*3.;}
  void SetRequireIsCharged(Bool_t b=kFALSE) {fRequireIsCharged=b;}

  // QA histograms
  void FillHistogramsBeforeCuts(TObject* obj) {return FillHistograms(obj,kFALSE);}
  void FillHistogramsAfterCuts(TObject* obj)  {return FillHistograms(obj,kTRUE);}
  void DrawHistograms(Bool_t drawLogScale=kTRUE);
  void SaveHistograms(const Char_t* dir = 0);
  void AddQAHistograms(TList *qaList) const;
  // QA histogram setter
  // please use indices from the enumeration below
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins);
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax);

  // indeces/counters for single selections
  enum { 
    kCutP=0,	// momentum
    kCutPt,	// pt
    kCutPx,	// px
    kCutPy,	// py
    kCutPz,	// pz
    kCutRapidity,// raptidity
    kCutEta,	// eta
    kCutPhi,	// phi
    kCutCharge,	// charge
    kNCuts=10,	// number of single selections
    kNStepQA=2,	// number of QA steps (before/after the cuts)
    kNHist=9	// number of QA histograms
  };


 private:
  TBits* SelectionBitMap(TObject* obj);
  void DefineHistograms(); 		// books histograms and TList
  void Initialise();			// sets everything to 0
  void FillHistograms(TObject* obj, Bool_t b);
					// Fills histograms before and after cuts
  Double_t fMomentumMin ;		// lower limit of accepted total momentum range
  Double_t fMomentumMax ;		// upper limit of accepted total momentum range
  Double_t fPtMin ;			// lower limit of accepted transverse momentum range
  Double_t fPtMax ;			// upper limit of accepted transverse momentum range
  Double_t fPxMin ;			// lower limit of accepted px range
  Double_t fPxMax ;			// upper limit of accepted px range
  Double_t fPyMin ;			// lower limit of accepted py range
  Double_t fPyMax ;			// upper limit of accepted py range
  Double_t fPzMin ;			// lower limit of accepted pz range
  Double_t fPzMax ;			// upper limit of accepted pz range
  Double_t fEtaMin ;			// lower limit of accepted pseudo-rapidity range
  Double_t fEtaMax ;			// upper limit of accepted pseudo-rapidity range
  Double_t fRapidityMin ;		// lower limit of accepted rapidity range
  Double_t fRapidityMax ;		// upper limit of accepted rapidity range
  Double_t fPhiMin ;			// lower limit of accepted phi range
  Double_t fPhiMax ;			// upper limit of accepted phi range
  Double_t fCharge ;			// electric charge
  Bool_t  fRequireIsCharged;		// accept charged particles only

  TH1F* fhCutStatistics;              // Histogram: statistics of what cuts the tracks did not survive
  TH2F* fhCutCorrelation;             // Histogram: 2d statistics plot

  TH1F* fhQA[kNHist][kNStepQA];		// QA Histograms
  TBits *fBitmap ; 				// stores single selection decisions

  // QA histogram setters
  Int_t fhNBinsMomentum;		// number of bins: momentum
  Int_t fhNBinsPt;				// number of bins: pt
  Int_t fhNBinsPx;				// number of bins: px
  Int_t fhNBinsPy;				// number of bins: py
  Int_t fhNBinsPz;				// number of bins: pz
  Int_t fhNBinsEta;				// number of bins: eta
  Int_t fhNBinsRapidity;			// number of bins: rapidity
  Int_t fhNBinsPhi;				// number of bins: phi
  Int_t fhNBinsCharge;			// number of bins: charge
  
  Double_t *fhBinLimMomentum;	// bin limits: momentum
  Double_t *fhBinLimPt;			// bin limits: pt
  Double_t *fhBinLimPx;		// bin limits: px
  Double_t *fhBinLimPy;		// bin limits: py
  Double_t *fhBinLimPz;		// bin limits: pz
  Double_t *fhBinLimEta;		// bin limits: eta
  Double_t *fhBinLimRapidity;	// bin limits: rapidity
  Double_t *fhBinLimPhi;		// bin limits: phi
  Double_t *fhBinLimCharge;	// bin limits: charge

  ClassDef(AliCFTrackKineCuts,1);
};

#endif
