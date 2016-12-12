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

// The class AliCFTrackIsPrimaryCut is designed to select reconstructed tracks
// with a small impact parameter and tracks which are (not) daughters of kink
// decays and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, reconstructed tracks, it is applicable on
// ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - min. and max. distance to main vertex in transverse plane (xy)
// - min. and max. longitudinal distance to main vertex (z)
// - min. and max. distance to main vertex as ellpise in xy - z plane
// - all above cuts on absolute values or in units of sigma (resolution)
// - min. and max. distance to main vertex in units of sigma (resolution)
// - max. transverse (xy) and longitudinal (z) impact parameter resolution
// - require that the dca calculation doesn't fail
// - accept or not accept daughter tracks of kink decays
//
// By default, the distance to 'vertex calculated from tracks' is used.
// Optionally the SPD (tracklet based) or TPC (TPC only tracks based) vertex
// can be used.
// Note: the distance to the TPC-vertex is already stored in the ESD,
// the distance to the SPD-vertex has to be re-calculated by propagating each
// track while executing this cut.
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

#ifndef ALICFTRACKISPRIMARYCUTS_H
#define ALICFTRACKISPRIMARYCUTS_H

#include "AliCFCutBase.h"
#include "AliAODTrack.h"
#include <TH2.h>
#include "AliESDtrackCuts.h"
class TBits;
class AliESDtrack;
class AliAODTrack;
class AliVEvent;

class AliCFTrackIsPrimaryCuts : public AliCFCutBase
{
 public :
  AliCFTrackIsPrimaryCuts() ;
  AliCFTrackIsPrimaryCuts(const Char_t* name, const Char_t* title) ;
  AliCFTrackIsPrimaryCuts(const AliCFTrackIsPrimaryCuts& c) ;
  AliCFTrackIsPrimaryCuts& operator=(const AliCFTrackIsPrimaryCuts& c) ;
  ~AliCFTrackIsPrimaryCuts();
  void Copy(TObject &c) const;

  Bool_t IsSelected(TObject* obj);
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  // cut value setter
  void UseSPDvertex(Bool_t b=kFALSE);
  void UseTPCvertex(Bool_t b=kFALSE);
  void SetMinDCAToVertexXY(Float_t dist=0.)             {fMinDCAToVertexXY = dist;}
  void SetMinDCAToVertexZ(Float_t dist=0.)              {fMinDCAToVertexZ = dist;}
  void SetMaxDCAToVertexXY(Float_t dist=1e10)           {fMaxDCAToVertexXY = dist;}
  void SetMaxDCAToVertexZ(Float_t dist=1e10)            {fMaxDCAToVertexZ = dist;}
  void SetDCAToVertex2D(Bool_t b=kFALSE)                {fDCAToVertex2D = b;}
  void SetAbsDCAToVertex(Bool_t b=kTRUE)                {fAbsDCAToVertex = b;}
  void SetMinNSigmaToVertex(Double_t sigmaMin=0.)	{fNSigmaToVertexMin = sigmaMin;}
  void SetMaxNSigmaToVertex(Double_t sigmaMax=1.e+10)	{fNSigmaToVertexMax = sigmaMax;}
  void SetMaxSigmaDCAxy(Double_t sigmaMax=1.e+10)	{fSigmaDCAxy = sigmaMax;}
  void SetMaxSigmaDCAz(Double_t sigmaMax=1.e+10)	{fSigmaDCAz = sigmaMax;}
  void SetRequireSigmaToVertex(Bool_t b=kFALSE)	        {fRequireSigmaToVertex=b;}
  void SetAcceptKinkDaughters(Bool_t b=kTRUE)	        {fAcceptKinkDaughters=b;}
  void SetAODType(Char_t type=AliAODTrack::kUndef)      {fAODType = type;}

  // QA histograms
  void DrawHistograms();
  void SaveHistograms(const Char_t* dir = 0);
  void AddQAHistograms(TList *qaList);
  // QA histogram setter
  // please use indices from the enumeration below
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins);
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax);
  virtual void SetRecEventInfo(const TObject* esd);

  // indeces/counters for single selections
  enum { 
    kCutNSigmaToVertex=0,	// tracks's distance to main vertex in units of sigma
    kCutRequireSigmaToVertex,	// calculation is successful
    kCutAcceptKinkDaughters,	// do (not) accept secondaries
    kDcaXY,			// controll histogram: dca in xy plane
    kDcaZ,			// controll histogram: dca along z axis
    kDcaXYnorm,			// controll histogram: normalised dca in xy plane
    kDcaZnorm,			// controll histogram: normalised dca along z axis
    kSigmaDcaXY,		// controll histogram: impact parameter resolution in transverse plane
    kSigmaDcaZ,			// controll histogram: impact parameter resolution along z axis
    kCutAODType,                // cut on AliAODTrack::fType
    kNCuts=9,			// number of single selections
    kNStepQA=2,			// number of QA steps (before/after the cuts)
    kNHist=9			// number of QA histograms
  };

 private:
  void SelectionBitMap(TObject* obj);
  void GetDCA(AliESDtrack* esdTrack);
  void GetDCA(AliAODTrack* aodTrack);
  void DefineHistograms(); 		// books histograms and TList
  void Initialise();			// sets everything to 0
  void FillHistograms(TObject* obj, Bool_t b);
					// Fills histograms before and after cuts
  AliVEvent*  fEvt;			// pointer to event, needed for SPD vertex and DCA in AODs
  Bool_t   fUseSPDvertex;		// flag: calculate dca to SPD-vertex, off by default
  Bool_t   fUseTPCvertex;		// flag: calculate dca to TPC-vertex, off by default
  Double_t fMinDCAToVertexXY;		// cut value: min distance to main vertex in transverse plane
  Double_t fMinDCAToVertexZ;		// cut value: min longitudinal distance to main vertex
  Double_t fMaxDCAToVertexXY;		// cut value: max distance to main vertex in transverse plane
  Double_t fMaxDCAToVertexZ;		// cut value: max longitudinal distance to main vertex
  Bool_t   fDCAToVertex2D;		// flag: cut on ellipse in xy - z plane
  Bool_t   fAbsDCAToVertex;		// flag: cut on absolute values or in units of sigma
  Double_t fNSigmaToVertexMin;		// cut value: min distance to main vertex in units of sigma
  Double_t fNSigmaToVertexMax;		// cut value: max distance to main vertex in units of sigma
  Double_t fSigmaDCAxy;			// cut value: impact parameter resolution in xy plane
  Double_t fSigmaDCAz;			// cut value: impact parameter resolution in z direction
  Double_t fDCA[6];			// impact parameters
  Bool_t   fRequireSigmaToVertex;	// require calculable distance to main vertex
  Char_t   fAODType;                    // type of AOD track (undef, primary, secondary, orphan)
                                        // applicable at AOD level only !

  TH2F* fhDcaXYvsDcaZ[2];		// Histogram: dca xy vs. z
  Bool_t fAcceptKinkDaughters;		// accepting kink daughters

  TH1F* fhCutStatistics;		// Histogram: statistics of what cuts the tracks did not survive
  TH2F* fhCutCorrelation;		// Histogram: 2d statistics plot

  TH1F* fhQA[kNHist][kNStepQA];		// QA Histograms
  TBits *fBitmap ; 			// stores single selection decisions

  // QA histogram setters
  Int_t fhNBinsNSigma;			// number of bins+1: dca in units of sigma
  Int_t fhNBinsRequireSigma;		// number of bins+1: require successful calcuation
  Int_t fhNBinsAcceptKink;		// number of bins+1: acceptkink daughters
  Int_t fhNBinsDcaXY;			// number of bins+1: dca in transverse plane
  Int_t fhNBinsDcaZ;			// number of bins+1: dca along beam axis
  Int_t fhNBinsDcaXYnorm;		// number of bins+1: normalised dca in transverse plane
  Int_t fhNBinsDcaZnorm;		// number of bins+1: normalised dca along beam axis
  Int_t fhNBinsSigmaDcaXY;		// number of bins+1: impact parameter resolution in transverse plane
  Int_t fhNBinsSigmaDcaZ;		// number of bins+1: impact parameter resolution along beam axis

  Double_t *fhBinLimNSigma; //[fhNBinsNSigma] bin limits: dca in units of sigma
  Double_t *fhBinLimRequireSigma;//[fhNBinsRequireSigma] bin limits: require successful calcuation
  Double_t *fhBinLimAcceptKink;//[fhNBinsAcceptKink] bin limits: acceptkink daughters
  Double_t *fhBinLimDcaXY;//[fhNBinsDcaXY] bin limits: dca in transverse plane
  Double_t *fhBinLimDcaZ; //[fhNBinsDcaZ] bin limits: dca along beam axis
  Double_t *fhBinLimDcaXYnorm; //[fhNBinsDcaXYnorm] bin limits: normalised dca in transverse plane
  Double_t *fhBinLimDcaZnorm;//[fhNBinsDcaZnorm] bin limits: normalised dca along beam axis
  Double_t *fhBinLimSigmaDcaXY; //[fhNBinsSigmaDcaXY] bin limits: impact parameter in transverse plane
  Double_t *fhBinLimSigmaDcaZ; //[fhNBinsSigmaDcaZ] bin limits: impact parameter along beam axis

  ClassDef(AliCFTrackIsPrimaryCuts,3);
};

#endif
