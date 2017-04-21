///
/// \file AliFemtoV0PurityBgdEstimator.h
///

#pragma once

#ifndef ALIFEMTOV0PURITYBGDESTIMATOR_H
#define ALIFEMTOV0PURITYBGDESTIMATOR_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoV0.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoV0TrackCutNSigmaFilter.h"
#include "AliFemtoV0Collection.h"
#include "AliFemtoV0SharedDaughterCut.h"

#include "AliFmPhysicalHelixD.h"
#include "TVector3.h"


/// \class AliFemtoV0PurityBgdEstimator
///
/// \brief Estimates the background to the V0 purity instead of
///        simply fitting it with a polynomial from an M_inv distribution
/// 
/// The background may be estimated by forming V0 candidates with daughters from different events
///   These V0s will pass all of the cuts imposed on the analysis V0s, but will not be real, and 
///   therefore should give a handle on the background.
///   This maybe is not the idea solution to forming mixed event V0s, but it seems that any other
///   solution would require work done much further up in the analysis tree, where I do not have
///   rights to change or add files
/// 
/// The user should create a typical femto analysis of the daughter pairs of the V0s of interest
///   My work involves (Anti)Lambda and K0Short particles, so this file will be tailored toward
///   those.  Impose the typical daughter cuts on the single particles, and AliFemtoSimpleAnalysis
///   (or similar) will form the pairs, and will fill the histograms here like any other 
///   correlation functions.
/// Note, however, that in order to be added to the histograms in this class, the V0 formed by the
///   pair must pass all of the usual cuts.  These will be included in a Pass function of the 
///   AliFemtoV0TrackCutNSigmaFilter class.  The user should create a AliFemtoV0TrackCutNSigmaFilter
///   object with all desired V0 cuts, and load it with the SetV0TrackCut method.
///   Also, the Minv cut in AliFemtoV0TrackCutNSigmaFilter should be disabled or set to loose values
/// 
/// This will not identically mimic the creation of V0s.  This is done much further up the analysis
///   tree.  However, this should do a good job.  Some cuts simply cannot be included.
///   If everything works as planned, fNumerator should match closely the invariant mass plot (built
///   last, directly before final invariant mass cut) obtained from the fMinvPurityAidHistoV0 
///   histogram in AliFemtoV0TrackCut, and will show the signal+background in the region  of
///   interest.  fDenominator will show just the background, and this may be used in the purity
///   calculation in place of simply fitting the background with a polynomial.
/// 
/// I have obtained the best results when using FilterBit 0 (instead of typical value of 7) and using
///   DCA from global track instead of TPC only, i.e. in the event reader, call
///   AliFemtoEventReaderAOD::SetFilterBit(0) and
///   AliFemtoEventReaderAOD::SetDCAglobalTrack(1)
/// 
///   NOTE:  This code is very much still in development.  This note will be excluded when the
///     author is comfortable with others using the code.  Of course, if you are brave, you can
///     still test things out before this note is gone.
///     Suggestions are welcome.
/// 
/// \author  Jesse Buxton, Ohio State University, <jesse.thomas.buxton@cern.ch>
///            
///
class AliFemtoV0PurityBgdEstimator : public AliFemtoCorrFctn {
public:

  /// Default Constructor
  AliFemtoV0PurityBgdEstimator();

  /// Constructor
  AliFemtoV0PurityBgdEstimator(const char* title,
                               const int nbins,
                               const float MinvLo,
                               const float MinvHi);

  AliFemtoV0PurityBgdEstimator(const AliFemtoV0PurityBgdEstimator& aCorrFctn);
  AliFemtoV0PurityBgdEstimator& operator=(const AliFemtoV0PurityBgdEstimator& aCorrFctn);
  virtual AliFemtoV0PurityBgdEstimator* Clone();

  virtual ~AliFemtoV0PurityBgdEstimator();

  void ClearV0Collections();
  void FillHistograms();
  void IsSameEvent();  //not great, but should work for now

  AliFmThreeVector<double> ShiftMomentumToDCA(const AliFmHelix& tHelix, const AliFmThreeVector<double>& tMomentum, double tPathLengthToDCA);
  void UseCorrectedDaughterHelices(const AliFemtoTrack* tTrackPos, const AliFemtoTrack* tTrackNeg);
  void BuildV0(AliFemtoPair* aPair);

  virtual AliFemtoString Report();
  virtual TList* GetOutputList();
  virtual void Finish();
  void Write();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  //inline functions
  TH1D* Numerator();
  TH1D* Denominator();
  TH1D* Ratio();

  void SetV0TrackCut(AliFemtoV0TrackCutNSigmaFilter* aCut);

protected:
  TString fTitle;
  int fNbinsMinv;
  double fMinvLow, fMinvHigh;
  TH1D* fNumeratorWithoutSharedDaughterCut;          // numerator - real pairs
  TH1D* fDenominatorWithoutSharedDaughterCut;        // denominator - mixed pairs
  TH1D* fRatioWithoutSharedDaughterCut;              // unnormalized ratio num/den
  TH1D* fNumerator;          // numerator - real pairs
  TH1D* fDenominator;        // denominator - mixed pairs
  TH1D* fRatio;              // unnormalized ratio num/den

  AliFemtoV0 *fFemtoV0;
  AliFemtoV0TrackCutNSigmaFilter* fFemtoV0TrackCut;

  AliFmThreeVector<double> fCurrentPrimVtx;
  AliFemtoV0Collection* fRealV0Collection;
  AliFemtoV0Collection* fMixedV0Collection;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0PurityBgdEstimator, 1);
  /// \endcond
#endif

};

inline TH1D* AliFemtoV0PurityBgdEstimator::Numerator(){return fNumerator;}
inline TH1D* AliFemtoV0PurityBgdEstimator::Denominator(){return fDenominator;}
inline TH1D* AliFemtoV0PurityBgdEstimator::Ratio(){return fRatio;}

inline void AliFemtoV0PurityBgdEstimator::SetV0TrackCut(AliFemtoV0TrackCutNSigmaFilter* aCut){fFemtoV0TrackCut = aCut;}


#endif
