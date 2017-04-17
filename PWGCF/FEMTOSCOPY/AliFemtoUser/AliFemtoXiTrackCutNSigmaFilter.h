///
/// \file AliFemtoXiTrackCutNSigmaFilter.h
///

#ifndef ALIFEMTOXITRACKCUTNSIGMAFILTER_H
#define ALIFEMTOXITRACKCUTNSIGMAFILTER_H

#pragma once

#include "AliFemtoTrackCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoNSigmaFilter.h"
#include "AliFemtoXiTrackCut.h"
#include "AliFemtoV0TrackCutNSigmaFilter.h"
#include "AliFemtoESDTrackCutNSigmaFilter.h"

/// \class AliFemtoXiTrackCutNSigmaFilter
/// \brief A class designed to help track cuts test particle calculated NSigma values
/// \See AliFemtoNSigmaFilter.h for detailed description of how to use filter
/// \
/// \It seems more appropriate for this class to inherit from
/// \AliFemtoV0TrackCutNSigmaFilter instead of AliFemtoXiTrackCut,
/// \but it must inherit from AliFemtoXiTrackCut in order to fit well in AliFemtoSimpleAnalysis etc.
/// \As such, this class is much longer than initially envisioned
/// \
/// \author Jesse Buxton <jesse.thomas.buxton@cern.ch>
///
class AliFemtoXiTrackCutNSigmaFilter : public AliFemtoXiTrackCut {
public:

  enum XiDaughterParticleType {
    kPion = 0,
    kKaon = 1,
    kProton = 2,
    kBacPion = 3
  };

  /// Default Constructor
  AliFemtoXiTrackCutNSigmaFilter();

  AliFemtoXiTrackCutNSigmaFilter(const AliFemtoXiTrackCutNSigmaFilter& aCut);
  AliFemtoXiTrackCutNSigmaFilter& operator=(const AliFemtoXiTrackCutNSigmaFilter& aCut);
  virtual AliFemtoXiTrackCutNSigmaFilter* Clone();
  /// Destructor - deletes filters and histograms
  virtual ~AliFemtoXiTrackCutNSigmaFilter();

  virtual bool PassV0(const AliFemtoXi* aXi);
  virtual bool Pass(const AliFemtoXi* aXi);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual TList *GetOutputList();  //include fMinvPurityAidHistoXi and fMinvPurityAidHistoV0 in the output list if they exist
  virtual AliFemtoParticleType Type(){return hbtXi;}

  void SetDaughterV0Filter(AliFemtoV0TrackCutNSigmaFilter* aCut);
  AliFemtoV0TrackCutNSigmaFilter* GetDaughterV0Filter();

  //----n sigma----
  void UpdateDaughterV0Filter();
  bool IsBacPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=false);

  void CreateCustomNSigmaFilter(XiDaughterParticleType aDaughterType);
  void CreateCustomBacPionNSigmaFilter();

  void AddTPCAndTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddBacPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);

  void AddTPCNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddBacPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);

  void AddTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddBacPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);

  //DO NOT set this to true unless you completely understand the consequences.  See AliFemtoNSigmaFilter.h for more information
  void SetOverrideImproperBacPionNSigmaFilter(bool aOverride);

  // Methods to help cut out misidentified V0s
  void CreateCustomV0Rejection(AliFemtoV0Type aV0TypeToReject);
  void AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos, double aNSigmaValueTOFPos,
                                                                  double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTPCNeg, double aNSigmaValueTOFNeg);
  void AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddTPCNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTPCPos,
                                                            double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTPCNeg);
  void AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, int aDaughterCharge, double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddTOFNSigmaCutToV0Rejection(AliFemtoV0Type aV0Type, double aMomMinPos, double aMomMaxPos, double aNSigmaValueTOFPos,
                                                            double aMomMinNeg, double aMomMaxNeg, double aNSigmaValueTOFNeg);

  double GetCTauXi(const AliFemtoXi* aXi);

 protected:
  bool fIsDaughterV0FilterUpdated;
  bool fIsDaughterBacPionFilterUpdated;
  AliFemtoV0TrackCutNSigmaFilter *fDaughterV0Filter;
  AliFemtoESDTrackCutNSigmaFilter *fDaughterBacPionFilter;


  bool fUseCustomBacPionNSigmaFilter;
  AliFemtoNSigmaFilter *fBacPionNSigmaFilter;

  TH1D* fCTauXi;  //Currently in development, if useful, will be added to AliFemtoV0TrackCut,
		  // or to the cut monitor


#ifdef __ROOT__
  ClassDef(AliFemtoXiTrackCutNSigmaFilter, 1);
#endif

};

inline void AliFemtoXiTrackCutNSigmaFilter::SetDaughterV0Filter(AliFemtoV0TrackCutNSigmaFilter* aFilter){fDaughterV0Filter = aFilter;}
inline AliFemtoV0TrackCutNSigmaFilter* AliFemtoXiTrackCutNSigmaFilter::GetDaughterV0Filter() {return fDaughterV0Filter;}

inline void AliFemtoXiTrackCutNSigmaFilter::SetOverrideImproperBacPionNSigmaFilter(bool aOverride)
{fBacPionNSigmaFilter->SetOverrideImproperConfig(aOverride);}

inline void AliFemtoXiTrackCutNSigmaFilter::CreateCustomBacPionNSigmaFilter()
{CreateCustomNSigmaFilter(kBacPion);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddBacPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{AddTPCAndTOFNSigmaCut(kBacPion, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);}


inline void AliFemtoXiTrackCutNSigmaFilter::AddBacPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{AddTPCNSigmaCut(kBacPion, aMomMin, aMomMax, aNSigmaValueTPC);}


inline void AliFemtoXiTrackCutNSigmaFilter::AddBacPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{AddTOFNSigmaCut(kBacPion, aMomMin, aMomMax, aNSigmaValueTOF);}


#endif

