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

  //----n sigma----
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=false);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=true);
  bool IsBacPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=false);

  void CreateCustomNSigmaFilter(XiDaughterParticleType aDaughterType);
  void CreateCustomPionNSigmaFilter();
  void CreateCustomKaonNSigmaFilter();
  void CreateCustomProtonNSigmaFilter();
  void CreateCustomBacPionNSigmaFilter();

  void AddTPCAndTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddBacPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);

  void AddTPCNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddBacPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);

  void AddTOFNSigmaCut(XiDaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddBacPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);

  //DO NOT set this to true unless you completely understand the consequences.  See AliFemtoNSigmaFilter.h for more information
  void SetOverrideImproperPionNSigmaFilter(bool aOverride);
  void SetOverrideImproperKaonNSigmaFilter(bool aOverride);
  void SetOverrideImproperProtonNSigmaFilter(bool aOverride);
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

  //-Somewhat misleading, the following are checking to see if V0 daughter is misidentified
  bool IsMisIDK0s(const AliFemtoXi* aXi);
  bool IsMisIDLambda(const AliFemtoXi* aXi);
  bool IsMisIDAntiLambda(const AliFemtoXi* aXi);

 protected:
  bool fUseCustomPionNSigmaFilter;
  bool fUseCustomKaonNSigmaFilter;
  bool fUseCustomProtonNSigmaFilter;
  bool fUseCustomBacPionNSigmaFilter;

  AliFemtoNSigmaFilter *fPionNSigmaFilter;
  AliFemtoNSigmaFilter *fKaonNSigmaFilter;
  AliFemtoNSigmaFilter *fProtonNSigmaFilter;
  AliFemtoNSigmaFilter *fBacPionNSigmaFilter;

  // Members to help cut out misidentified V0s
  bool fUseCustomK0sRejectionFilters;
  bool fUseCustomLambdaRejectionFilters;
  bool fUseCustomAntiLambdaRejectionFilters;

  vector<AliFemtoNSigmaFilter> fK0sRejectionFilters;
  vector<AliFemtoNSigmaFilter> fLambdaRejectionFilters;
  vector<AliFemtoNSigmaFilter> fAntiLambdaRejectionFilters;


#ifdef __ROOT__
  ClassDef(AliFemtoXiTrackCutNSigmaFilter, 1);
#endif

};


inline void AliFemtoXiTrackCutNSigmaFilter::SetOverrideImproperPionNSigmaFilter(bool aOverride)
{fPionNSigmaFilter->SetOverrideImproperConfig(aOverride);}

inline void AliFemtoXiTrackCutNSigmaFilter::SetOverrideImproperKaonNSigmaFilter(bool aOverride)
{fKaonNSigmaFilter->SetOverrideImproperConfig(aOverride);}

inline void AliFemtoXiTrackCutNSigmaFilter::SetOverrideImproperProtonNSigmaFilter(bool aOverride)
{fProtonNSigmaFilter->SetOverrideImproperConfig(aOverride);}

inline void AliFemtoXiTrackCutNSigmaFilter::SetOverrideImproperBacPionNSigmaFilter(bool aOverride)
{fBacPionNSigmaFilter->SetOverrideImproperConfig(aOverride);}



inline void AliFemtoXiTrackCutNSigmaFilter::CreateCustomPionNSigmaFilter()
{CreateCustomNSigmaFilter(kPion);}

inline void AliFemtoXiTrackCutNSigmaFilter::CreateCustomKaonNSigmaFilter()
{CreateCustomNSigmaFilter(kKaon);}

inline void AliFemtoXiTrackCutNSigmaFilter::CreateCustomProtonNSigmaFilter()
{CreateCustomNSigmaFilter(kProton);}

inline void AliFemtoXiTrackCutNSigmaFilter::CreateCustomBacPionNSigmaFilter()
{CreateCustomNSigmaFilter(kBacPion);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{AddTPCAndTOFNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{AddTPCAndTOFNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{AddTPCAndTOFNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddBacPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{AddTPCAndTOFNSigmaCut(kBacPion, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);}



inline void AliFemtoXiTrackCutNSigmaFilter::AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{AddTPCNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTPC);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{AddTPCNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTPC);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{AddTPCNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTPC);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddBacPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{AddTPCNSigmaCut(kBacPion, aMomMin, aMomMax, aNSigmaValueTPC);}



inline void AliFemtoXiTrackCutNSigmaFilter::AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{AddTOFNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTOF);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{AddTOFNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTOF);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{AddTOFNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTOF);}

inline void AliFemtoXiTrackCutNSigmaFilter::AddBacPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{AddTOFNSigmaCut(kBacPion, aMomMin, aMomMax, aNSigmaValueTOF);}


#endif

