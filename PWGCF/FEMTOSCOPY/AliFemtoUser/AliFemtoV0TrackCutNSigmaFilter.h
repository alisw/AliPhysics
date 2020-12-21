///
/// \file AliFemtoV0TrackCutNSigmaFilter.h
///

#ifndef ALIFEMTOV0TRACKCUTNSIGMAFILTER_H
#define ALIFEMTOV0TRACKCUTNSIGMAFILTER_H

#pragma once

#include "AliFemtoTrackCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoNSigmaFilter.h"

class TH1D;

#include <vector>

/// \class AliFemtoV0TrackCutNSigmaFilter
/// \brief A class designed to help track cuts test particle calculated NSigma values
/// \See AliFemtoNSigmaFilter.h for detailed description of how to use filter
///
/// \author Jesse Buxton <jesse.thomas.buxton@cern.ch>
///
class AliFemtoV0TrackCutNSigmaFilter : public AliFemtoV0TrackCut {
public:

  /// Default Constructor
  AliFemtoV0TrackCutNSigmaFilter();

  AliFemtoV0TrackCutNSigmaFilter(const AliFemtoV0TrackCutNSigmaFilter& aCut);
  AliFemtoV0TrackCutNSigmaFilter& operator=(const AliFemtoV0TrackCutNSigmaFilter& aCut);
  virtual AliFemtoV0TrackCutNSigmaFilter* Clone();
  /// Destructor - deletes filters and histograms
  virtual ~AliFemtoV0TrackCutNSigmaFilter();

  /// Main method of the filter. Returns true if V0 passes NSigma
  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();  //TODO
  virtual TList* ListSettings();
  virtual TList* AppendSettings(TList *settings, const TString &prefix = "") const;

  //----n sigma----
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=false);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, double nsigmacutTPC=3.0, double nsigmacutTOF=3.0, bool requireTOF=true);

  enum DaughterParticleType {
    kPion = 0,
    kKaon = 1,
    kProton = 2
  };

  void CreateCustomNSigmaFilter(DaughterParticleType aDaughterType);
  void CreateCustomPionNSigmaFilter();
  void CreateCustomKaonNSigmaFilter();
  void CreateCustomProtonNSigmaFilter();

  void AddTPCAndTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);

  void AddTPCNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC);

  void AddTOFNSigmaCut(DaughterParticleType aDaughterType, double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);
  void AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF);

  //DO NOT set these to true unless you completely understand the consequences.  See AliFemtoNSigmaFilter.h for more information
  void SetOverrideImproperPionNSigmaFilter(bool aOverride);
  void SetOverrideImproperKaonNSigmaFilter(bool aOverride);
  void SetOverrideImproperProtonNSigmaFilter(bool aOverride);

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

  bool IsMisIDK0s(const AliFemtoV0* aV0);
  bool IsMisIDLambda(const AliFemtoV0* aV0);
  bool IsMisIDAntiLambda(const AliFemtoV0* aV0);

  void SetCTauHistoV0(int aNbins=500, double aMin=0.0, double aMax=50.0);
  short GetParticleType();
  double GetCTau(const AliFemtoV0* aV0);
  TH1D* GetCTauHist();

  virtual TList *GetOutputList();

protected:

  bool fUseCustomPionNSigmaFilter;
  bool fUseCustomKaonNSigmaFilter;
  bool fUseCustomProtonNSigmaFilter;

  AliFemtoNSigmaFilter *fPionNSigmaFilter;
  AliFemtoNSigmaFilter *fKaonNSigmaFilter;
  AliFemtoNSigmaFilter *fProtonNSigmaFilter;

  // Members to help cut out misidentified V0s
  bool fUseCustomK0sRejectionFilters;
  bool fUseCustomLambdaRejectionFilters;
  bool fUseCustomAntiLambdaRejectionFilters;

  vector<AliFemtoNSigmaFilter> fK0sRejectionFilters;
  vector<AliFemtoNSigmaFilter> fLambdaRejectionFilters;
  vector<AliFemtoNSigmaFilter> fAntiLambdaRejectionFilters;

  bool fBuildCTau;
  TH1D* fCTau;  //Currently in development, if useful, will be added to AliFemtoV0TrackCut,
		// or to the cut monitor


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0TrackCutNSigmaFilter, 1);
  /// \endcond
#endif

};



inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperPionNSigmaFilter(bool aOverride)
{
  fPionNSigmaFilter->SetOverrideImproperConfig(aOverride);
}

inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperKaonNSigmaFilter(bool aOverride)
{
  fKaonNSigmaFilter->SetOverrideImproperConfig(aOverride);
}

inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperProtonNSigmaFilter(bool aOverride)
{
  fProtonNSigmaFilter->SetOverrideImproperConfig(aOverride);
}

inline TList* AliFemtoV0TrackCutNSigmaFilter::ListSettings()
{
  return AppendSettings(new TList());
}

inline void AliFemtoV0TrackCutNSigmaFilter::CreateCustomPionNSigmaFilter() {
  CreateCustomNSigmaFilter(kPion);
}

inline void AliFemtoV0TrackCutNSigmaFilter::CreateCustomKaonNSigmaFilter()
{
  CreateCustomNSigmaFilter(kKaon);
}

inline void AliFemtoV0TrackCutNSigmaFilter::CreateCustomProtonNSigmaFilter()
{
  CreateCustomNSigmaFilter(kProton);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddPionTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  AddTPCAndTOFNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddKaonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  AddTPCAndTOFNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddProtonTPCAndTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  AddTPCAndTOFNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTPC, aNSigmaValueTOF);
}


inline void AliFemtoV0TrackCutNSigmaFilter::AddPionTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  AddTPCNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTPC);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddKaonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  AddTPCNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTPC);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddProtonTPCNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTPC)
{
  AddTPCNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTPC);
}



inline void AliFemtoV0TrackCutNSigmaFilter::AddPionTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  AddTOFNSigmaCut(kPion, aMomMin, aMomMax, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddKaonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  AddTOFNSigmaCut(kKaon, aMomMin, aMomMax, aNSigmaValueTOF);
}

inline void AliFemtoV0TrackCutNSigmaFilter::AddProtonTOFNSigmaCut(double aMomMin, double aMomMax, double aNSigmaValueTOF)
{
  AddTOFNSigmaCut(kProton, aMomMin, aMomMax, aNSigmaValueTOF);
}

inline short AliFemtoV0TrackCutNSigmaFilter::GetParticleType() {return fParticleType;}

inline TH1D* AliFemtoV0TrackCutNSigmaFilter::GetCTauHist() {return fCTau;}

#endif
