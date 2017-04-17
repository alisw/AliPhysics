///
/// \file AliFemtoESDTrackCutNSigmaFilter.h
///

#ifndef ALIFEMTOESDTRACKCUTNSIGMAFILTER_H
#define ALIFEMTOESDTRACKCUTNSIGMAFILTER_H

#include "AliFemtoESDTrackCut.h"
#include "AliFemtoNSigmaFilter.h"

/// \class AliFemtoESDTrackCutNSigmaFilter
/// \brief A class designed to help track cuts test particle calculated NSigma values
/// \See AliFemtoNSigmaFilter.h for detailed description of how to use filter
///
/// \author Jesse Buxton <jesse.thomas.buxton@cern.ch>
///
class AliFemtoESDTrackCutNSigmaFilter : public AliFemtoESDTrackCut {

public:

  AliFemtoESDTrackCutNSigmaFilter();
  AliFemtoESDTrackCutNSigmaFilter(const AliFemtoESDTrackCutNSigmaFilter& aCut);
  AliFemtoESDTrackCutNSigmaFilter& operator=(const AliFemtoESDTrackCutNSigmaFilter& aCut);  //TODO
  virtual AliFemtoESDTrackCutNSigmaFilter* Clone();
  virtual ~AliFemtoESDTrackCutNSigmaFilter();  //TODO

  virtual bool Pass(const AliFemtoTrack* aTrack);

//  virtual AliFemtoString Report();  //TODO
//  virtual TList* ListSettings();  //TODO
//  virtual TList* AppendSettings(TList *settings, const TString &prefix = "") const;  //TODO



  enum ParticleType {
    kPion = 0,
    kKaon = 1,
    kProton = 2,
    kElectron = 3
  };

  //Used in IsKaonNSigma, IsPionNSigma, and IsProtonNSigma methods
  void CreateCustomNSigmaFilter(ParticleType aType);
  void AddTPCAndTOFNSigmaCut(ParticleType aType, double aMomMin, double aMomMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddTPCNSigmaCut(ParticleType aType, double aMomMin, double aMomMax, double aNSigmaValueTPC);
  void AddTOFNSigmaCut(ParticleType aType, double aMomMin, double aMomMax, double aNSigmaValueTOF);

  //DO NOT set these to true unless you completely understand the consequences.  See AliFemtoNSigmaFilter.h for more information
  void SetOverrideImproperPionNSigmaFilter(ParticleType aType, bool aOverride);

  void SetPionRejection(bool aReject);
  void SetUseCustomElectronRejection(bool aReject);

private:
  //Used in IsKaonNSigma, IsPionNSigma, and IsProtonNSigma methods
  bool fUseCustomPionNSigmaFilter;
  bool fUseCustomKaonNSigmaFilter;
  bool fUseCustomProtonNSigmaFilter;
  bool fUseCustomElectronNSigmaFilter;

  AliFemtoNSigmaFilter *fPionNSigmaFilter;
  AliFemtoNSigmaFilter *fKaonNSigmaFilter;
  AliFemtoNSigmaFilter *fProtonNSigmaFilter;
  AliFemtoNSigmaFilter *fElectronNSigmaFilter;

  bool fPionRejection;  // 26/02/2016


  //----n sigma----
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);
  bool IsElectronNSigma(float mom, float nsigmaTPCE, float nsigmaTOFE);
  bool IsProbableElectron(const AliFemtoTrack* aTrack);

#ifdef __ROOT__
  ClassDef(AliFemtoESDTrackCutNSigmaFilter, 1)
#endif

};

inline void AliFemtoESDTrackCutNSigmaFilter::SetPionRejection(bool aReject) {fPionRejection = aReject;}
inline void AliFemtoESDTrackCutNSigmaFilter::SetUseCustomElectronRejection(bool aReject) {fUseCustomElectronNSigmaFilter = aReject;}

#endif
