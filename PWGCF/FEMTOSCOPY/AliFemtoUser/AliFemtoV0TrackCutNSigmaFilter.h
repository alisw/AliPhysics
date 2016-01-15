/// \class AliFemtoV0TrackCutNSigmaFilter

#ifndef ALIFEMTOV0TRACKCUTNSIGMAFILTER_H
#define ALIFEMTOV0TRACKCUTNSIGMAFILTER_H

#include "AliFemtoTrackCut.h"
#include "AliFemtoV0TrackCut.h"

#include "AliFemtoNSigmaFilter.h"
class AliFemtoNSigmaFilter;

#include "TH1D.h"

/**
 * \class AliFemtoV0TrackCutNSigmaFilter

 */
class AliFemtoV0TrackCutNSigmaFilter : public AliFemtoV0TrackCut {
public:

  AliFemtoV0TrackCutNSigmaFilter();
  virtual ~AliFemtoV0TrackCutNSigmaFilter();

  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtV0;}

  //----n sigma----
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);  //tweaked 14/12/2015
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);  //tweaked 14/12/2015
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);  //tweaked 14/12/2015


// !!!!!----- 14/12/2015 ----------------------------------------------------------
  enum DaughterParticleType {kPion=0, kKaon=1, kProton=2};

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

  //-----The MinvHisto is built immediately before the (final) Minv cut, and thus may be used to calculate the purity of the V0 collection
  // However, the use MUST manually add this histogram to the output list of the analysis.
  //   i.e. add something similar to:  tOutputList->Add(p1cut->GetMinvHisto());
  //   where p1cut is a AliFemtoV0TrackCut object
  void SetMinvHisto(const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax);   //set the Minv histogram attributes and
														//set flag fBuildMinvHisto=true
  TH1D* GetMinvHisto();

  //DO NOT set these to true unless you completely understand the consequences.  See AliFemtoNSigmaFilter.h for more information
  void SetOverrideImproperPionNSigmaFilter(bool aOverride);
  void SetOverrideImproperKaonNSigmaFilter(bool aOverride);
  void SetOverrideImproperProtonNSigmaFilter(bool aOverride);

 private:

// !!!!!----- 14/12/2015 ----------------------------------------------------------
  bool fUseCustomPionNSigmaFilter;
  bool fUseCustomKaonNSigmaFilter;
  bool fUseCustomProtonNSigmaFilter;

  AliFemtoNSigmaFilter *fPionNSigmaFilter;
  AliFemtoNSigmaFilter *fKaonNSigmaFilter;
  AliFemtoNSigmaFilter *fProtonNSigmaFilter;

  bool fBuildMinvHisto;
  TH1D* fMinvHisto;


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0TrackCutNSigmaFilter, 1);
  /// \endcond
#endif

};

inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperPionNSigmaFilter(bool aOverride) {fPionNSigmaFilter->SetOverrideImproperConfig(aOverride);}
inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperKaonNSigmaFilter(bool aOverride) {fKaonNSigmaFilter->SetOverrideImproperConfig(aOverride);}
inline void AliFemtoV0TrackCutNSigmaFilter::SetOverrideImproperProtonNSigmaFilter(bool aOverride) {fProtonNSigmaFilter->SetOverrideImproperConfig(aOverride);}

#endif
